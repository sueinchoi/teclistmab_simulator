#===============================================================================
# Teclistamab PK-Response Association Analysis
#
# For patients with ≥4 weeks of treatment:
# 1. VGPR or better (≥VGPR) - Odds Ratio (Logistic Regression)
# 2. 2-month PFS status - Odds Ratio (Logistic Regression)
# 3. PFS - Hazard Ratio (Cox Regression)
#===============================================================================

library(mrgsolve)
library(tidyverse)
library(survival)
library(broom)
library(gridExtra)

set.seed(12345)

cat("==========================================================\n")
cat("    Teclistamab PK-Response Association Analysis\n")
cat("==========================================================\n\n")

#-------------------------------------------------------------------------------
# 1. mrgsolve Model Definition
#-------------------------------------------------------------------------------

model_code <- '
$PARAM @annotated
CL1  : 0.449   : Linear clearance (L/day)
CL2  : 0.547   : Time-dependent clearance component (L/day)
KDES : 0.0292  : Clearance decay rate constant (1/day)
V1   : 4.13    : Central volume (L)
V2   : 1.34    : Peripheral volume (L)
Q    : 0.039   : Intercompartmental clearance (L/day)
KA   : 0.133   : Absorption rate constant (1/day)
F1   : 0.718   : Bioavailability

$CMT @annotated
DEPOT  : SC depot compartment (mg)
CENT   : Central compartment (mg)
PERIPH : Peripheral compartment (mg)

$GLOBAL
#define CP (CENT/V1)

$ODE
double TDAY = SOLVERTIME / 24.0;
double CL = CL1 + CL2 * exp(-KDES * TDAY);
double KA_HR  = KA / 24.0;
double K10_HR = (CL / V1) / 24.0;
double K12_HR = (Q / V1) / 24.0;
double K21_HR = (Q / V2) / 24.0;

dxdt_DEPOT  = -KA_HR * DEPOT;
dxdt_CENT   = KA_HR * DEPOT * F1 - K10_HR * CENT - K12_HR * CENT + K21_HR * PERIPH;
dxdt_PERIPH = K12_HR * CENT - K21_HR * PERIPH;

$TABLE
double DV = CP;

$CAPTURE @annotated
DV : Plasma concentration (mg/L)
'

cat("Compiling mrgsolve model...\n")
mod <- mcode("teclistamab_response", model_code, compile = TRUE)

#-------------------------------------------------------------------------------
# 2. IIV Parameters
#-------------------------------------------------------------------------------

cv_to_omega <- function(cv) sqrt(log(1 + cv^2))

IIV_OMEGA <- list(
  CL1 = cv_to_omega(0.536),
  CL2 = cv_to_omega(1.07),
  V1  = cv_to_omega(0.488),
  KA  = cv_to_omega(0.452)
)

#-------------------------------------------------------------------------------
# 3. Load Data and Filter for ≥4 Weeks Treatment
#-------------------------------------------------------------------------------

cat("\nLoading data...\n")

# Load patient parameters
params <- read_csv("mrgsolve_params_full.csv", show_col_types = FALSE)
cat("  - Loaded", nrow(params), "patients\n")

# Load dosing data
dosing_all <- read_csv("mrgsolve_dosing_full.csv", show_col_types = FALSE) %>%
  filter(TIME >= 0) %>%
  arrange(ID, TIME)

# Calculate treatment duration per patient (max TIME in hours)
treatment_duration <- dosing_all %>%
  group_by(ID) %>%
  summarise(
    N_doses = n(),
    Max_TIME_hr = max(TIME),
    Treatment_days = max(TIME) / 24,
    .groups = "drop"
  )

cat("\n=== Treatment Duration Summary ===\n")
print(treatment_duration)

# Filter for ≥4 weeks treatment (28 days = 672 hours)
patients_4weeks <- treatment_duration %>%
  filter(Treatment_days >= 28)

cat("\n>>> Patients with ≥4 weeks treatment:", nrow(patients_4weeks), "out of", nrow(treatment_duration), "<<<\n")

if (nrow(patients_4weeks) == 0) {
  stop("No patients with ≥4 weeks of treatment found!")
}

# Get filtered dosing data
dosing_filtered <- dosing_all %>%
  filter(ID %in% patients_4weeks$ID)

# Get filtered params
params_filtered <- params %>%
  filter(ID %in% patients_4weeks$ID)

cat("\nFiltered patients IDs:", paste(patients_4weeks$ID, collapse = ", "), "\n")

# Load response data
response_data <- read_csv("response_data.csv", show_col_types = FALSE)
cat("  - Loaded response data for", nrow(response_data), "patients\n")

#-------------------------------------------------------------------------------
# 4. Monte Carlo Simulation for PK Metrics
#-------------------------------------------------------------------------------

cat("\n=== Running Monte Carlo Simulations (1000 per patient) ===\n")

N_MC <- 1000

# Function to calculate Cavg using trapezoidal rule
calc_cavg <- function(time, conc) {
  if (length(time) < 2) return(NA)
  auc <- sum(diff(time) * (head(conc, -1) + tail(conc, -1)) / 2)
  auc / (max(time) - min(time))
}

# Run simulations for each filtered patient
all_pk_results <- list()

for (i in seq_len(nrow(params_filtered))) {
  pt <- params_filtered[i, ]
  pt_id <- pt$ID

  cat(sprintf("  Patient %d (ID=%d)...\n", i, pt_id))

  # Get individual dosing
  pt_dosing <- dosing_filtered %>%
    filter(ID == pt_id) %>%
    select(TIME, AMT, CMT, EVID) %>%
    mutate(ID = 1)  # Reset to ID=1 for simulation

  if (nrow(pt_dosing) == 0) {
    cat("    -> No dosing data, skipping\n")
    next
  }

  # Find dose times for interval calculation
  dose_times <- sort(unique(pt_dosing$TIME))

  # Run MC simulations
  mc_results <- list()

  for (sim in 1:N_MC) {
    # Sample IIV
    eta_CL1 <- rnorm(1, 0, IIV_OMEGA$CL1)
    eta_CL2 <- rnorm(1, 0, IIV_OMEGA$CL2)
    eta_V1  <- rnorm(1, 0, IIV_OMEGA$V1)
    eta_KA  <- rnorm(1, 0, IIV_OMEGA$KA)

    # Apply IIV to parameters
    ind_params <- list(
      CL1 = pt$CL1 * exp(eta_CL1),
      CL2 = pt$CL2 * exp(eta_CL2),
      V1  = pt$V1 * exp(eta_V1),
      KA  = pt$KA * exp(eta_KA)
    )

    # Simulate
    max_time <- max(pt_dosing$TIME) + 168

    sim_out <- mod %>%
      data_set(as.data.frame(pt_dosing)) %>%
      param(ind_params) %>%
      mrgsim(end = max_time, delta = 1) %>%
      as_tibble()

    # Calculate PK metrics
    # 72hr and 120hr metrics
    d72 <- sim_out %>% filter(time <= 72)
    d120 <- sim_out %>% filter(time <= 120)

    Cmax_72hr <- if (nrow(d72) > 0) max(d72$DV, na.rm = TRUE) else NA
    Cavg_72hr <- if (nrow(d72) > 0) calc_cavg(d72$time, d72$DV) else NA
    Cmax_120hr <- if (nrow(d120) > 0) max(d120$DV, na.rm = TRUE) else NA
    Cavg_120hr <- if (nrow(d120) > 0) calc_cavg(d120$time, d120$DV) else NA

    # 1st dosing interval (dose 1 to dose 2)
    if (length(dose_times) >= 2) {
      d1_start <- dose_times[1]
      d1_end <- dose_times[2]
      d_dose1 <- sim_out %>% filter(time >= d1_start, time < d1_end)
      Cmax_dose1 <- if (nrow(d_dose1) > 0) max(d_dose1$DV, na.rm = TRUE) else NA
      Cavg_dose1 <- if (nrow(d_dose1) > 0) calc_cavg(d_dose1$time, d_dose1$DV) else NA
    } else {
      Cmax_dose1 <- Cavg_dose1 <- NA
    }

    # 3rd dosing interval (dose 3 to dose 4)
    if (length(dose_times) >= 4) {
      d3_start <- dose_times[3]
      d3_end <- dose_times[4]
      d_dose3 <- sim_out %>% filter(time >= d3_start, time < d3_end)
      Cmax_dose3 <- if (nrow(d_dose3) > 0) max(d_dose3$DV, na.rm = TRUE) else NA
      Cavg_dose3 <- if (nrow(d_dose3) > 0) calc_cavg(d_dose3$time, d_dose3$DV) else NA
    } else {
      Cmax_dose3 <- Cavg_dose3 <- NA
    }

    mc_results[[sim]] <- tibble(
      Cmax_72hr = Cmax_72hr,
      Cavg_72hr = Cavg_72hr,
      Cmax_120hr = Cmax_120hr,
      Cavg_120hr = Cavg_120hr,
      Cmax_dose1 = Cmax_dose1,
      Cavg_dose1 = Cavg_dose1,
      Cmax_dose3 = Cmax_dose3,
      Cavg_dose3 = Cavg_dose3
    )
  }

  # Summarize MC results (median)
  mc_df <- bind_rows(mc_results)

  pk_summary <- tibble(
    ID = pt_id,
    PID = pt$PID,
    Cmax_72hr = median(mc_df$Cmax_72hr, na.rm = TRUE),
    Cavg_72hr = median(mc_df$Cavg_72hr, na.rm = TRUE),
    Cmax_120hr = median(mc_df$Cmax_120hr, na.rm = TRUE),
    Cavg_120hr = median(mc_df$Cavg_120hr, na.rm = TRUE),
    Cmax_dose1 = median(mc_df$Cmax_dose1, na.rm = TRUE),
    Cavg_dose1 = median(mc_df$Cavg_dose1, na.rm = TRUE),
    Cmax_dose3 = median(mc_df$Cmax_dose3, na.rm = TRUE),
    Cavg_dose3 = median(mc_df$Cavg_dose3, na.rm = TRUE)
  )

  all_pk_results[[i]] <- pk_summary
}

pk_data <- bind_rows(all_pk_results)
cat("\n  -> PK simulation complete for", nrow(pk_data), "patients\n")

#-------------------------------------------------------------------------------
# 5. Merge PK Data with Response Data and Create Analysis Dataset
#-------------------------------------------------------------------------------

cat("\n=== Merging data ===\n")

# Merge PK data with response data
filtered_data <- pk_data %>%
  left_join(response_data, by = "PID") %>%
  left_join(patients_4weeks %>% select(ID, N_doses, Treatment_days), by = "ID") %>%
  mutate(
    # VGPR or better (sCR, CR, VGPR)
    VGPR_or_better = ifelse(RESP_CTX %in% c("sCR", "CR", "VGPR"), 1, 0),

    # 2-month (60 days) PFS status
    PFS_2month = case_when(
      DAYS_VS_PFS_CTX >= 60 ~ 1,
      VS_PFS_CTX == 0 ~ 1,
      VS_PFS_CTX == 1 & DAYS_VS_PFS_CTX < 60 ~ 0,
      TRUE ~ NA_real_
    ),

    # For Cox model
    PFS_event = VS_PFS_CTX,
    PFS_time = DAYS_VS_PFS_CTX
  )

# Summary of response outcomes
cat("\n=== Response Summary (patients with ≥4 weeks treatment) ===\n")
cat("Total patients:", nrow(filtered_data), "\n")
cat("\nResponse distribution:\n")
print(table(filtered_data$RESP_CTX, useNA = "ifany"))

cat("\nVGPR or better:\n")
print(table(filtered_data$VGPR_or_better, useNA = "ifany"))

cat("\n2-month PFS status:\n")
print(table(filtered_data$PFS_2month, useNA = "ifany"))

cat("\nPFS events:\n")
print(table(filtered_data$PFS_event, useNA = "ifany"))

# Save filtered data
write_csv(filtered_data, "pk_response_analysis_data.csv")
cat("\nSaved filtered analysis data to: pk_response_analysis_data.csv\n")

#-------------------------------------------------------------------------------
# 6. Define PK Metrics for Analysis
#-------------------------------------------------------------------------------

pk_metrics <- c("Cmax_72hr", "Cavg_72hr", "Cmax_120hr", "Cavg_120hr",
                "Cmax_dose1", "Cavg_dose1", "Cmax_dose3", "Cavg_dose3")

#-------------------------------------------------------------------------------
# 7. Logistic Regression Function (for OR)
#-------------------------------------------------------------------------------

run_logistic_analysis <- function(data, outcome_var, pk_var, outcome_name) {

  # Filter valid data
  model_data <- data %>%
    filter(!is.na(.data[[outcome_var]]) & !is.na(.data[[pk_var]]))

  if (nrow(model_data) < 5) return(NULL)
  if (length(unique(model_data[[outcome_var]])) < 2) return(NULL)

  # Standardize PK variable for comparable ORs
  model_data$pk_std <- scale(model_data[[pk_var]])[,1]

  # Fit logistic regression
  formula <- as.formula(paste(outcome_var, "~ pk_std"))

  tryCatch({
    fit <- glm(formula, data = model_data, family = binomial)

    # Extract results
    coef_summary <- summary(fit)$coefficients

    # Calculate OR and 95% CI
    beta <- coef_summary["pk_std", "Estimate"]
    se <- coef_summary["pk_std", "Std. Error"]
    p_value <- coef_summary["pk_std", "Pr(>|z|)"]

    or <- exp(beta)
    or_lower <- exp(beta - 1.96 * se)
    or_upper <- exp(beta + 1.96 * se)

    # Also calculate OR per unit (not standardized)
    pk_sd <- sd(model_data[[pk_var]], na.rm = TRUE)
    or_per_unit <- exp(beta / pk_sd * 0.1)  # OR per 0.1 unit increase

    tibble(
      Outcome = outcome_name,
      `PK Metric` = pk_var,
      N = nrow(model_data),
      `N Events` = sum(model_data[[outcome_var]] == 1),
      `OR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", or, or_lower, or_upper),
      OR = or,
      OR_lower = or_lower,
      OR_upper = or_upper,
      `p-value` = p_value,
      Sig = ifelse(p_value < 0.05, "*", "")
    )
  }, error = function(e) {
    NULL
  })
}

#-------------------------------------------------------------------------------
# 8. Cox Regression Function (for HR)
#-------------------------------------------------------------------------------

run_cox_analysis <- function(data, pk_var) {

  # Filter valid data
  model_data <- data %>%
    filter(!is.na(PFS_event) & !is.na(PFS_time) & !is.na(.data[[pk_var]])) %>%
    filter(PFS_time > 0)

  if (nrow(model_data) < 5) return(NULL)
  if (sum(model_data$PFS_event) < 2) return(NULL)

  # Standardize PK variable
  model_data$pk_std <- scale(model_data[[pk_var]])[,1]

  tryCatch({
    # Fit Cox model
    fit <- coxph(Surv(PFS_time, PFS_event) ~ pk_std, data = model_data)

    # Extract results
    fit_summary <- summary(fit)

    hr <- fit_summary$conf.int[1, "exp(coef)"]
    hr_lower <- fit_summary$conf.int[1, "lower .95"]
    hr_upper <- fit_summary$conf.int[1, "upper .95"]
    p_value <- fit_summary$coefficients[1, "Pr(>|z|)"]

    tibble(
      Outcome = "PFS",
      `PK Metric` = pk_var,
      N = nrow(model_data),
      `N Events` = sum(model_data$PFS_event),
      `HR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", hr, hr_lower, hr_upper),
      HR = hr,
      HR_lower = hr_lower,
      HR_upper = hr_upper,
      `p-value` = p_value,
      Sig = ifelse(p_value < 0.05, "*", "")
    )
  }, error = function(e) {
    NULL
  })
}

#-------------------------------------------------------------------------------
# 9. Run Analyses
#-------------------------------------------------------------------------------

cat("\n")
cat("==========================================================\n")
cat("             STATISTICAL ANALYSIS RESULTS\n")
cat("==========================================================\n")

# 9.1 VGPR or better (Logistic Regression - OR)
cat("\n\n--- 1. VGPR or Better (≥VGPR) - Odds Ratio ---\n")

vgpr_results <- map_dfr(pk_metrics, function(pk_var) {
  run_logistic_analysis(filtered_data, "VGPR_or_better", pk_var, "≥VGPR")
})

if (nrow(vgpr_results) > 0) {
  print(vgpr_results %>% select(-OR, -OR_lower, -OR_upper), n = 100)
} else {
  cat("Insufficient data for VGPR analysis\n")
}

# 9.2 2-month PFS (Logistic Regression - OR)
cat("\n\n--- 2. 2-Month PFS Status - Odds Ratio ---\n")

pfs2m_results <- map_dfr(pk_metrics, function(pk_var) {
  run_logistic_analysis(filtered_data, "PFS_2month", pk_var, "2-month PFS")
})

if (nrow(pfs2m_results) > 0) {
  print(pfs2m_results %>% select(-OR, -OR_lower, -OR_upper), n = 100)
} else {
  cat("Insufficient data for 2-month PFS analysis\n")
}

# 9.3 PFS (Cox Regression - HR)
cat("\n\n--- 3. Progression-Free Survival (PFS) - Hazard Ratio ---\n")

pfs_results <- map_dfr(pk_metrics, function(pk_var) {
  run_cox_analysis(filtered_data, pk_var)
})

if (nrow(pfs_results) > 0) {
  print(pfs_results %>% select(-HR, -HR_lower, -HR_upper), n = 100)
} else {
  cat("Insufficient data for PFS analysis\n")
}

#-------------------------------------------------------------------------------
# 10. Combine and Save Results
#-------------------------------------------------------------------------------

# Combine all results
all_results <- bind_rows(
  vgpr_results %>% mutate(Analysis = "Logistic (OR)"),
  pfs2m_results %>% mutate(Analysis = "Logistic (OR)"),
  pfs_results %>%
    rename(`OR (95% CI)` = `HR (95% CI)`, OR = HR, OR_lower = HR_lower, OR_upper = HR_upper) %>%
    mutate(Analysis = "Cox (HR)")
)

write_csv(all_results, "pk_response_statistical_results.csv")
cat("\n\nSaved all results to: pk_response_statistical_results.csv\n")

#-------------------------------------------------------------------------------
# 11. Create Forest Plots
#-------------------------------------------------------------------------------

cat("\n=== Creating Forest Plots ===\n")

create_forest_plot <- function(results, title, estimate_type = "OR") {

  if (nrow(results) == 0) return(NULL)

  plot_data <- results %>%
    mutate(
      estimate = if (estimate_type == "OR") OR else HR,
      lower = if (estimate_type == "OR") OR_lower else HR_lower,
      upper = if (estimate_type == "OR") OR_upper else HR_upper,
      label = sprintf("%.2f (%.2f-%.2f)", estimate, lower, upper),
      pk_metric = factor(`PK Metric`, levels = rev(pk_metrics))
    )

  ggplot(plot_data, aes(x = estimate, y = pk_metric)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, color = "#3498db") +
    geom_point(size = 3, color = "#e74c3c") +
    geom_text(aes(label = label, x = max(upper) * 1.5), hjust = 0, size = 3) +
    geom_text(aes(label = ifelse(Sig == "*", "*", ""), x = estimate),
              vjust = -1, size = 5, color = "red") +
    scale_x_log10() +
    labs(
      title = title,
      x = estimate_type,
      y = "PK Metric"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.text.y = element_text(size = 10)
    )
}

# Forest plot for VGPR
if (nrow(vgpr_results) > 0) {
  p_vgpr <- create_forest_plot(vgpr_results, "VGPR or Better - Odds Ratio", "OR")
  if (!is.null(p_vgpr)) {
    ggsave("forest_VGPR_OR.png", p_vgpr, width = 10, height = 6, dpi = 200)
    cat("Saved: forest_VGPR_OR.png\n")
  }
}

# Forest plot for 2-month PFS
if (nrow(pfs2m_results) > 0) {
  p_pfs2m <- create_forest_plot(pfs2m_results, "2-Month PFS - Odds Ratio", "OR")
  if (!is.null(p_pfs2m)) {
    ggsave("forest_2monthPFS_OR.png", p_pfs2m, width = 10, height = 6, dpi = 200)
    cat("Saved: forest_2monthPFS_OR.png\n")
  }
}

# Forest plot for PFS HR
if (nrow(pfs_results) > 0) {
  pfs_plot_data <- pfs_results %>%
    rename(OR = HR, OR_lower = HR_lower, OR_upper = HR_upper)
  p_pfs <- create_forest_plot(pfs_plot_data, "Progression-Free Survival - Hazard Ratio", "HR")
  if (!is.null(p_pfs)) {
    ggsave("forest_PFS_HR.png", p_pfs, width = 10, height = 6, dpi = 200)
    cat("Saved: forest_PFS_HR.png\n")
  }
}

#-------------------------------------------------------------------------------
# 12. Summary Boxplots by Response
#-------------------------------------------------------------------------------

cat("\nCreating boxplots by response...\n")

# Boxplot for VGPR
create_response_boxplot <- function(data, response_var, response_label, pk_var) {
  plot_data <- data %>%
    filter(!is.na(.data[[response_var]]) & !is.na(.data[[pk_var]])) %>%
    mutate(Response = factor(ifelse(.data[[response_var]] == 1, "Yes", "No"),
                              levels = c("No", "Yes")))

  if (nrow(plot_data) < 4) return(NULL)

  ggplot(plot_data, aes(x = Response, y = .data[[pk_var]], fill = Response)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.15, alpha = 0.6, size = 2.5) +
    scale_fill_manual(values = c("No" = "#e74c3c", "Yes" = "#27ae60")) +
    labs(title = response_label, x = "", y = pk_var) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(size = 10, face = "bold"))
}

# Key metrics
key_metrics <- c("Cavg_72hr", "Cavg_120hr", "Cmax_72hr", "Cmax_120hr")

# VGPR boxplots
vgpr_plots <- list()
for (pk_var in key_metrics) {
  p <- create_response_boxplot(filtered_data, "VGPR_or_better", "≥VGPR", pk_var)
  if (!is.null(p)) {
    p <- p + labs(title = pk_var, x = "≥VGPR")
    vgpr_plots[[pk_var]] <- p
  }
}

if (length(vgpr_plots) > 0) {
  combined <- arrangeGrob(grobs = vgpr_plots, ncol = 2, nrow = 2,
                           top = "PK Metrics by VGPR Response")
  ggsave("boxplot_VGPR_response.png", combined, width = 10, height = 10, dpi = 200)
  cat("Saved: boxplot_VGPR_response.png\n")
}

#-------------------------------------------------------------------------------
# 13. Final Summary
#-------------------------------------------------------------------------------

cat("\n")
cat("==========================================================\n")
cat("                    FINAL SUMMARY\n")
cat("==========================================================\n\n")

cat("Analysis Population:\n")
cat("  - Patients with ≥4 weeks treatment:", nrow(filtered_data), "\n")
cat("  - Patients with response data:", sum(!is.na(filtered_data$RESP_CTX)), "\n")
cat("  - Patients with PFS data:", sum(!is.na(filtered_data$PFS_time)), "\n")

cat("\nResponse Outcomes:\n")
cat("  - VGPR or better (≥VGPR):", sum(filtered_data$VGPR_or_better == 1, na.rm = TRUE),
    "/", sum(!is.na(filtered_data$VGPR_or_better)), "\n")
cat("  - 2-month PFS achieved:", sum(filtered_data$PFS_2month == 1, na.rm = TRUE),
    "/", sum(!is.na(filtered_data$PFS_2month)), "\n")
cat("  - PFS events:", sum(filtered_data$PFS_event == 1, na.rm = TRUE),
    "/", sum(!is.na(filtered_data$PFS_event)), "\n")

# Significant findings
cat("\n>>> Significant Associations (p < 0.05) <<<\n\n")

sig_results <- all_results %>% filter(Sig == "*")

if (nrow(sig_results) > 0) {
  print(sig_results %>% select(Outcome, `PK Metric`, N, `N Events`, `OR (95% CI)`, `p-value`), n = 100)
} else {
  cat("No significant associations found (p < 0.05)\n")
}

cat("\n==========================================================\n")
cat("                    OUTPUT FILES\n")
cat("==========================================================\n")
cat("  - pk_response_analysis_data.csv (filtered analysis data)\n")
cat("  - pk_response_statistical_results.csv (all statistical results)\n")
cat("  - forest_VGPR_OR.png (VGPR forest plot)\n")
cat("  - forest_2monthPFS_OR.png (2-month PFS forest plot)\n")
cat("  - forest_PFS_HR.png (PFS hazard ratio forest plot)\n")
cat("  - boxplot_VGPR_response.png (VGPR boxplots)\n")
cat("==========================================================\n")
cat("                 Analysis Complete!\n")
cat("==========================================================\n")
