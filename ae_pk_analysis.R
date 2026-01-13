#===============================================================================
# Teclistamab PK-AE Association Analysis
#
# Individual PK simulations with Monte Carlo (1000 iterations)
# Compare PK metrics by Adverse Event groups
# Boxplots + Median (IQR) + p-values
#
# AE categories:
# - CRS (all grade)
# - Neurotoxicity (all grade)
# - Hematological toxicity (Grade 3+): Neutropenia, Thrombocytopenia, Lymphopenia
# - Infection (all grade)
# - Neuropathy (all grade)
# - Bilirubinemia (all grade)
# - Liver enzyme elevation (all grade)
# - Creatinine elevation (all grade)
# - Secondary malignancy (all grade)
# - Psychiatric dysfunction (all grade)
#===============================================================================

library(mrgsolve)
library(tidyverse)
library(gridExtra)

set.seed(12345)

cat("==========================================================\n")
cat("    Teclistamab PK-AE Association Analysis\n")
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
mod <- mcode("teclistamab_ae", model_code, compile = TRUE)

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
# 3. Load Data
#-------------------------------------------------------------------------------

cat("\nLoading data files...\n")

# Patient parameters
params <- read_csv("mrgsolve_params_full.csv", show_col_types = FALSE)
cat("  - Loaded", nrow(params), "patients from params file\n")

# Dosing data
dosing_all <- read_csv("mrgsolve_dosing_full.csv", show_col_types = FALSE) %>%
  filter(TIME >= 0) %>%
  arrange(ID, TIME)
cat("  - Loaded", nrow(dosing_all), "dosing records\n")

# AE data (CRS, Neuro from existing file)
ae_crs_neuro <- read_csv("ae_pk_metrics_final.csv", show_col_types = FALSE) %>%
  select(ID, PID, CRS_any, CRS_gr2, Neuro_any, Neuro_gr2,
         Neutropenia_gr3, Thrombocytopenia_gr3, Lymphopenia_gr3, Anemia_gr2, Other_any)
cat("  - Loaded CRS/Neuro AE data for", nrow(ae_crs_neuro), "patients\n")

# Detailed AE data (from Excel extraction)
ae_detailed <- read_csv("ae_detailed_data.csv", show_col_types = FALSE)
cat("  - Loaded detailed AE data for", nrow(ae_detailed), "patients\n")

# Merge AE data
ae_combined <- ae_crs_neuro %>%
  left_join(ae_detailed, by = "PID") %>%
  mutate(
    # Create binary indicators for all-grade AE
    Infection_any = ifelse(AE_INFECTION_GR > 0, 1, 0),
    Neuropathy_any = ifelse(AE_NEUROPATHY_LOAD_GR > 0, 1, 0),
    Bilirubinemia_any = ifelse(AE_BILIRUBINEMIA_GR > 0, 1, 0),
    LiverEnzyme_any = ifelse(AE_LIVER_ENZYME_ELEVATION_GR > 0, 1, 0),
    Creatinine_any = ifelse(AE_CREATININE_ELEVATION_GR > 0, 1, 0),
    SecondaryMalig_any = ifelse(AE_SECONDARY_MALIGNANCY_GR > 0, 1, 0),
    Psychiatric_any = ifelse(AE_PSYCHIATRIC_DYSFUNCTION_GR > 0, 1, 0),
    # Hematological Grade 3+ (from detailed data)
    Neutropenia_gr3_detail = ifelse(AE_NEUTROPENIA_GR >= 3, 1, 0),
    Thrombocytopenia_gr3_detail = ifelse(AE_TROMBOCYTOPENIA_GR >= 3, 1, 0),
    Lymphopenia_gr3_detail = ifelse(AE_LYMPHOPENIA >= 3, 1, 0)
  )

cat("  - Combined AE data ready\n")

#-------------------------------------------------------------------------------
# 4. Monte Carlo Simulation Function
#-------------------------------------------------------------------------------

run_patient_simulation <- function(pt_id, params_df, dosing_df, n_sim = 1000, mod) {
  pt_params <- params_df %>% filter(ID == pt_id)
  if (nrow(pt_params) == 0) return(NULL)

  pt_dosing <- dosing_df %>%
    filter(ID == pt_id) %>%
    arrange(TIME) %>%
    select(time = TIME, amt = AMT, cmt = CMT, evid = EVID)

  if (nrow(pt_dosing) == 0) return(NULL)

  dose_times <- pt_dosing$time
  max_time <- max(pt_dosing$time) + 21 * 24
  sim_times <- seq(0, max_time, by = 1)

  all_sim_results <- vector("list", n_sim)

  for (sim_i in 1:n_sim) {
    eta_CL1 <- rnorm(1, 0, IIV_OMEGA$CL1)
    eta_CL2 <- rnorm(1, 0, IIV_OMEGA$CL2)
    eta_V1  <- rnorm(1, 0, IIV_OMEGA$V1)
    eta_KA  <- rnorm(1, 0, IIV_OMEGA$KA)

    ind_CL1 <- pt_params$CL1 * exp(eta_CL1)
    ind_CL2 <- pt_params$CL2 * exp(eta_CL2)
    ind_V1  <- pt_params$V1 * exp(eta_V1)
    ind_KA  <- pt_params$KA * exp(eta_KA)

    obs_data <- tibble(ID = 1, time = sim_times, amt = 0, cmt = 0, evid = 0)
    dosing_data <- pt_dosing %>% mutate(ID = 1)
    sim_data <- bind_rows(dosing_data, obs_data) %>% arrange(time, desc(evid))

    mod_i <- mod %>%
      param(CL1 = ind_CL1, CL2 = ind_CL2, KDES = pt_params$KDES,
            V1 = ind_V1, V2 = pt_params$V2, Q = pt_params$Q,
            KA = ind_KA, F1 = pt_params$F1)

    out <- mod_i %>%
      data_set(sim_data) %>%
      mrgsim(carry_out = "amt,evid") %>%
      as_tibble() %>%
      filter(evid == 0) %>%
      select(time, DV)

    all_sim_results[[sim_i]] <- out %>% mutate(SIM = sim_i)
  }

  combined <- bind_rows(all_sim_results)
  list(simulations = combined, dose_times = dose_times, pt_info = pt_params)
}

#-------------------------------------------------------------------------------
# 5. PK Metrics Calculation
#-------------------------------------------------------------------------------

calc_metrics_up_to_hour <- function(sim_data, end_hour) {
  sim_data %>%
    filter(time >= 0 & time <= end_hour) %>%
    group_by(SIM) %>%
    summarise(
      Cmax = max(DV, na.rm = TRUE),
      AUC = {
        auc <- 0; dv <- DV; t <- time
        if (length(dv) > 1) {
          for (i in 2:length(dv)) {
            auc <- auc + (t[i] - t[i-1]) / 24 * (dv[i] + dv[i-1]) / 2
          }
        }
        auc
      },
      .groups = "drop"
    ) %>%
    mutate(Cavg = AUC / (end_hour / 24))
}

calc_metrics_dosing_interval <- function(sim_data, dose_times, interval_num) {
  if (length(dose_times) < interval_num) {
    return(tibble(SIM = unique(sim_data$SIM), Cmax = NA_real_, Cavg = NA_real_))
  }
  start_time <- dose_times[interval_num]
  end_time <- if (length(dose_times) > interval_num) dose_times[interval_num + 1] else start_time + 168

  sim_data %>%
    filter(time >= start_time & time <= end_time) %>%
    group_by(SIM) %>%
    summarise(
      Cmax = max(DV, na.rm = TRUE),
      AUC = {
        auc <- 0; dv <- DV; t <- time
        if (length(dv) > 1) {
          for (i in 2:length(dv)) {
            auc <- auc + (t[i] - t[i-1]) / 24 * (dv[i] + dv[i-1]) / 2
          }
        }
        auc
      },
      interval_hours = max(time) - min(time),
      .groups = "drop"
    ) %>%
    mutate(Cavg = ifelse(interval_hours > 0, AUC / (interval_hours / 24), NA_real_))
}

#-------------------------------------------------------------------------------
# 6. Run Analysis for All Patients
#-------------------------------------------------------------------------------

n_simulations <- 1000
cat("\n=== Running Monte Carlo Simulation (n=", n_simulations, ") ===\n\n")

all_pk_results <- list()

for (i in 1:nrow(params)) {
  pt_id <- params$ID[i]
  cat(sprintf("Patient %d/%d (ID: %d)...", i, nrow(params), pt_id))

  mc_result <- run_patient_simulation(pt_id, params, dosing_all, n_sim = n_simulations, mod = mod)

  if (is.null(mc_result)) {
    cat(" Skipped\n")
    next
  }

  sim_data <- mc_result$simulations
  dose_times <- mc_result$dose_times
  pt_info <- mc_result$pt_info

  m_72hr <- calc_metrics_up_to_hour(sim_data, 72)
  m_120hr <- calc_metrics_up_to_hour(sim_data, 120)
  m_dose1 <- calc_metrics_dosing_interval(sim_data, dose_times, 1)
  m_dose3 <- calc_metrics_dosing_interval(sim_data, dose_times, 3)

  pt_summary <- tibble(
    ID = pt_id,
    PID = pt_info$PID,
    WT = pt_info$WT,
    N_doses = length(dose_times),
    Cmax_72hr = median(m_72hr$Cmax, na.rm = TRUE),
    Cavg_72hr = median(m_72hr$Cavg, na.rm = TRUE),
    Cmax_120hr = median(m_120hr$Cmax, na.rm = TRUE),
    Cavg_120hr = median(m_120hr$Cavg, na.rm = TRUE),
    Cmax_dose1 = median(m_dose1$Cmax, na.rm = TRUE),
    Cavg_dose1 = median(m_dose1$Cavg, na.rm = TRUE),
    Cmax_dose3 = median(m_dose3$Cmax, na.rm = TRUE),
    Cavg_dose3 = median(m_dose3$Cavg, na.rm = TRUE)
  )

  all_pk_results[[i]] <- pt_summary
  cat(sprintf(" Cmax_72hr=%.4f, Cavg_72hr=%.4f\n", pt_summary$Cmax_72hr, pt_summary$Cavg_72hr))
}

pk_results <- bind_rows(all_pk_results)
cat("\n=== Completed:", nrow(pk_results), "patients ===\n")

#-------------------------------------------------------------------------------
# 7. Merge PK Results with AE Data
#-------------------------------------------------------------------------------

cat("\nMerging PK results with AE data...\n")

analysis_data <- pk_results %>%
  left_join(ae_combined, by = c("ID", "PID"))

write_csv(analysis_data, "pk_ae_merged_results.csv")
cat("Saved merged data to: pk_ae_merged_results.csv\n")

#-------------------------------------------------------------------------------
# 8. Define AE Groups and PK Metrics
#-------------------------------------------------------------------------------

# AE variables to analyze
ae_vars <- list(
  "CRS (all grade)" = "CRS_any",
  "Neurotoxicity (all grade)" = "Neuro_any",
  "Neutropenia (Gr3+)" = "Neutropenia_gr3",
  "Thrombocytopenia (Gr3+)" = "Thrombocytopenia_gr3",
  "Lymphopenia (Gr3+)" = "Lymphopenia_gr3",
  "Infection (all grade)" = "Infection_any",
  "Neuropathy (all grade)" = "Neuropathy_any",
  "Bilirubinemia (all grade)" = "Bilirubinemia_any",
  "Liver Enzyme Elevation (all grade)" = "LiverEnzyme_any",
  "Creatinine Elevation (all grade)" = "Creatinine_any",
  "Secondary Malignancy (all grade)" = "SecondaryMalig_any",
  "Psychiatric Dysfunction (all grade)" = "Psychiatric_any"
)

pk_metrics <- c("Cmax_72hr", "Cavg_72hr", "Cmax_120hr", "Cavg_120hr",
                "Cmax_dose1", "Cavg_dose1", "Cmax_dose3", "Cavg_dose3")

#-------------------------------------------------------------------------------
# 9. Statistical Analysis Function
#-------------------------------------------------------------------------------

analyze_ae_pk <- function(data, ae_var, ae_label, pk_var) {
  valid_data <- data %>%
    filter(!is.na(.data[[ae_var]]) & !is.na(.data[[pk_var]]))

  if (nrow(valid_data) < 4) return(NULL)

  ae_yes <- valid_data %>% filter(.data[[ae_var]] == 1) %>% pull(.data[[pk_var]])
  ae_no <- valid_data %>% filter(.data[[ae_var]] == 0) %>% pull(.data[[pk_var]])

  if (length(ae_yes) < 1 || length(ae_no) < 1) return(NULL)

  median_iqr <- function(x) {
    if (length(x) == 0) return("N/A")
    sprintf("%.4f (%.4f-%.4f)", median(x), quantile(x, 0.25), quantile(x, 0.75))
  }

  test_result <- tryCatch({
    if (length(ae_yes) >= 2 && length(ae_no) >= 2) {
      wilcox.test(ae_yes, ae_no)
    } else {
      list(p.value = NA)
    }
  }, error = function(e) list(p.value = NA))

  tibble(
    `Adverse Event` = ae_label,
    `PK Metric` = pk_var,
    `AE+ (N)` = length(ae_yes),
    `AE- (N)` = length(ae_no),
    `AE+ Median (IQR)` = median_iqr(ae_yes),
    `AE- Median (IQR)` = median_iqr(ae_no),
    `p-value` = test_result$p.value,
    `Sig` = ifelse(!is.na(test_result$p.value) && test_result$p.value < 0.05, "*", "")
  )
}

#-------------------------------------------------------------------------------
# 10. Run All Statistical Comparisons
#-------------------------------------------------------------------------------

cat("\n=== Statistical Analysis ===\n")

stat_results <- map_dfr(names(ae_vars), function(ae_label) {
  ae_var <- ae_vars[[ae_label]]
  map_dfr(pk_metrics, function(pk_var) {
    analyze_ae_pk(analysis_data, ae_var, ae_label, pk_var)
  })
})

# Filter out NA results
stat_results <- stat_results %>% filter(!is.na(`p-value`) | `AE+ (N)` > 0)

write_csv(stat_results, "ae_pk_statistical_results.csv")
cat("Saved statistical results to: ae_pk_statistical_results.csv\n")

#-------------------------------------------------------------------------------
# 10b. Logistic Regression (OR) for CRS and Neurotoxicity
#-------------------------------------------------------------------------------

cat("\n=== Logistic Regression (OR) for CRS and Neurotoxicity ===\n")

run_ae_logistic <- function(data, ae_var, ae_label, pk_var) {
  model_data <- data %>%
    filter(!is.na(.data[[ae_var]]) & !is.na(.data[[pk_var]]))

  if (nrow(model_data) < 5) return(NULL)
  if (length(unique(model_data[[ae_var]])) < 2) return(NULL)

  # Standardize PK variable
  model_data$pk_std <- scale(model_data[[pk_var]])[,1]

  tryCatch({
    fit <- glm(as.formula(paste(ae_var, "~ pk_std")),
               data = model_data, family = binomial)

    coef_summary <- summary(fit)$coefficients

    beta <- coef_summary["pk_std", "Estimate"]
    se <- coef_summary["pk_std", "Std. Error"]
    p_value <- coef_summary["pk_std", "Pr(>|z|)"]

    or <- exp(beta)
    or_lower <- exp(beta - 1.96 * se)
    or_upper <- exp(beta + 1.96 * se)

    tibble(
      `Adverse Event` = ae_label,
      `PK Metric` = pk_var,
      N = nrow(model_data),
      `AE+` = sum(model_data[[ae_var]] == 1),
      `OR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", or, or_lower, or_upper),
      OR = or,
      `p-value` = p_value,
      Sig = ifelse(p_value < 0.05, "*", "")
    )
  }, error = function(e) NULL)
}

# Run OR analysis for CRS and Neurotoxicity
ae_or_vars <- list(
  "CRS (all grade)" = "CRS_any",
  "CRS (Grade 2+)" = "CRS_gr2",
  "Neurotoxicity (all grade)" = "Neuro_any",
  "Neurotoxicity (Grade 2+)" = "Neuro_gr2"
)

or_results <- map_dfr(names(ae_or_vars), function(ae_label) {
  ae_var <- ae_or_vars[[ae_label]]
  map_dfr(pk_metrics, function(pk_var) {
    run_ae_logistic(analysis_data, ae_var, ae_label, pk_var)
  })
})

if (nrow(or_results) > 0) {
  cat("\n--- Logistic Regression OR Results ---\n")
  print(or_results %>% select(-OR), n = 100)

  write_csv(or_results, "ae_pk_or_results.csv")
  cat("\nSaved: ae_pk_or_results.csv\n")
}

#-------------------------------------------------------------------------------
# 10c. Categorical Analysis (Median Split) for CRS and Neurotoxicity
#-------------------------------------------------------------------------------

cat("\n=== Categorical Analysis (Median Split) ===\n")

# Create median-split categories
analysis_data_cat <- analysis_data
for (pk_var in pk_metrics) {
  med_val <- median(analysis_data[[pk_var]], na.rm = TRUE)
  cat_var <- paste0(pk_var, "_cat")
  analysis_data_cat[[cat_var]] <- ifelse(analysis_data[[pk_var]] >= med_val, "High", "Low")
  analysis_data_cat[[cat_var]] <- factor(analysis_data_cat[[cat_var]], levels = c("Low", "High"))
}

run_ae_categorical_logistic <- function(data, ae_var, ae_label, pk_var, orig_data) {
  cat_var <- paste0(pk_var, "_cat")

  model_data <- data %>%
    filter(!is.na(.data[[ae_var]]) & !is.na(.data[[cat_var]]))

  if (nrow(model_data) < 5) return(NULL)
  if (length(unique(model_data[[ae_var]])) < 2) return(NULL)
  if (length(unique(model_data[[cat_var]])) < 2) return(NULL)

  med_val <- median(orig_data[[pk_var]], na.rm = TRUE)

  tryCatch({
    fit <- glm(as.formula(paste(ae_var, "~", cat_var)),
               data = model_data, family = binomial)

    coef_summary <- summary(fit)$coefficients

    beta <- coef_summary[2, "Estimate"]
    se <- coef_summary[2, "Std. Error"]
    p_value <- coef_summary[2, "Pr(>|z|)"]

    or <- exp(beta)
    or_lower <- exp(beta - 1.96 * se)
    or_upper <- exp(beta + 1.96 * se)

    # Count by group
    high_ae <- sum(model_data[[cat_var]] == "High" & model_data[[ae_var]] == 1)
    high_n <- sum(model_data[[cat_var]] == "High")
    low_ae <- sum(model_data[[cat_var]] == "Low" & model_data[[ae_var]] == 1)
    low_n <- sum(model_data[[cat_var]] == "Low")

    tibble(
      `Adverse Event` = ae_label,
      `PK Metric` = pk_var,
      `Median Cutoff` = round(med_val, 4),
      `High (AE+/N)` = sprintf("%d/%d", high_ae, high_n),
      `Low (AE+/N)` = sprintf("%d/%d", low_ae, low_n),
      `OR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", or, or_lower, or_upper),
      OR = or,
      `p-value` = p_value,
      Sig = ifelse(p_value < 0.05, "*", "")
    )
  }, error = function(e) NULL)
}

cat_or_results <- map_dfr(names(ae_or_vars), function(ae_label) {
  ae_var <- ae_or_vars[[ae_label]]
  map_dfr(pk_metrics, function(pk_var) {
    run_ae_categorical_logistic(analysis_data_cat, ae_var, ae_label, pk_var, analysis_data)
  })
})

if (nrow(cat_or_results) > 0) {
  cat("\n--- Categorical OR Results (Median Split) ---\n")
  print(cat_or_results %>% select(-OR), n = 100)

  write_csv(cat_or_results, "ae_pk_categorical_or_results.csv")
  cat("\nSaved: ae_pk_categorical_or_results.csv\n")
}

#-------------------------------------------------------------------------------
# 10d. Predictive Performance for CRS and Neurotoxicity
#-------------------------------------------------------------------------------

cat("\n=== Predictive Performance for AE ===\n")

library(pROC)

# Collect all p-values
ae_all_pvals <- bind_rows(
  stat_results %>%
    filter(`Adverse Event` %in% c("CRS (all grade)", "Neurotoxicity (all grade)")) %>%
    select(`PK Metric`, `p-value`) %>%
    mutate(Analysis = "Wilcoxon"),
  or_results %>%
    filter(`Adverse Event` %in% c("CRS (all grade)", "Neurotoxicity (all grade)")) %>%
    select(`PK Metric`, `p-value`) %>%
    mutate(Analysis = "OR_continuous"),
  cat_or_results %>%
    filter(`Adverse Event` %in% c("CRS (all grade)", "Neurotoxicity (all grade)")) %>%
    select(`PK Metric`, `p-value`) %>%
    mutate(Analysis = "OR_categorical")
)

# Rank metrics
ae_sig_summary <- ae_all_pvals %>%
  group_by(`PK Metric`) %>%
  summarise(
    Mean_pvalue = mean(`p-value`, na.rm = TRUE),
    Min_pvalue = min(`p-value`, na.rm = TRUE),
    N_significant = sum(`p-value` < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Min_pvalue)

cat("\n--- AE PK Metric Significance Ranking ---\n")
print(ae_sig_summary, n = 10)

# Top 2 metrics
ae_top_metrics <- ae_sig_summary %>% head(2) %>% pull(`PK Metric`)
cat("\n>>> Top 2 Most Significant PK Metrics for AE:", paste(ae_top_metrics, collapse = ", "), "<<<\n")

# ROC Analysis
cat("\n--- ROC Analysis for AE ---\n")

ae_roc_results <- list()

for (pk_var in ae_top_metrics) {
  cat(sprintf("\n%s:\n", pk_var))

  # ROC for CRS
  crs_data <- analysis_data %>%
    filter(!is.na(CRS_any) & !is.na(.data[[pk_var]]))

  if (nrow(crs_data) >= 5 && length(unique(crs_data$CRS_any)) == 2) {
    roc_crs <- tryCatch({
      roc(crs_data$CRS_any, crs_data[[pk_var]], quiet = TRUE)
    }, error = function(e) NULL)

    if (!is.null(roc_crs)) {
      auc_crs <- auc(roc_crs)
      ci_crs <- ci.auc(roc_crs)
      coords_crs <- coords(roc_crs, "best", ret = c("threshold", "sensitivity", "specificity"))

      cat(sprintf("  CRS: AUC = %.3f (%.3f-%.3f)\n", auc_crs, ci_crs[1], ci_crs[3]))
      cat(sprintf("       Optimal cutoff = %.4f (Sens=%.2f, Spec=%.2f)\n",
                  coords_crs$threshold, coords_crs$sensitivity, coords_crs$specificity))

      ae_roc_results[[paste0(pk_var, "_CRS")]] <- list(
        metric = pk_var, outcome = "CRS",
        auc = as.numeric(auc_crs), auc_lower = ci_crs[1], auc_upper = ci_crs[3],
        cutoff = coords_crs$threshold,
        sensitivity = coords_crs$sensitivity,
        specificity = coords_crs$specificity
      )
    }
  }

  # ROC for Neurotoxicity
  neuro_data <- analysis_data %>%
    filter(!is.na(Neuro_any) & !is.na(.data[[pk_var]]))

  if (nrow(neuro_data) >= 5 && length(unique(neuro_data$Neuro_any)) == 2) {
    roc_neuro <- tryCatch({
      roc(neuro_data$Neuro_any, neuro_data[[pk_var]], quiet = TRUE)
    }, error = function(e) NULL)

    if (!is.null(roc_neuro)) {
      auc_neuro <- auc(roc_neuro)
      ci_neuro <- ci.auc(roc_neuro)
      coords_neuro <- coords(roc_neuro, "best", ret = c("threshold", "sensitivity", "specificity"))

      cat(sprintf("  Neurotoxicity: AUC = %.3f (%.3f-%.3f)\n", auc_neuro, ci_neuro[1], ci_neuro[3]))
      cat(sprintf("                 Optimal cutoff = %.4f (Sens=%.2f, Spec=%.2f)\n",
                  coords_neuro$threshold, coords_neuro$sensitivity, coords_neuro$specificity))

      ae_roc_results[[paste0(pk_var, "_Neuro")]] <- list(
        metric = pk_var, outcome = "Neurotoxicity",
        auc = as.numeric(auc_neuro), auc_lower = ci_neuro[1], auc_upper = ci_neuro[3],
        cutoff = coords_neuro$threshold,
        sensitivity = coords_neuro$sensitivity,
        specificity = coords_neuro$specificity
      )
    }
  }
}

# Summary table
if (length(ae_roc_results) > 0) {
  ae_roc_summary <- map_dfr(ae_roc_results, ~ tibble(
    `PK Metric` = .x$metric,
    `Adverse Event` = .x$outcome,
    `AUC (95% CI)` = sprintf("%.3f (%.3f-%.3f)", .x$auc, .x$auc_lower, .x$auc_upper),
    `Optimal Cutoff` = round(.x$cutoff, 4),
    Sensitivity = round(.x$sensitivity, 2),
    Specificity = round(.x$specificity, 2)
  ))

  cat("\n--- AE Predictive Performance Summary ---\n")
  print(ae_roc_summary, n = 100)
  write_csv(ae_roc_summary, "ae_pk_roc_results.csv")
  cat("\nSaved: ae_pk_roc_results.csv\n")
}

#-------------------------------------------------------------------------------
# 11. Create Boxplots
#-------------------------------------------------------------------------------

cat("\n=== Creating Boxplots ===\n")

create_ae_boxplot <- function(data, ae_var, ae_label, pk_var) {
  plot_data <- data %>%
    filter(!is.na(.data[[ae_var]]) & !is.na(.data[[pk_var]])) %>%
    mutate(AE_Status = factor(ifelse(.data[[ae_var]] == 1, "Yes", "No"), levels = c("No", "Yes")))

  if (nrow(plot_data) < 4) return(NULL)

  # Get statistics
  stats <- stat_results %>%
    filter(`Adverse Event` == ae_label & `PK Metric` == pk_var)

  p_val <- if (nrow(stats) > 0 && !is.na(stats$`p-value`[1])) stats$`p-value`[1] else NA

  p_label <- if (!is.na(p_val)) {
    if (p_val < 0.001) "p < 0.001"
    else if (p_val < 0.01) sprintf("p = %.3f", p_val)
    else sprintf("p = %.2f", p_val)
  } else "p = NA"

  y_max <- max(plot_data[[pk_var]], na.rm = TRUE)

  ggplot(plot_data, aes(x = AE_Status, y = .data[[pk_var]], fill = AE_Status)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21) +
    geom_jitter(width = 0.15, alpha = 0.6, size = 2.5) +
    scale_fill_manual(values = c("No" = "#3498db", "Yes" = "#e74c3c")) +
    labs(title = ae_label, x = "", y = pk_var) +
    annotate("text", x = 1.5, y = y_max * 1.15, label = p_label, size = 3.5, fontface = "bold") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(size = 9, face = "bold"),
          axis.title.y = element_text(size = 9))
}

# Create plots for each PK metric
for (pk_var in pk_metrics) {
  cat(sprintf("Creating boxplots for %s...\n", pk_var))

  plots <- list()
  for (ae_label in names(ae_vars)) {
    ae_var <- ae_vars[[ae_label]]
    p <- create_ae_boxplot(analysis_data, ae_var, ae_label, pk_var)
    if (!is.null(p)) plots[[ae_label]] <- p
  }

  if (length(plots) > 0) {
    n_plots <- length(plots)
    ncol <- 4
    nrow <- ceiling(n_plots / ncol)

    combined_plot <- arrangeGrob(grobs = plots, ncol = ncol, nrow = nrow,
                                  top = paste0(pk_var, " by Adverse Event Status"))
    ggsave(sprintf("boxplot_%s_by_AE.png", pk_var), combined_plot,
           width = 16, height = 4 * nrow, dpi = 200)
  }
}

#-------------------------------------------------------------------------------
# 12. Summary Figure for Key Metrics (CRS focus)
#-------------------------------------------------------------------------------

cat("\nCreating summary figure for CRS...\n")

key_pk_metrics <- c("Cavg_72hr", "Cavg_120hr", "Cmax_72hr", "Cmax_120hr")
crs_plots <- list()

for (pk_var in key_pk_metrics) {
  p <- create_ae_boxplot(analysis_data, "CRS_any", "CRS (all grade)", pk_var)
  if (!is.null(p)) {
    p <- p + labs(title = pk_var, x = "CRS")
    crs_plots[[pk_var]] <- p
  }
}

if (length(crs_plots) > 0) {
  crs_combined <- arrangeGrob(grobs = crs_plots, ncol = 2, nrow = 2,
                               top = "PK Metrics by CRS Status (All Grade)")
  ggsave("boxplot_CRS_summary.png", crs_combined, width = 10, height = 10, dpi = 200)
}

#-------------------------------------------------------------------------------
# 13. Print Final Summary
#-------------------------------------------------------------------------------

cat("\n")
cat("==========================================================\n")
cat("                    FINAL SUMMARY\n")
cat("==========================================================\n\n")

# Significant findings
sig_results <- stat_results %>% filter(`Sig` == "*")

if (nrow(sig_results) > 0) {
  cat(">>> Significant Associations (p < 0.05) <<<\n\n")
  print(sig_results %>% select(-Sig), n = 100)
} else {
  cat("No significant associations found (p < 0.05)\n")
}

cat("\n")
cat("==========================================================\n")
cat("              COMPLETE STATISTICAL RESULTS\n")
cat("==========================================================\n\n")

# Print full results table grouped by AE
for (ae_label in unique(stat_results$`Adverse Event`)) {
  cat(sprintf("\n--- %s ---\n", ae_label))
  ae_subset <- stat_results %>%
    filter(`Adverse Event` == ae_label) %>%
    select(`PK Metric`, `AE+ (N)`, `AE- (N)`, `AE+ Median (IQR)`, `AE- Median (IQR)`, `p-value`)
  print(ae_subset, n = 100)
}

cat("\n==========================================================\n")
cat("                    OUTPUT FILES\n")
cat("==========================================================\n")
cat("  - pk_ae_merged_results.csv (merged PK and AE data)\n")
cat("  - ae_pk_statistical_results.csv (Wilcoxon test results)\n")
cat("  - ae_pk_or_results.csv (logistic regression OR results)\n")
cat("  - ae_pk_categorical_or_results.csv (median-split categorical OR)\n")
cat("  - ae_pk_roc_results.csv (ROC/AUC predictive performance)\n")
cat("  - boxplot_*_by_AE.png (boxplots for each PK metric)\n")
cat("  - boxplot_CRS_summary.png (CRS summary)\n")
cat("==========================================================\n")
cat("                 Analysis Complete!\n")
cat("==========================================================\n")
