#===============================================================================
# Teclistamab PK-Response Association Analysis
#
# For patients with ≥4 weeks of treatment:
# 1. VGPR or better (≥VGPR) - Odds Ratio (Logistic Regression)
# 2. 2-month PFS status - Odds Ratio (Logistic Regression)
# 3. PFS - Hazard Ratio (Cox Regression)
#
# Additional analyses:
# - Boxplots with Wilcoxon p-values
# - Categorical analysis (median split) with OR/HR
# - Predictive performance evaluation
#
# NOTE: Uses PK metrics from ae_pk_analysis.R (pk_ae_merged_results.csv)
#       Run ae_pk_analysis.R first to generate the PK data
#===============================================================================

library(tidyverse)
library(survival)
library(broom)
library(gridExtra)
library(pROC)

set.seed(12345)

cat("==========================================================\n")
cat("    Teclistamab PK-Response Association Analysis\n")
cat("==========================================================\n\n")

#-------------------------------------------------------------------------------
# 1. Load PK Data from Previous Analysis
#-------------------------------------------------------------------------------

cat("Loading data...\n")

# Load PK metrics from ae_pk_analysis.R output
if (!file.exists("pk_ae_merged_results.csv")) {
  stop("pk_ae_merged_results.csv not found! Run ae_pk_analysis.R first.")
}

pk_data <- read_csv("pk_ae_merged_results.csv", show_col_types = FALSE)
cat("  - Loaded PK data for", nrow(pk_data), "patients from ae_pk_analysis.R\n")

# Load dosing data to calculate treatment duration
dosing_all <- read_csv("mrgsolve_dosing_full.csv", show_col_types = FALSE) %>%
  filter(TIME >= 0) %>%
  arrange(ID, TIME)

# Load response data
response_data <- read_csv("response_data.csv", show_col_types = FALSE)
cat("  - Loaded response data for", nrow(response_data), "patients\n")

#-------------------------------------------------------------------------------
# 2. Calculate Treatment Duration and Filter for ≥4 Weeks
#-------------------------------------------------------------------------------

cat("\n=== Filtering for ≥4 weeks treatment ===\n")

# Calculate treatment duration per patient (max TIME in hours)
treatment_duration <- dosing_all %>%
  group_by(ID) %>%
  summarise(
    N_doses = n(),
    Max_TIME_hr = max(TIME),
    Treatment_days = max(TIME) / 24,
    .groups = "drop"
  )

cat("\nTreatment Duration Summary:\n")
print(treatment_duration)

# Filter for ≥4 weeks treatment (28 days = 672 hours)
patients_4weeks <- treatment_duration %>%
  filter(Treatment_days >= 28)

cat("\n>>> Patients with ≥4 weeks treatment:", nrow(patients_4weeks), "out of", nrow(treatment_duration), "<<<\n")

if (nrow(patients_4weeks) == 0) {
  stop("No patients with ≥4 weeks of treatment found!")
}

cat("\nFiltered patients IDs:", paste(patients_4weeks$ID, collapse = ", "), "\n")

#-------------------------------------------------------------------------------
# 3. Merge Data and Create Analysis Dataset
#-------------------------------------------------------------------------------

cat("\n=== Creating analysis dataset ===\n")

# Filter PK data for ≥4 weeks patients and merge with response data
filtered_data <- pk_data %>%
  filter(ID %in% patients_4weeks$ID) %>%
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
# 4. Define PK Metrics for Analysis
#-------------------------------------------------------------------------------

pk_metrics <- c("Cmax_72hr", "Cavg_72hr", "Cmax_120hr", "Cavg_120hr",
                "Cmax_dose1", "Cavg_dose1", "Cmax_dose3", "Cavg_dose3")

#-------------------------------------------------------------------------------
# 5. Logistic Regression Function (for OR)
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
# 6. Cox Regression Function (for HR)
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
# 7. Run Analyses
#-------------------------------------------------------------------------------

cat("\n")
cat("==========================================================\n")
cat("             STATISTICAL ANALYSIS RESULTS\n")
cat("==========================================================\n")

# 7.1 VGPR or better (Logistic Regression - OR)
cat("\n\n--- 1. VGPR or Better (≥VGPR) - Odds Ratio ---\n")

vgpr_results <- map_dfr(pk_metrics, function(pk_var) {
  run_logistic_analysis(filtered_data, "VGPR_or_better", pk_var, "≥VGPR")
})

if (nrow(vgpr_results) > 0) {
  print(vgpr_results %>% select(-OR, -OR_lower, -OR_upper), n = 100)
} else {
  cat("Insufficient data for VGPR analysis\n")
}

# 7.2 2-month PFS (Logistic Regression - OR)
cat("\n\n--- 2. 2-Month PFS Status - Odds Ratio ---\n")

pfs2m_results <- map_dfr(pk_metrics, function(pk_var) {
  run_logistic_analysis(filtered_data, "PFS_2month", pk_var, "2-month PFS")
})

if (nrow(pfs2m_results) > 0) {
  print(pfs2m_results %>% select(-OR, -OR_lower, -OR_upper), n = 100)
} else {
  cat("Insufficient data for 2-month PFS analysis\n")
}

# 7.3 PFS (Cox Regression - HR)
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
# 8. Boxplots with Wilcoxon p-values for Response Outcomes
#-------------------------------------------------------------------------------

cat("\n\n=== Boxplots with Wilcoxon p-values ===\n")

# Function to create boxplot with p-value (like ae_pk_analysis.R)
create_response_boxplot_pval <- function(data, outcome_var, outcome_label, pk_var) {
  plot_data <- data %>%
    filter(!is.na(.data[[outcome_var]]) & !is.na(.data[[pk_var]])) %>%
    mutate(Response = factor(ifelse(.data[[outcome_var]] == 1, "Yes", "No"),
                              levels = c("No", "Yes")))

  if (nrow(plot_data) < 4) return(list(plot = NULL, stats = NULL))

  yes_vals <- plot_data %>% filter(Response == "Yes") %>% pull(.data[[pk_var]])
  no_vals <- plot_data %>% filter(Response == "No") %>% pull(.data[[pk_var]])

  if (length(yes_vals) < 1 || length(no_vals) < 1) return(list(plot = NULL, stats = NULL))

  # Wilcoxon test
  p_val <- tryCatch({
    wilcox.test(yes_vals, no_vals)$p.value
  }, error = function(e) NA)

  # Median (IQR)
  median_iqr <- function(x) {
    sprintf("%.4f (%.4f-%.4f)", median(x), quantile(x, 0.25), quantile(x, 0.75))
  }

  p_label <- if (!is.na(p_val)) {
    if (p_val < 0.001) "p < 0.001"
    else if (p_val < 0.01) sprintf("p = %.3f", p_val)
    else sprintf("p = %.2f", p_val)
  } else "p = NA"

  y_max <- max(plot_data[[pk_var]], na.rm = TRUE)

  p <- ggplot(plot_data, aes(x = Response, y = .data[[pk_var]], fill = Response)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21) +
    geom_jitter(width = 0.15, alpha = 0.6, size = 2.5) +
    scale_fill_manual(values = c("No" = "#e74c3c", "Yes" = "#27ae60")) +
    labs(title = pk_var, x = outcome_label, y = pk_var) +
    annotate("text", x = 1.5, y = y_max * 1.15, label = p_label, size = 3.5, fontface = "bold") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(size = 10, face = "bold"))

  stats <- tibble(
    Outcome = outcome_label,
    `PK Metric` = pk_var,
    `Yes (N)` = length(yes_vals),
    `No (N)` = length(no_vals),
    `Yes Median (IQR)` = median_iqr(yes_vals),
    `No Median (IQR)` = median_iqr(no_vals),
    `Wilcoxon p-value` = p_val,
    Sig = ifelse(!is.na(p_val) && p_val < 0.05, "*", "")
  )

  list(plot = p, stats = stats)
}

# Create boxplots for VGPR
cat("\n--- VGPR or Better Boxplots ---\n")
vgpr_boxplot_results <- map(pk_metrics, function(pk_var) {
  create_response_boxplot_pval(filtered_data, "VGPR_or_better", "≥VGPR", pk_var)
})

vgpr_boxplot_stats <- map_dfr(vgpr_boxplot_results, ~ .x$stats)
vgpr_plots <- compact(map(vgpr_boxplot_results, ~ .x$plot))

if (nrow(vgpr_boxplot_stats) > 0) {
  print(vgpr_boxplot_stats, n = 100)
}

if (length(vgpr_plots) > 0) {
  combined <- arrangeGrob(grobs = vgpr_plots, ncol = 4, nrow = 2,
                          top = "PK Metrics by VGPR Response (with Wilcoxon p-values)")
  ggsave("boxplot_VGPR_wilcoxon.png", combined, width = 16, height = 8, dpi = 200)
  cat("Saved: boxplot_VGPR_wilcoxon.png\n")
}

# Create boxplots for 2-month PFS
cat("\n--- 2-Month PFS Boxplots ---\n")
pfs2m_boxplot_results <- map(pk_metrics, function(pk_var) {
  create_response_boxplot_pval(filtered_data, "PFS_2month", "2-month PFS", pk_var)
})

pfs2m_boxplot_stats <- map_dfr(pfs2m_boxplot_results, ~ .x$stats)
pfs2m_plots <- compact(map(pfs2m_boxplot_results, ~ .x$plot))

if (nrow(pfs2m_boxplot_stats) > 0) {
  print(pfs2m_boxplot_stats, n = 100)
}

if (length(pfs2m_plots) > 0) {
  combined <- arrangeGrob(grobs = pfs2m_plots, ncol = 4, nrow = 2,
                          top = "PK Metrics by 2-Month PFS Status (with Wilcoxon p-values)")
  ggsave("boxplot_2monthPFS_wilcoxon.png", combined, width = 16, height = 8, dpi = 200)
  cat("Saved: boxplot_2monthPFS_wilcoxon.png\n")
}

#-------------------------------------------------------------------------------
# 9. Categorical Analysis (Median Split) - OR/HR
#-------------------------------------------------------------------------------

cat("\n\n=== Categorical Analysis (Median Split) ===\n")

# Create median-split categories for each PK metric
filtered_data_cat <- filtered_data
for (pk_var in pk_metrics) {
  med_val <- median(filtered_data[[pk_var]], na.rm = TRUE)
  cat_var <- paste0(pk_var, "_cat")
  filtered_data_cat[[cat_var]] <- ifelse(filtered_data[[pk_var]] >= med_val, "High", "Low")
  filtered_data_cat[[cat_var]] <- factor(filtered_data_cat[[cat_var]], levels = c("Low", "High"))
}

# Logistic regression for categorical variables
run_categorical_logistic <- function(data, outcome_var, pk_var, outcome_name) {
  cat_var <- paste0(pk_var, "_cat")

  model_data <- data %>%
    filter(!is.na(.data[[outcome_var]]) & !is.na(.data[[cat_var]]))

  if (nrow(model_data) < 5) return(NULL)
  if (length(unique(model_data[[outcome_var]])) < 2) return(NULL)
  if (length(unique(model_data[[cat_var]])) < 2) return(NULL)

  med_val <- median(filtered_data[[pk_var]], na.rm = TRUE)

  tryCatch({
    fit <- glm(as.formula(paste(outcome_var, "~", cat_var)),
               data = model_data, family = binomial)

    coef_summary <- summary(fit)$coefficients

    beta <- coef_summary[2, "Estimate"]
    se <- coef_summary[2, "Std. Error"]
    p_value <- coef_summary[2, "Pr(>|z|)"]

    or <- exp(beta)
    or_lower <- exp(beta - 1.96 * se)
    or_upper <- exp(beta + 1.96 * se)

    # Count by group
    high_yes <- sum(model_data[[cat_var]] == "High" & model_data[[outcome_var]] == 1)
    high_n <- sum(model_data[[cat_var]] == "High")
    low_yes <- sum(model_data[[cat_var]] == "Low" & model_data[[outcome_var]] == 1)
    low_n <- sum(model_data[[cat_var]] == "Low")

    tibble(
      Outcome = outcome_name,
      `PK Metric` = pk_var,
      `Median Cutoff` = round(med_val, 4),
      `High (events/N)` = sprintf("%d/%d", high_yes, high_n),
      `Low (events/N)` = sprintf("%d/%d", low_yes, low_n),
      `OR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", or, or_lower, or_upper),
      OR = or,
      `p-value` = p_value,
      Sig = ifelse(p_value < 0.05, "*", "")
    )
  }, error = function(e) NULL)
}

# Cox regression for categorical variables
run_categorical_cox <- function(data, pk_var) {
  cat_var <- paste0(pk_var, "_cat")

  model_data <- data %>%
    filter(!is.na(PFS_event) & !is.na(PFS_time) & !is.na(.data[[cat_var]])) %>%
    filter(PFS_time > 0)

  if (nrow(model_data) < 5) return(NULL)
  if (sum(model_data$PFS_event) < 2) return(NULL)
  if (length(unique(model_data[[cat_var]])) < 2) return(NULL)

  med_val <- median(filtered_data[[pk_var]], na.rm = TRUE)

  tryCatch({
    fit <- coxph(as.formula(paste("Surv(PFS_time, PFS_event) ~", cat_var)), data = model_data)
    fit_summary <- summary(fit)

    hr <- fit_summary$conf.int[1, "exp(coef)"]
    hr_lower <- fit_summary$conf.int[1, "lower .95"]
    hr_upper <- fit_summary$conf.int[1, "upper .95"]
    p_value <- fit_summary$coefficients[1, "Pr(>|z|)"]

    # Count by group
    high_events <- sum(model_data[[cat_var]] == "High" & model_data$PFS_event == 1)
    high_n <- sum(model_data[[cat_var]] == "High")
    low_events <- sum(model_data[[cat_var]] == "Low" & model_data$PFS_event == 1)
    low_n <- sum(model_data[[cat_var]] == "Low")

    tibble(
      Outcome = "PFS",
      `PK Metric` = pk_var,
      `Median Cutoff` = round(med_val, 4),
      `High (events/N)` = sprintf("%d/%d", high_events, high_n),
      `Low (events/N)` = sprintf("%d/%d", low_events, low_n),
      `HR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", hr, hr_lower, hr_upper),
      HR = hr,
      `p-value` = p_value,
      Sig = ifelse(p_value < 0.05, "*", "")
    )
  }, error = function(e) NULL)
}

# Run categorical analysis for VGPR
cat("\n--- VGPR (Categorical - Median Split) ---\n")
vgpr_cat_results <- map_dfr(pk_metrics, function(pk_var) {
  run_categorical_logistic(filtered_data_cat, "VGPR_or_better", pk_var, "≥VGPR")
})
if (nrow(vgpr_cat_results) > 0) {
  print(vgpr_cat_results %>% select(-OR), n = 100)
}

# Run categorical analysis for 2-month PFS
cat("\n--- 2-Month PFS (Categorical - Median Split) ---\n")
pfs2m_cat_results <- map_dfr(pk_metrics, function(pk_var) {
  run_categorical_logistic(filtered_data_cat, "PFS_2month", pk_var, "2-month PFS")
})
if (nrow(pfs2m_cat_results) > 0) {
  print(pfs2m_cat_results %>% select(-OR), n = 100)
}

# Run categorical analysis for PFS (Cox)
cat("\n--- PFS (Categorical - Median Split, Cox HR) ---\n")
pfs_cat_results <- map_dfr(pk_metrics, function(pk_var) {
  run_categorical_cox(filtered_data_cat, pk_var)
})
if (nrow(pfs_cat_results) > 0) {
  print(pfs_cat_results %>% select(-HR), n = 100)
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
  pfs_plot_data <- pfs_results 
  p_pfs <- create_forest_plot(pfs_plot_data, "Progression-Free Survival - Hazard Ratio", "HR")
  if (!is.null(p_pfs)) {
    ggsave("forest_PFS_HR.png", p_pfs, width = 10, height = 6, dpi = 200)
    cat("Saved: forest_PFS_HR.png\n")
  }
}

#-------------------------------------------------------------------------------
# 12. Predictive Performance Evaluation
#-------------------------------------------------------------------------------

cat("\n\n=== Predictive Performance Evaluation ===\n")

# Collect all p-values from different analyses
all_pvals <- bind_rows(
  # Continuous logistic regression
  vgpr_results %>% select(`PK Metric`, `p-value`) %>% mutate(Analysis = "VGPR_continuous"),
  pfs2m_results %>% select(`PK Metric`, `p-value`) %>% mutate(Analysis = "PFS2m_continuous"),
  # Wilcoxon tests
  vgpr_boxplot_stats %>% select(`PK Metric`, `p-value` = `Wilcoxon p-value`) %>% mutate(Analysis = "VGPR_wilcoxon"),
  pfs2m_boxplot_stats %>% select(`PK Metric`, `p-value` = `Wilcoxon p-value`) %>% mutate(Analysis = "PFS2m_wilcoxon"),
  # Categorical
  vgpr_cat_results %>% select(`PK Metric`, `p-value`) %>% mutate(Analysis = "VGPR_categorical"),
  pfs2m_cat_results %>% select(`PK Metric`, `p-value`) %>% mutate(Analysis = "PFS2m_categorical")
)

# Find most significant metrics (lowest average p-value or most frequent significance)
sig_summary <- all_pvals %>%
  group_by(`PK Metric`) %>%
  summarise(
    Mean_pvalue = mean(`p-value`, na.rm = TRUE),
    Min_pvalue = min(`p-value`, na.rm = TRUE),
    N_significant = sum(`p-value` < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Min_pvalue)

cat("\n--- PK Metric Significance Ranking ---\n")
print(sig_summary, n = 10)

# Select top 2 most significant metrics
top_metrics <- sig_summary %>% head(2) %>% pull(`PK Metric`)
cat("\n>>> Top 2 Most Significant PK Metrics:", paste(top_metrics, collapse = ", "), "<<<\n")

# ROC Analysis for top metrics
cat("\n--- ROC Analysis for Top Metrics ---\n")

roc_results <- list()

for (pk_var in top_metrics) {
  cat(sprintf("\n%s:\n", pk_var))

  # ROC for VGPR
  vgpr_data <- filtered_data %>%
    filter(!is.na(VGPR_or_better) & !is.na(.data[[pk_var]]))

  if (nrow(vgpr_data) >= 5 && length(unique(vgpr_data$VGPR_or_better)) == 2) {
    roc_vgpr <- tryCatch({
      roc(vgpr_data$VGPR_or_better, vgpr_data[[pk_var]], quiet = TRUE)
    }, error = function(e) NULL)

    if (!is.null(roc_vgpr)) {
      auc_vgpr <- auc(roc_vgpr)
      ci_vgpr <- ci.auc(roc_vgpr)

      # Find optimal cutoff (Youden)
      coords_vgpr <- coords(roc_vgpr, "best", ret = c("threshold", "sensitivity", "specificity"))

      cat(sprintf("  VGPR: AUC = %.3f (%.3f-%.3f)\n", auc_vgpr, ci_vgpr[1], ci_vgpr[3]))
      cat(sprintf("         Optimal cutoff = %.4f (Sens=%.2f, Spec=%.2f)\n",
                  coords_vgpr$threshold, coords_vgpr$sensitivity, coords_vgpr$specificity))

      roc_results[[paste0(pk_var, "_VGPR")]] <- list(
        metric = pk_var, outcome = "VGPR",
        auc = as.numeric(auc_vgpr), auc_lower = ci_vgpr[1], auc_upper = ci_vgpr[3],
        cutoff = coords_vgpr$threshold,
        sensitivity = coords_vgpr$sensitivity,
        specificity = coords_vgpr$specificity
      )
    }
  }

  # ROC for 2-month PFS
  pfs2m_data <- filtered_data %>%
    filter(!is.na(PFS_2month) & !is.na(.data[[pk_var]]))

  if (nrow(pfs2m_data) >= 5 && length(unique(pfs2m_data$PFS_2month)) == 2) {
    roc_pfs2m <- tryCatch({
      roc(pfs2m_data$PFS_2month, pfs2m_data[[pk_var]], quiet = TRUE)
    }, error = function(e) NULL)

    if (!is.null(roc_pfs2m)) {
      auc_pfs2m <- auc(roc_pfs2m)
      ci_pfs2m <- ci.auc(roc_pfs2m)
      coords_pfs2m <- coords(roc_pfs2m, "best", ret = c("threshold", "sensitivity", "specificity"))

      cat(sprintf("  2-mo PFS: AUC = %.3f (%.3f-%.3f)\n", auc_pfs2m, ci_pfs2m[1], ci_pfs2m[3]))
      cat(sprintf("            Optimal cutoff = %.4f (Sens=%.2f, Spec=%.2f)\n",
                  coords_pfs2m$threshold, coords_pfs2m$sensitivity, coords_pfs2m$specificity))

      roc_results[[paste0(pk_var, "_PFS2m")]] <- list(
        metric = pk_var, outcome = "2-month PFS",
        auc = as.numeric(auc_pfs2m), auc_lower = ci_pfs2m[1], auc_upper = ci_pfs2m[3],
        cutoff = coords_pfs2m$threshold,
        sensitivity = coords_pfs2m$sensitivity,
        specificity = coords_pfs2m$specificity
      )
    }
  }
}

# Summary table
if (length(roc_results) > 0) {
  roc_summary <- map_dfr(roc_results, ~ tibble(
    `PK Metric` = .x$metric,
    Outcome = .x$outcome,
    `AUC (95% CI)` = sprintf("%.3f (%.3f-%.3f)", .x$auc, .x$auc_lower, .x$auc_upper),
    `Optimal Cutoff` = round(.x$cutoff, 4),
    Sensitivity = round(.x$sensitivity, 2),
    Specificity = round(.x$specificity, 2)
  ))

  cat("\n--- Predictive Performance Summary ---\n")
  print(roc_summary, n = 100)
  write_csv(roc_summary, "pk_response_roc_results.csv")
  cat("\nSaved: pk_response_roc_results.csv\n")
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
cat("  - pk_response_statistical_results.csv (continuous OR/HR results)\n")
cat("  - pk_response_roc_results.csv (ROC/AUC predictive performance)\n")
cat("  - boxplot_VGPR_wilcoxon.png (VGPR boxplots with p-values)\n")
cat("  - boxplot_2monthPFS_wilcoxon.png (2-month PFS boxplots with p-values)\n")
cat("  - forest_VGPR_OR.png (VGPR forest plot)\n")
cat("  - forest_2monthPFS_OR.png (2-month PFS forest plot)\n")
cat("  - forest_PFS_HR.png (PFS hazard ratio forest plot)\n")
cat("==========================================================\n")
cat("                 Analysis Complete!\n")
cat("==========================================================\n")
