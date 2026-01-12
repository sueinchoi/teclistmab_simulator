#===============================================================================
# Teclistamab PK-Response Association Analysis
#
# For patients with ≥4 weeks of treatment:
# 1. VGPR or better (≥VGPR) - Odds Ratio (Logistic Regression)
# 2. 2-month PFS status - Odds Ratio (Logistic Regression)
# 3. PFS - Hazard Ratio (Cox Regression)
#===============================================================================

library(tidyverse)
library(survival)
library(broom)
library(gridExtra)

set.seed(12345)

cat("==========================================================\n")
cat("    Teclistamab PK-Response Association Analysis\n")
cat("==========================================================\n\n")

#-------------------------------------------------------------------------------
# 1. Load Data
#-------------------------------------------------------------------------------

cat("Loading data...\n")

# PK metrics data (from previous analysis)
pk_data <- read_csv("pk_ae_merged_results.csv", show_col_types = FALSE)
cat("  - Loaded PK data for", nrow(pk_data), "patients\n")

# Response data
response_data <- read_csv("response_data.csv", show_col_types = FALSE)
cat("  - Loaded response data for", nrow(response_data), "patients\n")

# Merge data
analysis_data <- pk_data %>%
  left_join(response_data, by = "PID") %>%
  mutate(
    # VGPR or better (sCR, CR, VGPR)
    VGPR_or_better = ifelse(RESP_CTX %in% c("sCR", "CR", "VGPR"), 1, 0),

    # 2-month (60 days) PFS status: 1 = alive without progression at 2 months
    PFS_2month = ifelse(DAYS_VS_PFS_CTX >= 60 | (VS_PFS_CTX == 0 & DAYS_VS_PFS_CTX < 60),
                        ifelse(DAYS_VS_PFS_CTX >= 60 | VS_PFS_CTX == 0, 1, 0), 0),

    # For Cox model: event and time
    PFS_event = VS_PFS_CTX,  # 1 = event (progression/death), 0 = censored
    PFS_time = DAYS_VS_PFS_CTX,

    # Treatment duration (N_TEC_DOSES as proxy for weeks)
    Treatment_weeks = N_TEC_DOSES  # Approximately 1 dose per week
  )

# Correct 2-month PFS calculation
analysis_data <- analysis_data %>%
  mutate(
    PFS_2month = case_when(
      DAYS_VS_PFS_CTX >= 60 ~ 1,  # Survived 2 months without event or with event after 2 months
      VS_PFS_CTX == 0 ~ 1,        # Censored before 60 days but no event
      VS_PFS_CTX == 1 & DAYS_VS_PFS_CTX < 60 ~ 0,  # Event before 60 days
      TRUE ~ NA_real_
    )
  )

cat("\n=== Data Summary ===\n")
cat("Total patients:", nrow(analysis_data), "\n")

#-------------------------------------------------------------------------------
# 2. Filter: Patients with ≥4 weeks of treatment
#-------------------------------------------------------------------------------

cat("\n=== Filtering for ≥4 weeks treatment ===\n")

# Filter patients with at least 4 doses (approximately 4 weeks)
filtered_data <- analysis_data %>%
  filter(N_doses >= 4 | Treatment_weeks >= 4)

cat("Patients with ≥4 weeks treatment:", nrow(filtered_data), "\n")

# Check if we have enough data
if (nrow(filtered_data) < 10) {
  cat("\nNote: Using N_doses >= 3 as alternative criteria\n")
  filtered_data <- analysis_data %>%
    filter(N_doses >= 3)
  cat("Patients with ≥3 doses:", nrow(filtered_data), "\n")
}

# Summary of response outcomes
cat("\n=== Response Summary (filtered patients) ===\n")
cat("Response distribution:\n")
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
# 3. Define PK Metrics for Analysis
#-------------------------------------------------------------------------------

pk_metrics <- c("Cmax_72hr", "Cavg_72hr", "Cmax_120hr", "Cavg_120hr",
                "Cmax_dose1", "Cavg_dose1", "Cmax_dose3", "Cavg_dose3")

#-------------------------------------------------------------------------------
# 4. Logistic Regression Function (for OR)
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
# 5. Cox Regression Function (for HR)
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
# 6. Run Analyses
#-------------------------------------------------------------------------------

cat("\n")
cat("==========================================================\n")
cat("             STATISTICAL ANALYSIS RESULTS\n")
cat("==========================================================\n")

# 6.1 VGPR or better (Logistic Regression - OR)
cat("\n\n--- 1. VGPR or Better (≥VGPR) - Odds Ratio ---\n")

vgpr_results <- map_dfr(pk_metrics, function(pk_var) {
  run_logistic_analysis(filtered_data, "VGPR_or_better", pk_var, "≥VGPR")
})

if (nrow(vgpr_results) > 0) {
  print(vgpr_results %>% select(-OR, -OR_lower, -OR_upper), n = 100)
} else {
  cat("Insufficient data for VGPR analysis\n")
}

# 6.2 2-month PFS (Logistic Regression - OR)
cat("\n\n--- 2. 2-Month PFS Status - Odds Ratio ---\n")

pfs2m_results <- map_dfr(pk_metrics, function(pk_var) {
  run_logistic_analysis(filtered_data, "PFS_2month", pk_var, "2-month PFS")
})

if (nrow(pfs2m_results) > 0) {
  print(pfs2m_results %>% select(-OR, -OR_lower, -OR_upper), n = 100)
} else {
  cat("Insufficient data for 2-month PFS analysis\n")
}

# 6.3 PFS (Cox Regression - HR)
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
# 7. Combine and Save Results
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
# 8. Create Forest Plots
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
# 9. Summary Boxplots by Response
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
# 10. Final Summary
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
