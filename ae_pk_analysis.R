#===============================================================================
# Teclistamab PK-AE Association Analysis
#
# Compare PK metrics by Adverse Event groups
# Boxplots + Median (IQR) + p-values
#
# PREREQUISITE: Run pk_simulation.R first to generate pk_ae_merged_results.csv
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

library(tidyverse)
library(gridExtra)

cat("==========================================================\n")
cat("    Teclistamab PK-AE Association Analysis\n")
cat("==========================================================\n\n")

#-------------------------------------------------------------------------------
# 1. Load Pre-computed PK Simulation Results
#-------------------------------------------------------------------------------

cat("Loading pre-computed PK simulation results...\n")

# Check if pre-computed results exist
if (!file.exists("pk_ae_merged_results.csv")) {
  stop("ERROR: pk_ae_merged_results.csv not found!\n",
       "Please run pk_simulation.R first to generate PK simulation results.")
}

analysis_data <- read_csv("pk_ae_merged_results.csv", show_col_types = FALSE)
cat("  - Loaded", nrow(analysis_data), "patients from pk_ae_merged_results.csv\n")

#-------------------------------------------------------------------------------
# 2. Define AE Groups and PK Metrics
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
# 3. Statistical Analysis Function
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
# 4. Run All Statistical Comparisons
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
# 5. Logistic Regression (OR) for CRS and Neurotoxicity
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
# 6. Categorical Analysis (Median Split) for CRS and Neurotoxicity
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
# 7. Predictive Performance for CRS and Neurotoxicity
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
# 8. Create Boxplots
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
# 9. Summary Figure for Key Metrics (CRS focus)
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
# 10. Print Final Summary
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
