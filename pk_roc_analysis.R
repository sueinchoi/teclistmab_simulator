#===============================================================================
# Teclistamab PK ROC Analysis: Cavg_120hr Cut-off Evaluation
#
# Generate:
# 1. Cut-off comparison table with performance metrics
# 2. ROC curves for Response (≥VGPR), CRS (any grade), CRS (Grade 2+)
# 3. Optimal cut-off selection based on Youden's Index
#
# PREREQUISITE: Run pk_simulation.R first
#===============================================================================

library(tidyverse)
library(pROC)
library(gridExtra)

cat("==========================================================\n")
cat("    Teclistamab PK ROC Analysis: Cavg_120hr\n")
cat("==========================================================\n\n")

#-------------------------------------------------------------------------------
# 1. Load Data
#-------------------------------------------------------------------------------

cat("Loading data...\n")

pk_ae_data <- read_csv("pk_ae_merged_results.csv", show_col_types = FALSE)
response_data <- read_csv("response_data.csv", show_col_types = FALSE)

dosing_all <- read_csv("mrgsolve_dosing_full.csv", show_col_types = FALSE) %>%
  filter(TIME >= 0) %>%
  arrange(ID, TIME)

dose_counts <- dosing_all %>%
  group_by(ID) %>%
  summarise(N_doses_actual = n(), .groups = "drop")

#-------------------------------------------------------------------------------
# 2. Prepare Analysis Data
#-------------------------------------------------------------------------------

analysis_data <- pk_ae_data %>%
  left_join(response_data, by = "PID") %>%
  left_join(dose_counts, by = "ID") %>%
  mutate(
    VGPR_or_better = ifelse(RESP_CTX %in% c("sCR", "CR", "VGPR"), 1, 0),
    Has_6doses = N_doses_actual >= 6
  )

cat("Total patients:", nrow(analysis_data), "\n")
cat("Patients with ≥6 doses:", sum(analysis_data$Has_6doses, na.rm = TRUE), "\n")

#-------------------------------------------------------------------------------
# 3. Function to Calculate Performance Metrics at Each Cut-off
#-------------------------------------------------------------------------------

calc_cutoff_metrics <- function(data, outcome_var, pk_var, cutoffs) {

  results <- map_dfr(cutoffs, function(cutoff) {
    # Predicted positive if pk_var >= cutoff
    pred_pos <- data[[pk_var]] >= cutoff
    actual_pos <- data[[outcome_var]] == 1

    # Confusion matrix
    TP <- sum(pred_pos & actual_pos, na.rm = TRUE)
    TN <- sum(!pred_pos & !actual_pos, na.rm = TRUE)
    FP <- sum(pred_pos & !actual_pos, na.rm = TRUE)
    FN <- sum(!pred_pos & actual_pos, na.rm = TRUE)

    # Metrics
    sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), NA)
    specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), NA)
    ppv <- ifelse((TP + FP) > 0, TP / (TP + FP), NA)
    npv <- ifelse((TN + FN) > 0, TN / (TN + FN), NA)
    accuracy <- ifelse((TP + TN + FP + FN) > 0, (TP + TN) / (TP + TN + FP + FN), NA)

    # Youden's J Index
    youden_j <- sensitivity + specificity - 1

    # F1 Score
    f1 <- ifelse((sensitivity + ppv) > 0, 2 * (ppv * sensitivity) / (ppv + sensitivity), NA)

    # Likelihood Ratios
    lr_pos <- ifelse(specificity < 1, sensitivity / (1 - specificity), Inf)
    lr_neg <- ifelse(sensitivity < 1, (1 - sensitivity) / specificity, Inf)

    tibble(
      Cutoff = cutoff,
      Sensitivity = round(sensitivity, 3),
      Specificity = round(specificity, 3),
      PPV = round(ppv, 3),
      NPV = round(npv, 3),
      Accuracy = round(accuracy, 3),
      `Youden's J` = round(youden_j, 3),
      `F1 Score` = round(f1, 3),
      `LR+` = round(lr_pos, 3),
      `LR-` = round(lr_neg, 3)
    )
  })

  return(results)
}

#-------------------------------------------------------------------------------
# 4. Response (≥VGPR) Analysis - Patients with ≥6 doses
#-------------------------------------------------------------------------------

cat("\n=== Response (≥VGPR) Cut-off Analysis ===\n")

# Filter for ≥6 doses
response_data_filtered <- analysis_data %>%
  filter(Has_6doses & !is.na(VGPR_or_better) & !is.na(Cavg_120hr))

cat("Patients for response analysis:", nrow(response_data_filtered), "\n")
cat("≥VGPR events:", sum(response_data_filtered$VGPR_or_better), "\n")

# Define cut-offs to evaluate
cutoffs_response <- seq(0.15, 0.45, by = 0.05)

# Calculate metrics
response_metrics <- calc_cutoff_metrics(
  response_data_filtered,
  "VGPR_or_better",
  "Cavg_120hr",
  cutoffs_response
)

# Find optimal cut-off by Youden's J
optimal_response <- response_metrics %>%
  filter(!is.na(`Youden's J`)) %>%
  arrange(desc(`Youden's J`)) %>%
  head(1)

cat("\n--- Response (≥VGPR) Cut-off Comparison ---\n")
print(response_metrics, n = 20)
cat(sprintf("\n>>> Optimal Cut-off (Youden's J): %.2f <<<\n", optimal_response$Cutoff))

# ROC curve
roc_response <- roc(response_data_filtered$VGPR_or_better,
                    response_data_filtered$Cavg_120hr, quiet = TRUE)
auc_response <- auc(roc_response)
ci_response <- ci.auc(roc_response)

# Get optimal threshold from ROC
coords_response <- coords(roc_response, "best", ret = c("threshold", "sensitivity", "specificity"))

cat(sprintf("ROC AUC: %.3f (%.3f-%.3f)\n", auc_response, ci_response[1], ci_response[3]))
cat(sprintf("ROC Optimal Threshold: %.3f\n", coords_response$threshold))

#-------------------------------------------------------------------------------
# 5. CRS (any grade) Analysis - All patients
#-------------------------------------------------------------------------------

cat("\n=== CRS (any grade) Cut-off Analysis ===\n")

crs_data <- analysis_data %>%
  filter(!is.na(CRS_any) & !is.na(Cavg_120hr))

cat("Patients for CRS analysis:", nrow(crs_data), "\n")
cat("CRS events:", sum(crs_data$CRS_any), "\n")

# Define cut-offs
cutoffs_crs <- seq(0.15, 0.55, by = 0.05)

# Calculate metrics
crs_metrics <- calc_cutoff_metrics(crs_data, "CRS_any", "Cavg_120hr", cutoffs_crs)

# Find optimal cut-off
optimal_crs <- crs_metrics %>%
  filter(!is.na(`Youden's J`)) %>%
  arrange(desc(`Youden's J`)) %>%
  head(1)

cat("\n--- CRS (any grade) Cut-off Comparison ---\n")
print(crs_metrics, n = 20)
cat(sprintf("\n>>> Optimal Cut-off (Youden's J): %.2f <<<\n", optimal_crs$Cutoff))

# ROC curve
roc_crs <- roc(crs_data$CRS_any, crs_data$Cavg_120hr, quiet = TRUE)
auc_crs <- auc(roc_crs)
ci_crs <- ci.auc(roc_crs)
coords_crs <- coords(roc_crs, "best", ret = c("threshold", "sensitivity", "specificity"))

cat(sprintf("ROC AUC: %.3f (%.3f-%.3f)\n", auc_crs, ci_crs[1], ci_crs[3]))
cat(sprintf("ROC Optimal Threshold: %.3f\n", coords_crs$threshold))

#-------------------------------------------------------------------------------
# 6. CRS (Grade 2+) Analysis
#-------------------------------------------------------------------------------

cat("\n=== CRS (Grade 2+) Cut-off Analysis ===\n")

crs_gr2_data <- analysis_data %>%
  filter(!is.na(CRS_gr2) & !is.na(Cavg_120hr))

cat("Patients for CRS Gr2+ analysis:", nrow(crs_gr2_data), "\n")
cat("CRS Gr2+ events:", sum(crs_gr2_data$CRS_gr2), "\n")

# Calculate metrics
crs_gr2_metrics <- calc_cutoff_metrics(crs_gr2_data, "CRS_gr2", "Cavg_120hr", cutoffs_crs)

# Find optimal cut-off
optimal_crs_gr2 <- crs_gr2_metrics %>%
  filter(!is.na(`Youden's J`)) %>%
  arrange(desc(`Youden's J`)) %>%
  head(1)

cat("\n--- CRS (Grade 2+) Cut-off Comparison ---\n")
print(crs_gr2_metrics, n = 20)
cat(sprintf("\n>>> Optimal Cut-off (Youden's J): %.2f <<<\n", optimal_crs_gr2$Cutoff))

# ROC curve
roc_crs_gr2 <- roc(crs_gr2_data$CRS_gr2, crs_gr2_data$Cavg_120hr, quiet = TRUE)
auc_crs_gr2 <- auc(roc_crs_gr2)
ci_crs_gr2 <- ci.auc(roc_crs_gr2)
coords_crs_gr2 <- coords(roc_crs_gr2, "best", ret = c("threshold", "sensitivity", "specificity"))

cat(sprintf("ROC AUC: %.3f (%.3f-%.3f)\n", auc_crs_gr2, ci_crs_gr2[1], ci_crs_gr2[3]))
cat(sprintf("ROC Optimal Threshold: %.3f\n", coords_crs_gr2$threshold))

#-------------------------------------------------------------------------------
# 7. Create ROC Curve Plots
#-------------------------------------------------------------------------------

cat("\n=== Creating ROC Curves ===\n")

# Function to create ROC plot
create_roc_plot <- function(roc_obj, auc_val, ci_vals, optimal_coords, title, color) {

  # Extract ROC data
  roc_df <- data.frame(
    specificity = roc_obj$specificities,
    sensitivity = roc_obj$sensitivities
  )

  # AUC label
  auc_label <- sprintf("AUC = %.3f (%.3f-%.3f)", auc_val, ci_vals[1], ci_vals[3])
  optimal_label <- sprintf("Optimal: %.3f", optimal_coords$threshold)

  ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
    geom_line(color = color, size = 1.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    geom_point(
      data = data.frame(
        x = 1 - optimal_coords$specificity,
        y = optimal_coords$sensitivity
      ),
      aes(x = x, y = y),
      color = "red", size = 4
    ) +
    annotate("text", x = 0.6, y = 0.2, label = auc_label, size = 3.5, hjust = 0) +
    annotate("point", x = 0.55, y = 0.1, color = "red", size = 3) +
    annotate("text", x = 0.6, y = 0.1, label = optimal_label, size = 3.5, hjust = 0) +
    labs(
      title = title,
      x = "1 - Specificity (FPR)",
      y = "Sensitivity (TPR)"
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      panel.grid.minor = element_blank()
    )
}

# Create individual ROC plots
p_roc_response <- create_roc_plot(
  roc_response, auc_response, ci_response, coords_response,
  "A. Response (≥VGPR)", "#2ecc71"
)

p_roc_crs <- create_roc_plot(
  roc_crs, auc_crs, ci_crs, coords_crs,
  "B. CRS (any grade)", "#3498db"
)

p_roc_crs_gr2 <- create_roc_plot(
  roc_crs_gr2, auc_crs_gr2, ci_crs_gr2, coords_crs_gr2,
  "C. CRS (Grade 2+)", "#e74c3c"
)

# Combine ROC plots
roc_combined <- grid.arrange(
  p_roc_response, p_roc_crs, p_roc_crs_gr2,
  ncol = 3, nrow = 1
)

ggsave("roc_curves_combined.png", roc_combined, width = 14, height = 5, dpi = 300)
cat("Saved: roc_curves_combined.png\n")

# Save individual plots
ggsave("roc_response.png", p_roc_response, width = 5, height = 5, dpi = 300)
ggsave("roc_crs_any.png", p_roc_crs, width = 5, height = 5, dpi = 300)
ggsave("roc_crs_gr2.png", p_roc_crs_gr2, width = 5, height = 5, dpi = 300)

#-------------------------------------------------------------------------------
# 8. Save Cut-off Comparison Tables
#-------------------------------------------------------------------------------

cat("\n=== Saving Tables ===\n")

# Add outcome column and combine
response_metrics_out <- response_metrics %>%
  mutate(Outcome = "Response (≥VGPR)", .before = 1)

crs_metrics_out <- crs_metrics %>%
  mutate(Outcome = "CRS (any grade)", .before = 1)

crs_gr2_metrics_out <- crs_gr2_metrics %>%
  mutate(Outcome = "CRS (Grade 2+)", .before = 1)

# Combined table
all_metrics <- bind_rows(response_metrics_out, crs_metrics_out, crs_gr2_metrics_out)

write_csv(all_metrics, "cutoff_metrics_all.csv")
cat("Saved: cutoff_metrics_all.csv\n")

# Individual tables
write_csv(response_metrics, "cutoff_metrics_response.csv")
write_csv(crs_metrics, "cutoff_metrics_crs_any.csv")
write_csv(crs_gr2_metrics, "cutoff_metrics_crs_gr2.csv")
cat("Saved: cutoff_metrics_response.csv\n")
cat("Saved: cutoff_metrics_crs_any.csv\n")
cat("Saved: cutoff_metrics_crs_gr2.csv\n")

#-------------------------------------------------------------------------------
# 9. Summary Table (Optimal Cut-offs)
#-------------------------------------------------------------------------------

cat("\n=== Optimal Cut-off Summary ===\n")

optimal_summary <- bind_rows(
  optimal_response %>% mutate(Outcome = "Response (≥VGPR)",
                               AUC = sprintf("%.3f (%.3f-%.3f)", auc_response, ci_response[1], ci_response[3])),
  optimal_crs %>% mutate(Outcome = "CRS (any grade)",
                          AUC = sprintf("%.3f (%.3f-%.3f)", auc_crs, ci_crs[1], ci_crs[3])),
  optimal_crs_gr2 %>% mutate(Outcome = "CRS (Grade 2+)",
                              AUC = sprintf("%.3f (%.3f-%.3f)", auc_crs_gr2, ci_crs_gr2[1], ci_crs_gr2[3]))
) %>%
  select(Outcome, Cutoff, AUC, Sensitivity, Specificity, `Youden's J`, PPV, NPV, Accuracy, `F1 Score`, `LR+`, `LR-`)

cat("\n")
print(optimal_summary, width = Inf)

write_csv(optimal_summary, "optimal_cutoff_summary.csv")
cat("\nSaved: optimal_cutoff_summary.csv\n")

#-------------------------------------------------------------------------------
# 10. Print Final Summary
#-------------------------------------------------------------------------------

cat("\n")
cat("==========================================================\n")
cat("              OPTIMAL CUT-OFF SUMMARY\n")
cat("==========================================================\n\n")

cat(sprintf("Response (≥VGPR):\n"))
cat(sprintf("  Optimal Cut-off: %.2f μg/mL\n", optimal_response$Cutoff))
cat(sprintf("  AUC: %.3f (%.3f-%.3f)\n", auc_response, ci_response[1], ci_response[3]))
cat(sprintf("  Sensitivity: %.3f, Specificity: %.3f\n", optimal_response$Sensitivity, optimal_response$Specificity))
cat(sprintf("  Youden's J: %.3f\n\n", optimal_response$`Youden's J`))

cat(sprintf("CRS (any grade):\n"))
cat(sprintf("  Optimal Cut-off: %.2f μg/mL\n", optimal_crs$Cutoff))
cat(sprintf("  AUC: %.3f (%.3f-%.3f)\n", auc_crs, ci_crs[1], ci_crs[3]))
cat(sprintf("  Sensitivity: %.3f, Specificity: %.3f\n", optimal_crs$Sensitivity, optimal_crs$Specificity))
cat(sprintf("  Youden's J: %.3f\n\n", optimal_crs$`Youden's J`))

cat(sprintf("CRS (Grade 2+):\n"))
cat(sprintf("  Optimal Cut-off: %.2f μg/mL\n", optimal_crs_gr2$Cutoff))
cat(sprintf("  AUC: %.3f (%.3f-%.3f)\n", auc_crs_gr2, ci_crs_gr2[1], ci_crs_gr2[3]))
cat(sprintf("  Sensitivity: %.3f, Specificity: %.3f\n", optimal_crs_gr2$Sensitivity, optimal_crs_gr2$Specificity))
cat(sprintf("  Youden's J: %.3f\n", optimal_crs_gr2$`Youden's J`))

cat("\n==========================================================\n")
cat("                    OUTPUT FILES\n")
cat("==========================================================\n")
cat("  Tables:\n")
cat("    - cutoff_metrics_all.csv (all outcomes combined)\n")
cat("    - cutoff_metrics_response.csv\n")
cat("    - cutoff_metrics_crs_any.csv\n")
cat("    - cutoff_metrics_crs_gr2.csv\n")
cat("    - optimal_cutoff_summary.csv\n")
cat("  Figures:\n")
cat("    - roc_curves_combined.png (3-panel ROC)\n")
cat("    - roc_response.png\n")
cat("    - roc_crs_any.png\n")
cat("    - roc_crs_gr2.png\n")
cat("==========================================================\n")
cat("                 Analysis Complete!\n")
cat("==========================================================\n")
