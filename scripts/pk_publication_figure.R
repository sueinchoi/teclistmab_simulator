#===============================================================================
# Teclistamab PK Publication Figure
#
# 4-Panel Figure:
# A: Cavg_120hr by CRS Grade and ICANS Grade (boxplots)
# B: Cavg_120hr by ≥VGPR and 2-month PFS (boxplots)
# C: Therapeutic window (Event Rate vs Cut-off)
# D: Event rates by exposure category (bar plot)
#
# PREREQUISITE: Run pk_simulation.R first
#===============================================================================

library(tidyverse)
library(gridExtra)
library(grid)
library(ggpubr)

cat("==========================================================\n")
cat("    Teclistamab PK Publication Figure\n")
cat("==========================================================\n\n")

#-------------------------------------------------------------------------------
# 1. Load Data
#-------------------------------------------------------------------------------

cat("Loading data...\n")

# Load PK + AE merged data
pk_ae_data <- read_csv("output/tables/pk_ae_merged_results.csv", show_col_types = FALSE)
cat("  - Loaded PK+AE data for", nrow(pk_ae_data), "patients\n")

# Load response data
response_data <- read_csv("data/response_data.csv", show_col_types = FALSE)
cat("  - Loaded response data\n")

# Load dosing data
dosing_all <- read_csv("output/tables/mrgsolve_dosing_full.csv", show_col_types = FALSE) %>%
  filter(TIME >= 0) %>%
  arrange(ID, TIME)

# Calculate number of doses per patient
dose_counts <- dosing_all %>%
  group_by(ID) %>%
  summarise(N_doses_actual = n(), .groups = "drop")

#-------------------------------------------------------------------------------
# 2. Prepare Analysis Data
#-------------------------------------------------------------------------------

cat("\n=== Preparing data ===\n")

# Merge all data
analysis_data <- pk_ae_data %>%
  left_join(response_data, by = "PID") %>%
  left_join(dose_counts, by = "ID") %>%
  mutate(
    # CRS Grade (0, 1, ≥2)
    CRS_Grade = case_when(
      CRS_any == 0 ~ "0",
      CRS_any == 1 & CRS_gr2 == 0 ~ "1",
      CRS_gr2 == 1 ~ "≥2"
    ),
    CRS_Grade = factor(CRS_Grade, levels = c("0", "1", "≥2")),

    # ICANS Grade (0, 1, ≥2)
    ICANS_Grade = case_when(
      Neuro_any == 0 ~ "0",
      Neuro_any == 1 & Neuro_gr2 == 0 ~ "1",
      Neuro_gr2 == 1 ~ "≥2"
    ),
    ICANS_Grade = factor(ICANS_Grade, levels = c("0", "1", "≥2")),

    # Response (≥VGPR vs <VGPR)
    Response_VGPR = case_when(
      RESP_CTX %in% c("sCR", "CR", "VGPR") ~ "≥VGPR",
      TRUE ~ "<VGPR"
    ),
    Response_VGPR = factor(Response_VGPR, levels = c("≥VGPR", "<VGPR")),

    # 2-month PFS (Progressed vs No progression)
    PFS_2month_status = case_when(
      VS_PFS_CTX == 1 & DAYS_VS_PFS_CTX <= 60 ~ "Progressed",
      TRUE ~ "No progression"
    ),
    PFS_2month_status = factor(PFS_2month_status, levels = c("Progressed", "No progression")),

    # At least 6 doses
    Has_6doses = N_doses_actual >= 6,

    # Binary outcomes for rate calculation
    VGPR_or_better = ifelse(RESP_CTX %in% c("sCR", "CR", "VGPR"), 1, 0)
  )

# Summary
cat("CRS Grade distribution:\n")
print(table(analysis_data$CRS_Grade, useNA = "ifany"))
cat("\nICANS Grade distribution:\n")
print(table(analysis_data$ICANS_Grade, useNA = "ifany"))

#-------------------------------------------------------------------------------
# 3. Panel A: CRS and ICANS by Grade
#-------------------------------------------------------------------------------

cat("\n=== Creating Panel A ===\n")

# Color palette
grade_colors <- c("0" = "#2ecc71", "1" = "#f39c12", "≥2" = "#e74c3c")

# CRS boxplot
crs_plot_data <- analysis_data %>%
  filter(!is.na(CRS_Grade) & !is.na(Cavg_120hr))

# Calculate p-values for CRS
crs_pval <- tryCatch({
  kruskal.test(Cavg_120hr ~ CRS_Grade, data = crs_plot_data)$p.value
}, error = function(e) NA)

crs_pval_text <- ifelse(!is.na(crs_pval) & crs_pval < 0.05,
                        sprintf("p = %.3f*", crs_pval),
                        sprintf("p = %.3f", crs_pval))

p_crs <- ggplot(crs_plot_data, aes(x = CRS_Grade, y = Cavg_120hr, fill = CRS_Grade)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 21, outlier.size = 2) +
  scale_fill_manual(values = grade_colors) +
  labs(x = "CRS Grade", y = "Cavg Day 5 (μg/mL)") +
  annotate("text", x = 2, y = max(crs_plot_data$Cavg_120hr, na.rm = TRUE) * 1.1,
           label = crs_pval_text, size = 3.5,
           color = ifelse(!is.na(crs_pval) & crs_pval < 0.05, "red", "black")) +
  scale_y_continuous(limits = c(0, max(crs_plot_data$Cavg_120hr, na.rm = TRUE) * 1.2)) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 10),
    plot.margin = margin(5, 5, 5, 5)
  )

# Add sample sizes to x-axis labels
crs_n <- crs_plot_data %>% count(CRS_Grade)
p_crs <- p_crs +
  scale_x_discrete(labels = function(x) {
    n_vals <- crs_n$n[match(x, crs_n$CRS_Grade)]
    paste0(x, "\n(n=", n_vals, ")")
  })

# ICANS boxplot
icans_plot_data <- analysis_data %>%
  filter(!is.na(ICANS_Grade) & !is.na(Cavg_120hr))

# Calculate p-values for ICANS
icans_pval <- tryCatch({
  kruskal.test(Cavg_120hr ~ ICANS_Grade, data = icans_plot_data)$p.value
}, error = function(e) NA)

icans_pval_text <- ifelse(!is.na(icans_pval) & icans_pval < 0.05,
                          sprintf("p = %.3f*", icans_pval),
                          sprintf("p = %.3f", icans_pval))

p_icans <- ggplot(icans_plot_data, aes(x = ICANS_Grade, y = Cavg_120hr, fill = ICANS_Grade)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 21, outlier.size = 2) +
  scale_fill_manual(values = grade_colors) +
  labs(x = "ICANS Grade", y = "") +
  annotate("text", x = 2, y = max(icans_plot_data$Cavg_120hr, na.rm = TRUE) * 1.1,
           label = icans_pval_text, size = 3.5,
           color = ifelse(!is.na(icans_pval) & icans_pval < 0.05, "red", "black")) +
  scale_y_continuous(limits = c(0, max(icans_plot_data$Cavg_120hr, na.rm = TRUE) * 1.2)) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 10),
    plot.margin = margin(5, 5, 5, 5)
  )

# Add sample sizes
icans_n <- icans_plot_data %>% count(ICANS_Grade)
p_icans <- p_icans +
  scale_x_discrete(labels = function(x) {
    n_vals <- icans_n$n[match(x, icans_n$ICANS_Grade)]
    paste0(x, "\n(n=", n_vals, ")")
  })

# Combine Panel A
panel_A <- ggarrange(p_crs, p_icans, ncol = 2, nrow = 1, widths = c(1, 1))

#-------------------------------------------------------------------------------
# 4. Panel B: Response by ≥VGPR and 2-month PFS
#-------------------------------------------------------------------------------

cat("\n=== Creating Panel B ===\n")

# Response colors
response_colors <- c("≥VGPR" = "#2ecc71", "<VGPR" = "#e74c3c")
pfs_colors <- c("Progressed" = "#e74c3c", "No progression" = "#2ecc71")

# VGPR boxplot (≥6 doses only)
vgpr_plot_data <- analysis_data %>%
  filter(Has_6doses & !is.na(Response_VGPR) & !is.na(Cavg_120hr))

# Calculate p-value
vgpr_pval <- tryCatch({
  wilcox.test(Cavg_120hr ~ Response_VGPR, data = vgpr_plot_data)$p.value
}, error = function(e) NA)

vgpr_pval_text <- ifelse(!is.na(vgpr_pval) & vgpr_pval < 0.05,
                         sprintf("p = %.3f*", vgpr_pval),
                         sprintf("p = %.3f", vgpr_pval))

p_vgpr <- ggplot(vgpr_plot_data, aes(x = Response_VGPR, y = Cavg_120hr, fill = Response_VGPR)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 21, outlier.size = 2) +
  scale_fill_manual(values = response_colors) +
  labs(x = "Response", y = "Cavg Day 5 (μg/mL)") +
  annotate("text", x = 1.5, y = max(vgpr_plot_data$Cavg_120hr, na.rm = TRUE) * 1.1,
           label = vgpr_pval_text, size = 3.5,
           color = ifelse(!is.na(vgpr_pval) & vgpr_pval < 0.05, "red", "black")) +
  scale_y_continuous(limits = c(0, max(vgpr_plot_data$Cavg_120hr, na.rm = TRUE) * 1.2)) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 10),
    plot.margin = margin(5, 5, 5, 5)
  )

# Add sample sizes
vgpr_n <- vgpr_plot_data %>% count(Response_VGPR)
p_vgpr <- p_vgpr +
  scale_x_discrete(labels = function(x) {
    n_vals <- vgpr_n$n[match(x, vgpr_n$Response_VGPR)]
    paste0(x, "\n(n=", n_vals, ")")
  })

# 2-month PFS boxplot (≥6 doses only)
pfs_plot_data <- analysis_data %>%
  filter(Has_6doses & !is.na(PFS_2month_status) & !is.na(Cavg_120hr))

# Calculate p-value
pfs_pval <- tryCatch({
  wilcox.test(Cavg_120hr ~ PFS_2month_status, data = pfs_plot_data)$p.value
}, error = function(e) NA)

pfs_pval_text <- ifelse(!is.na(pfs_pval) & pfs_pval < 0.05,
                        sprintf("p = %.3f*", pfs_pval),
                        sprintf("p = %.3f", pfs_pval))

p_pfs <- ggplot(pfs_plot_data, aes(x = PFS_2month_status, y = Cavg_120hr, fill = PFS_2month_status)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 21, outlier.size = 2) +
  scale_fill_manual(values = pfs_colors) +
  labs(x = "2 months PFS", y = "") +
  annotate("text", x = 1.5, y = max(pfs_plot_data$Cavg_120hr, na.rm = TRUE) * 1.1,
           label = pfs_pval_text, size = 3.5,
           color = ifelse(!is.na(pfs_pval) & pfs_pval < 0.05, "red", "black")) +
  scale_y_continuous(limits = c(0, max(pfs_plot_data$Cavg_120hr, na.rm = TRUE) * 1.2)) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 10),
    plot.margin = margin(5, 5, 5, 5)
  )

# Add sample sizes
pfs_n <- pfs_plot_data %>% count(PFS_2month_status)
p_pfs <- p_pfs +
  scale_x_discrete(labels = function(x) {
    n_vals <- pfs_n$n[match(x, pfs_n$PFS_2month_status)]
    paste0(x, "\n(n=", n_vals, ")")
  })

# Combine Panel B
panel_B <- ggarrange(p_vgpr, p_pfs, ncol = 2, nrow = 1, widths = c(1, 1))

#-------------------------------------------------------------------------------
# 5. Panel C: Therapeutic Window (Event Rate vs Cut-off)
#    Using ROC-based optimal cut-offs
#-------------------------------------------------------------------------------

cat("\n=== Creating Panel C ===\n")

library(pROC)

# Define cut-off range
cavg_values <- analysis_data$Cavg_120hr
cutoff_min <- quantile(cavg_values, 0.05, na.rm = TRUE)
cutoff_max <- quantile(cavg_values, 0.95, na.rm = TRUE)
cutoffs <- seq(cutoff_min, cutoff_max, length.out = 50)

# Calculate event rates at each cut-off
event_rate_results <- map_dfr(cutoffs, function(cutoff) {
  # CRS analysis (all patients)
  crs_high <- analysis_data %>% filter(Cavg_120hr >= cutoff)

  crs_rate_high <- mean(crs_high$CRS_any, na.rm = TRUE)

  # VGPR analysis (≥6 doses only)
  vgpr_data <- analysis_data %>% filter(Has_6doses == TRUE)
  vgpr_high <- vgpr_data %>% filter(Cavg_120hr >= cutoff)

  vgpr_rate_high <- mean(vgpr_high$VGPR_or_better, na.rm = TRUE)

  tibble(
    Cutoff = cutoff,
    CRS_rate = crs_rate_high * 100,
    VGPR_rate = vgpr_rate_high * 100
  )
})

#--- ROC-based optimal cut-offs (matching pk_roc_analysis.R) ---
cat("\n--- Finding ROC-based optimal cut-offs ---\n")

# 1. Response (≥VGPR): Use Youden's Index (≥6 doses only) - LOWER BOUND
response_data_roc <- analysis_data %>%
  filter(Has_6doses & !is.na(VGPR_or_better) & !is.na(Cavg_120hr))

roc_response <- roc(response_data_roc$VGPR_or_better,
                    response_data_roc$Cavg_120hr, quiet = TRUE)
coords_response <- coords(roc_response, "best", ret = c("threshold", "sensitivity", "specificity"),
                          best.method = "youden")

optimal_lower <- coords_response$threshold
cat(sprintf("Response cut-off (Youden): %.4f\n", optimal_lower))

# 2. CRS (Grade 2+): Use Accuracy (all patients) - UPPER BOUND
crs_gr2_data_roc <- analysis_data %>%
  filter(!is.na(CRS_gr2) & !is.na(Cavg_120hr))

# Calculate metrics at different cut-offs to find best Accuracy for CRS Gr2+
cutoffs_crs_gr2 <- seq(0.15, 0.55, by = 0.05)
crs_gr2_metrics <- map_dfr(cutoffs_crs_gr2, function(cutoff) {
  pred_pos <- crs_gr2_data_roc$Cavg_120hr >= cutoff
  actual_pos <- crs_gr2_data_roc$CRS_gr2 == 1

  TP <- sum(pred_pos & actual_pos, na.rm = TRUE)
  TN <- sum(!pred_pos & !actual_pos, na.rm = TRUE)
  FP <- sum(pred_pos & !actual_pos, na.rm = TRUE)
  FN <- sum(!pred_pos & actual_pos, na.rm = TRUE)

  accuracy <- ifelse((TP + TN + FP + FN) > 0, (TP + TN) / (TP + TN + FP + FN), NA)

  tibble(Cutoff = cutoff, Accuracy = accuracy)
})

optimal_crs_gr2 <- crs_gr2_metrics %>%
  filter(!is.na(Accuracy)) %>%
  arrange(desc(Accuracy)) %>%
  head(1)

optimal_upper <- optimal_crs_gr2$Cutoff
cat(sprintf("CRS Grade 2+ cut-off (Accuracy): %.4f\n", optimal_upper))

# Ensure valid range
if (optimal_lower >= optimal_upper) {
  cat("Warning: Lower >= Upper, swapping...\n")
  temp <- optimal_lower
  optimal_lower <- optimal_upper
  optimal_upper <- temp
}

cat(sprintf("Final Optimal Range: %.4f - %.4f\n", optimal_lower, optimal_upper))

# Create long format for plotting
plot_data_c <- event_rate_results %>%
  pivot_longer(cols = c(VGPR_rate, CRS_rate), names_to = "Outcome", values_to = "Rate") %>%
  mutate(
    Outcome = case_when(
      Outcome == "VGPR_rate" ~ "Response Rate (≥VGPR)",
      Outcome == "CRS_rate" ~ "CRS Rate (any Grade)"
    )
  )

panel_C <- ggplot(plot_data_c, aes(x = Cutoff, y = Rate, color = Outcome)) +
  annotate("rect", xmin = optimal_lower, xmax = optimal_upper,
           ymin = -Inf, ymax = Inf, fill = "#2ecc71", alpha = 0.2) +
  geom_line(size = 1.2) +
  geom_point(size = 1.5) +
  geom_vline(xintercept = optimal_lower, linetype = "dashed", color = "darkgreen", size = 0.7) +
  geom_vline(xintercept = optimal_upper, linetype = "dashed", color = "darkgreen", size = 0.7) +
  scale_color_manual(values = c("Response Rate (≥VGPR)" = "#2ecc71", "CRS Rate (any Grade)" = "#e74c3c")) +
  annotate("text", x = (optimal_lower + optimal_upper) / 2, y = 10,
           label = "Optimal Range", color = "darkgreen", fontface = "bold", size = 3.5) +
  # Add cut-off value annotations
  annotate("text", x = optimal_lower, y = 100,
           label = sprintf("%.3f", optimal_lower),
           color = "darkgreen", size = 3, hjust = -0.1) +
  annotate("text", x = optimal_upper, y = 100,
           label = sprintf("%.3f", optimal_upper),
           color = "darkgreen", size = 3, hjust = 1.1) +
  labs(
    x = "Cavg Day 5 Cutoff (μg/mL)",
    y = "Event Rate (%)",
    color = ""
  ) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 20)) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    plot.margin = margin(5, 10, 5, 10)
  )

#-------------------------------------------------------------------------------
# 6. Panel D: Event Rates by Exposure Category (Bar Plot)
#-------------------------------------------------------------------------------

cat("\n=== Creating Panel D ===\n")

# Categorize patients
analysis_data <- analysis_data %>%
  mutate(
    Exposure_Category = case_when(
      Cavg_120hr < optimal_lower ~ "Low",
      Cavg_120hr >= optimal_lower & Cavg_120hr <= optimal_upper ~ "Medium",
      Cavg_120hr > optimal_upper ~ "High"
    ),
    Exposure_Category = factor(Exposure_Category, levels = c("Low", "Medium", "High"))
  )

# Calculate event rates by category
crs_by_cat <- analysis_data %>%
  group_by(Exposure_Category) %>%
  summarise(N = n(), Rate = mean(CRS_any, na.rm = TRUE) * 100, .groups = "drop") %>%
  mutate(Outcome = "CRS (any Grade)")

vgpr_by_cat <- analysis_data %>%
  filter(Has_6doses) %>%
  group_by(Exposure_Category) %>%
  summarise(N = n(), Rate = mean(VGPR_or_better, na.rm = TRUE) * 100, .groups = "drop") %>%
  mutate(Outcome = "Response (≥VGPR)")

bar_data <- bind_rows(vgpr_by_cat, crs_by_cat)

# Create bar labels
bar_data <- bar_data %>%
  mutate(
    label = sprintf("%.0f%%", Rate),
    x_label = paste0(Exposure_Category, "\n(n=", N, ")")
  )

panel_D <- ggplot(bar_data, aes(x = Exposure_Category, y = Rate, fill = Outcome)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = label), position = position_dodge(width = 0.8),
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("Response (≥VGPR)" = "#2ecc71", "CRS (any Grade)" = "#e74c3c")) +
  labs(
    x = "Cavg Day 5 Group",
    y = "Event Rate (%)",
    fill = ""
  ) +
  scale_y_continuous(limits = c(0, 115), breaks = seq(0, 100, 20)) +
  scale_x_discrete(labels = function(x) {
    sapply(x, function(cat) {
      n_vgpr <- vgpr_by_cat$N[vgpr_by_cat$Exposure_Category == cat]
      range_text <- case_when(
        cat == "Low" ~ sprintf("(<%.2f)", optimal_lower),
        cat == "Medium" ~ sprintf("(%.2f-%.2f)", optimal_lower, optimal_upper),
        cat == "High" ~ sprintf("(>%.2f)", optimal_upper)
      )
      paste0(cat, "\n", range_text)
    })
  }) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    plot.margin = margin(5, 10, 5, 10)
  )

#-------------------------------------------------------------------------------
# 7. Combine All Panels
#-------------------------------------------------------------------------------

cat("\n=== Combining panels ===\n")

# Add panel labels to individual plots using labs(tag) for top-left positioning
# Panel A (CRS and ICANS combined) - add label to first subplot
p_crs <- p_crs + labs(tag = "A") +
  theme(plot.tag = element_text(face = "bold", size = 14),
        plot.tag.position = c(0, 1))

p_vgpr <- p_vgpr + labs(tag = "B") +
  theme(plot.tag = element_text(face = "bold", size = 14),
        plot.tag.position = c(0, 1))

panel_C <- panel_C + labs(tag = "C") +
  theme(plot.tag = element_text(face = "bold", size = 14),
        plot.tag.position = c(0, 1))

panel_D <- panel_D + labs(tag = "D") +
  theme(plot.tag = element_text(face = "bold", size = 14),
        plot.tag.position = c(0, 1))

# Recreate panel A and B with updated plots
panel_A <- ggarrange(p_crs, p_icans, ncol = 2, nrow = 1, widths = c(1, 1))
panel_B <- ggarrange(p_vgpr, p_pfs, ncol = 2, nrow = 1, widths = c(1, 1))

# Combine into 2x2 layout
final_figure <- ggarrange(
  panel_A, panel_B,
  panel_C, panel_D,
  ncol = 2, nrow = 2,
  heights = c(1, 1.2)
)

# Save figure
ggsave("output/figures/publication_figure_4panel.png", final_figure, width = 12, height = 10, dpi = 300)
ggsave("output/figures/publication_figure_4panel.pdf", final_figure, width = 12, height = 10)

cat("\n==========================================================\n")
cat("                    OUTPUT FILES\n")
cat("==========================================================\n")
cat("  - publication_figure_4panel.png (300 DPI)\n")
cat("  - publication_figure_4panel.pdf\n")
cat(sprintf("\nOptimal Range: %.4f - %.4f μg/mL\n", optimal_lower, optimal_upper))
cat("==========================================================\n")
cat("                 Figure Complete!\n")
cat("==========================================================\n")
 