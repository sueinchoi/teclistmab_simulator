#===============================================================================
# Teclistamab PK Cut-off Analysis: Cavg_120hr
#
# Analyze different Cavg_120hr cut-offs for:
# - Response (≥VGPR) rate
# - CRS event rate
#
# Find optimal therapeutic window (maximize efficacy, minimize toxicity)
# Categorize patients into Low/Optimal/High exposure groups
#
# PREREQUISITE: Run pk_simulation.R first to generate pk_ae_merged_results.csv
#===============================================================================

library(tidyverse)
library(gridExtra)

cat("==========================================================\n")
cat("    Teclistamab Cavg_120hr Cut-off Analysis\n")
cat("==========================================================\n\n")

#-------------------------------------------------------------------------------
# 1. Load Data
#-------------------------------------------------------------------------------

cat("Loading data...\n")

# Load PK + AE merged data
if (!file.exists("pk_ae_merged_results.csv")) {
  stop("ERROR: pk_ae_merged_results.csv not found!\n",
       "Please run pk_simulation.R first.")
}

pk_ae_data <- read_csv("pk_ae_merged_results.csv", show_col_types = FALSE)
cat("  - Loaded PK+AE data for", nrow(pk_ae_data), "patients\n")

# Load response data
response_data <- read_csv("response_data.csv", show_col_types = FALSE)
cat("  - Loaded response data for", nrow(response_data), "patients\n")

# Load dosing data to filter for ≥6 doses
dosing_all <- read_csv("mrgsolve_dosing_full.csv", show_col_types = FALSE) %>%
  filter(TIME >= 0) %>%
  arrange(ID, TIME)

# Calculate number of doses per patient
dose_counts <- dosing_all %>%
  group_by(ID) %>%
  summarise(N_doses_actual = n(), .groups = "drop")

#-------------------------------------------------------------------------------
# 2. Merge All Data
#-------------------------------------------------------------------------------

cat("\n=== Merging data ===\n")

# Merge PK + AE + Response data
analysis_data <- pk_ae_data %>%
  left_join(response_data, by = "PID") %>%
  left_join(dose_counts, by = "ID") %>%
  mutate(
    # VGPR or better (sCR, CR, VGPR)
    VGPR_or_better = ifelse(RESP_CTX %in% c("sCR", "CR", "VGPR"), 1, 0),
    # At least 6 doses
    Has_6doses = ifelse(N_doses_actual >= 6, TRUE, FALSE)
  )

cat("  - Merged data:", nrow(analysis_data), "patients\n")

# Summary
cat("\n--- Data Summary ---\n")
cat("Total patients:", nrow(analysis_data), "\n")
cat("Patients with ≥6 doses:", sum(analysis_data$Has_6doses, na.rm = TRUE), "\n")
cat("CRS any grade events:", sum(analysis_data$CRS_any, na.rm = TRUE), "\n")
cat("VGPR+ responses (≥6 doses):", sum(analysis_data$VGPR_or_better[analysis_data$Has_6doses], na.rm = TRUE), "\n")

#-------------------------------------------------------------------------------
# 3. Define Cut-off Range for Cavg_120hr
#-------------------------------------------------------------------------------

cat("\n=== Cavg_120hr Distribution ===\n")

cavg_values <- analysis_data$Cavg_120hr
cat("Min:", round(min(cavg_values, na.rm = TRUE), 4), "\n")
cat("Max:", round(max(cavg_values, na.rm = TRUE), 4), "\n")
cat("Median:", round(median(cavg_values, na.rm = TRUE), 4), "\n")
cat("Mean:", round(mean(cavg_values, na.rm = TRUE), 4), "\n")

# Define cut-off range (from 10th to 90th percentile)
cutoff_min <- quantile(cavg_values, 0.05, na.rm = TRUE)
cutoff_max <- quantile(cavg_values, 0.95, na.rm = TRUE)
cutoffs <- seq(cutoff_min, cutoff_max, length.out = 50)

cat("\nCut-off range:", round(cutoff_min, 4), "to", round(cutoff_max, 4), "\n")

#-------------------------------------------------------------------------------
# 4. Calculate Event Rates at Each Cut-off
#-------------------------------------------------------------------------------

cat("\n=== Calculating Event Rates at Each Cut-off ===\n")

# For CRS: use all patients (CRS typically occurs early)
# For VGPR: use patients with ≥6 doses only (need adequate treatment)

event_rate_results <- map_dfr(cutoffs, function(cutoff) {
  # CRS analysis (all patients)
  crs_high <- analysis_data %>% filter(Cavg_120hr >= cutoff)
  crs_low <- analysis_data %>% filter(Cavg_120hr < cutoff)

  crs_rate_high <- mean(crs_high$CRS_any, na.rm = TRUE)
  crs_rate_low <- mean(crs_low$CRS_any, na.rm = TRUE)

  # VGPR analysis (≥6 doses only)
  vgpr_data <- analysis_data %>% filter(Has_6doses == TRUE)
  vgpr_high <- vgpr_data %>% filter(Cavg_120hr >= cutoff)
  vgpr_low <- vgpr_data %>% filter(Cavg_120hr < cutoff)

  vgpr_rate_high <- mean(vgpr_high$VGPR_or_better, na.rm = TRUE)
  vgpr_rate_low <- mean(vgpr_low$VGPR_or_better, na.rm = TRUE)

  tibble(
    Cutoff = cutoff,
    # Above cutoff
    N_above = nrow(crs_high),
    CRS_rate_above = crs_rate_high * 100,
    VGPR_rate_above = vgpr_rate_high * 100,
    # Below cutoff
    N_below = nrow(crs_low),
    CRS_rate_below = crs_rate_low * 100,
    VGPR_rate_below = vgpr_rate_low * 100
  )
})

#-------------------------------------------------------------------------------
# 5. Create Event Rate vs Cut-off Plot
#-------------------------------------------------------------------------------

cat("\n=== Creating Cut-off Analysis Plots ===\n")

# Plot 1: VGPR Rate vs Cut-off (High Cavg_120hr)
p1 <- ggplot(event_rate_results, aes(x = Cutoff)) +
  geom_line(aes(y = VGPR_rate_above), color = "#2ecc71", size = 1.5) +
  geom_point(aes(y = VGPR_rate_above), color = "#2ecc71", size = 2) +
  labs(
    title = "VGPR+ Rate by Cavg_120hr Cut-off",
    subtitle = "Patients with Cavg_120hr ≥ cut-off (among ≥6 doses)",
    x = "Cavg_120hr Cut-off (mg/L)",
    y = "VGPR+ Rate (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# Plot 2: CRS Rate vs Cut-off (High Cavg_120hr)
p2 <- ggplot(event_rate_results, aes(x = Cutoff)) +
  geom_line(aes(y = CRS_rate_above), color = "#e74c3c", size = 1.5) +
  geom_point(aes(y = CRS_rate_above), color = "#e74c3c", size = 2) +
  labs(
    title = "CRS Rate by Cavg_120hr Cut-off",
    subtitle = "Patients with Cavg_120hr ≥ cut-off (all patients)",
    x = "Cavg_120hr Cut-off (mg/L)",
    y = "CRS Rate (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# Combined plot: Both rates vs cut-off
plot_data_long <- event_rate_results %>%
  select(Cutoff, VGPR_rate_above, CRS_rate_above) %>%
  pivot_longer(
    cols = c(VGPR_rate_above, CRS_rate_above),
    names_to = "Outcome",
    values_to = "Rate"
  ) %>%
  mutate(
    Outcome = case_when(
      Outcome == "VGPR_rate_above" ~ "VGPR+ Response",
      Outcome == "CRS_rate_above" ~ "CRS Event"
    )
  )

p3 <- ggplot(plot_data_long, aes(x = Cutoff, y = Rate, color = Outcome)) +
  geom_line(size = 1.5) +
  geom_point(size = 2) +
  scale_color_manual(values = c("CRS Event" = "#e74c3c", "VGPR+ Response" = "#2ecc71")) +
  labs(
    title = "VGPR+ Response & CRS Rate by Cavg_120hr Cut-off",
    subtitle = "Finding optimal therapeutic window",
    x = "Cavg_120hr Cut-off (mg/L)",
    y = "Event Rate (%)",
    color = "Outcome"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave("cutoff_VGPR_rate.png", p1, width = 8, height = 6, dpi = 200)
ggsave("cutoff_CRS_rate.png", p2, width = 8, height = 6, dpi = 200)
ggsave("cutoff_combined_rates.png", p3, width = 10, height = 7, dpi = 200)

cat("Saved: cutoff_VGPR_rate.png\n")
cat("Saved: cutoff_CRS_rate.png\n")
cat("Saved: cutoff_combined_rates.png\n")

#-------------------------------------------------------------------------------
# 6. Find Optimal Range Using Benefit-Risk Analysis
#-------------------------------------------------------------------------------

cat("\n=== Finding Optimal Range (Benefit-Risk Analysis) ===\n")

library(pROC)

# Strategy: Find range where VGPR is maximized and CRS is acceptable
# 1. Upper bound: Where CRS rate stays below safety threshold
# 2. Lower bound: Where VGPR rate is acceptable
# 3. Net Clinical Benefit: VGPR_rate - CRS_rate

# Define thresholds
CRS_THRESHOLD <- 90   # Maximum acceptable CRS rate (%)
VGPR_MIN_THRESHOLD <- 50  # Minimum acceptable VGPR rate (%)

#--- Method 1: Constraint-based approach ---
cat("\n--- Method 1: Constraint-Based Approach ---\n")

# Calculate event rates at each Cavg_120hr value (within-bin analysis)
# Create finer bins for more precise analysis
n_bins <- 20
cavg_sorted <- sort(analysis_data$Cavg_120hr)
bin_breaks <- quantile(cavg_sorted, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)
bin_breaks <- unique(bin_breaks)  # Remove duplicates

# Calculate event rates for each bin
bin_analysis <- map_dfr(1:(length(bin_breaks) - 1), function(i) {
  lower <- bin_breaks[i]
  upper <- bin_breaks[i + 1]

  # All patients in this bin (for CRS)
  bin_all <- analysis_data %>%
    filter(Cavg_120hr >= lower & Cavg_120hr < upper)

  # Patients with ≥6 doses in this bin (for VGPR)
  bin_vgpr <- analysis_data %>%
    filter(Has_6doses == TRUE & Cavg_120hr >= lower & Cavg_120hr < upper)

  tibble(
    Bin = i,
    Lower = lower,
    Upper = upper,
    Midpoint = (lower + upper) / 2,
    N_all = nrow(bin_all),
    N_vgpr = nrow(bin_vgpr),
    CRS_rate = ifelse(nrow(bin_all) > 0, mean(bin_all$CRS_any, na.rm = TRUE) * 100, NA),
    VGPR_rate = ifelse(nrow(bin_vgpr) > 0, mean(bin_vgpr$VGPR_or_better, na.rm = TRUE) * 100, NA)
  )
})

cat("\n--- Event Rates by Cavg_120hr Bin ---\n")
print(bin_analysis %>% select(Bin, Lower, Upper, N_all, CRS_rate, N_vgpr, VGPR_rate), n = 25)

# Find upper bound: highest Cavg_120hr where CRS < threshold
upper_candidates <- bin_analysis %>%
  filter(!is.na(CRS_rate) & CRS_rate < CRS_THRESHOLD) %>%
  arrange(desc(Upper))

if (nrow(upper_candidates) > 0) {
  optimal_upper <- upper_candidates$Upper[1]
  cat(sprintf("\nUpper bound (CRS < %d%%): %.4f mg/L\n", CRS_THRESHOLD, optimal_upper))
} else {
  # If all bins have CRS >= threshold, use median
  optimal_upper <- median(cavg_values, na.rm = TRUE)
  cat(sprintf("\nNo bin with CRS < %d%%, using median: %.4f mg/L\n", CRS_THRESHOLD, optimal_upper))
}

# Find lower bound: lowest Cavg_120hr where VGPR is still acceptable
# Use cumulative approach: VGPR rate for patients >= cutoff
lower_candidates <- event_rate_results %>%
  filter(VGPR_rate_above >= VGPR_MIN_THRESHOLD & Cutoff <= optimal_upper) %>%
  arrange(Cutoff)

if (nrow(lower_candidates) > 0) {
  optimal_lower <- lower_candidates$Cutoff[1]
  cat(sprintf("Lower bound (VGPR >= %d%% above): %.4f mg/L\n", VGPR_MIN_THRESHOLD, optimal_lower))
} else {
  optimal_lower <- quantile(cavg_values, 0.25, na.rm = TRUE)
  cat(sprintf("No cutoff with VGPR >= %d%%, using Q1: %.4f mg/L\n", VGPR_MIN_THRESHOLD, optimal_lower))
}

#--- Method 2: Net Clinical Benefit ---
cat("\n--- Method 2: Net Clinical Benefit Score ---\n")

# Calculate Net Benefit = VGPR_rate - CRS_rate for patients within range [lower, upper]
# Test different ranges
range_analysis <- expand_grid(
  Lower = quantile(cavg_values, seq(0.1, 0.5, by = 0.1), na.rm = TRUE),
  Upper = quantile(cavg_values, seq(0.5, 0.9, by = 0.1), na.rm = TRUE)
) %>%
  filter(Lower < Upper) %>%
  rowwise() %>%
  mutate(
    # Patients within range
    N_range_all = sum(analysis_data$Cavg_120hr >= Lower & analysis_data$Cavg_120hr <= Upper, na.rm = TRUE),
    N_range_vgpr = sum(analysis_data$Cavg_120hr >= Lower & analysis_data$Cavg_120hr <= Upper &
                        analysis_data$Has_6doses == TRUE, na.rm = TRUE),

    # Event rates within range
    CRS_in_range = {
      pts <- analysis_data %>% filter(Cavg_120hr >= Lower & Cavg_120hr <= Upper)
      if (nrow(pts) > 0) mean(pts$CRS_any, na.rm = TRUE) * 100 else NA_real_
    },
    VGPR_in_range = {
      pts <- analysis_data %>% filter(Cavg_120hr >= Lower & Cavg_120hr <= Upper & Has_6doses == TRUE)
      if (nrow(pts) > 0) mean(pts$VGPR_or_better, na.rm = TRUE) * 100 else NA_real_
    },

    # Net Clinical Benefit (higher is better)
    Net_Benefit = VGPR_in_range - CRS_in_range,

    # Benefit-Risk Ratio (higher is better, avoid division by zero)
    Benefit_Risk_Ratio = VGPR_in_range / (CRS_in_range + 1)
  ) %>%
  ungroup() %>%
  filter(!is.na(Net_Benefit) & N_range_all >= 3)

# Find best range by Net Benefit (with constraint that CRS < threshold)
best_range <- range_analysis %>%
  filter(CRS_in_range < CRS_THRESHOLD) %>%
  arrange(desc(Net_Benefit)) %>%
  head(1)

if (nrow(best_range) > 0) {
  cat(sprintf("\nBest Range by Net Benefit (CRS < %d%%):\n", CRS_THRESHOLD))
  cat(sprintf("  Range: %.4f - %.4f mg/L\n", best_range$Lower, best_range$Upper))
  cat(sprintf("  N in range: %d (all), %d (≥6 doses)\n", best_range$N_range_all, best_range$N_range_vgpr))
  cat(sprintf("  VGPR rate: %.1f%%\n", best_range$VGPR_in_range))
  cat(sprintf("  CRS rate: %.1f%%\n", best_range$CRS_in_range))
  cat(sprintf("  Net Benefit: %.1f\n", best_range$Net_Benefit))

  # Use Method 2 result if it gives better range
  if (best_range$Net_Benefit > (optimal_lower - optimal_upper)) {
    optimal_lower_m2 <- best_range$Lower
    optimal_upper_m2 <- best_range$Upper
  }
}

#--- Final Optimal Range Selection ---
cat("\n--- Final Optimal Range ---\n")

# Ensure upper bound respects CRS constraint
# Check CRS rate at the proposed upper bound
crs_at_upper <- bin_analysis %>%
  filter(Upper <= optimal_upper) %>%
  summarise(max_crs = max(CRS_rate, na.rm = TRUE))

if (!is.na(crs_at_upper$max_crs) && crs_at_upper$max_crs >= 100) {
  # Find the last bin before CRS reaches 100%
  safe_bins <- bin_analysis %>%
    filter(CRS_rate < 100) %>%
    arrange(desc(Upper))

  if (nrow(safe_bins) > 0) {
    optimal_upper <- safe_bins$Upper[1]
    cat(sprintf("Adjusted upper bound (before CRS=100%%): %.4f mg/L\n", optimal_upper))
  }
}

# Validate bounds
if (optimal_lower >= optimal_upper) {
  cat("Warning: Lower bound >= Upper bound, adjusting...\n")
  optimal_lower <- quantile(cavg_values, 0.33, na.rm = TRUE)
  optimal_upper <- quantile(cavg_values, 0.67, na.rm = TRUE)
}

cat(sprintf("\n>>> OPTIMAL RANGE: %.4f - %.4f mg/L <<<\n", optimal_lower, optimal_upper))

#-------------------------------------------------------------------------------
# 7. Categorize Patients: Low / Optimal / High
#-------------------------------------------------------------------------------

cat("\n=== Categorizing Patients ===\n")

analysis_data <- analysis_data %>%
  mutate(
    Exposure_Category = case_when(
      Cavg_120hr < optimal_lower ~ "Low",
      Cavg_120hr >= optimal_lower & Cavg_120hr <= optimal_upper ~ "Optimal",
      Cavg_120hr > optimal_upper ~ "High"
    ),
    Exposure_Category = factor(Exposure_Category, levels = c("Low", "Optimal", "High"))
  )

cat("\nExposure Category Distribution:\n")
print(table(analysis_data$Exposure_Category))

#-------------------------------------------------------------------------------
# 8. Calculate Event Rates by Category
#-------------------------------------------------------------------------------

cat("\n=== Event Rates by Exposure Category ===\n")

# CRS rates (all patients)
crs_by_category <- analysis_data %>%
  group_by(Exposure_Category) %>%
  summarise(
    N = n(),
    CRS_events = sum(CRS_any, na.rm = TRUE),
    CRS_rate = mean(CRS_any, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cat("\n--- CRS Event Rate by Exposure Category ---\n")
print(crs_by_category)

# VGPR rates (≥6 doses only)
vgpr_by_category <- analysis_data %>%
  filter(Has_6doses == TRUE) %>%
  group_by(Exposure_Category) %>%
  summarise(
    N = n(),
    VGPR_events = sum(VGPR_or_better, na.rm = TRUE),
    VGPR_rate = mean(VGPR_or_better, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cat("\n--- VGPR+ Rate by Exposure Category (≥6 doses) ---\n")
print(vgpr_by_category)

# Combined summary
summary_by_category <- crs_by_category %>%
  left_join(
    vgpr_by_category %>% select(Exposure_Category, N_vgpr = N, VGPR_events, VGPR_rate),
    by = "Exposure_Category"
  )

cat("\n--- Combined Summary ---\n")
print(summary_by_category)

write_csv(summary_by_category, "pk_cutoff_category_summary.csv")
cat("\nSaved: pk_cutoff_category_summary.csv\n")

#-------------------------------------------------------------------------------
# 9. Create Bar Plot for Event Rates by Category
#-------------------------------------------------------------------------------

cat("\n=== Creating Category Bar Plots ===\n")

# Prepare data for plotting
bar_data <- bind_rows(
  crs_by_category %>%
    select(Exposure_Category, Rate = CRS_rate, N) %>%
    mutate(Outcome = "CRS Event"),
  vgpr_by_category %>%
    select(Exposure_Category, Rate = VGPR_rate, N) %>%
    mutate(Outcome = "VGPR+ Response")
)

# Bar plot
p_bar <- ggplot(bar_data, aes(x = Exposure_Category, y = Rate, fill = Outcome)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", Rate)),
            position = position_dodge(width = 0.8), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("CRS Event" = "#e74c3c", "VGPR+ Response" = "#2ecc71")) +
  labs(
    title = "Event Rates by Cavg_120hr Exposure Category",
    subtitle = sprintf("Optimal Range: %.4f - %.4f mg/L", optimal_lower, optimal_upper),
    x = "Exposure Category",
    y = "Event Rate (%)",
    fill = "Outcome"
  ) +
  ylim(0, max(bar_data$Rate, na.rm = TRUE) * 1.15) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave("cutoff_category_barplot.png", p_bar, width = 10, height = 7, dpi = 200)
cat("Saved: cutoff_category_barplot.png\n")

# Separate plots for each outcome
p_crs_bar <- ggplot(crs_by_category, aes(x = Exposure_Category, y = CRS_rate, fill = Exposure_Category)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", CRS_rate, N)), vjust = -0.3, size = 4) +
  scale_fill_manual(values = c("Low" = "#3498db", "Optimal" = "#2ecc71", "High" = "#e74c3c")) +
  labs(
    title = "CRS Event Rate by Exposure Category",
    subtitle = sprintf("Low: <%.4f | Optimal: %.4f-%.4f | High: >%.4f",
                       optimal_lower, optimal_lower, optimal_upper, optimal_upper),
    x = "Cavg_120hr Category",
    y = "CRS Rate (%)"
  ) +
  ylim(0, max(crs_by_category$CRS_rate, na.rm = TRUE) * 1.25) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "none")

p_vgpr_bar <- ggplot(vgpr_by_category, aes(x = Exposure_Category, y = VGPR_rate, fill = Exposure_Category)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", VGPR_rate, N)), vjust = -0.3, size = 4) +
  scale_fill_manual(values = c("Low" = "#3498db", "Optimal" = "#2ecc71", "High" = "#e74c3c")) +
  labs(
    title = "VGPR+ Response Rate by Exposure Category",
    subtitle = sprintf("Low: <%.4f | Optimal: %.4f-%.4f | High: >%.4f (≥6 doses)",
                       optimal_lower, optimal_lower, optimal_upper, optimal_upper),
    x = "Cavg_120hr Category",
    y = "VGPR+ Rate (%)"
  ) +
  ylim(0, max(vgpr_by_category$VGPR_rate, na.rm = TRUE) * 1.25) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "none")

ggsave("cutoff_CRS_barplot.png", p_crs_bar, width = 8, height = 6, dpi = 200)
ggsave("cutoff_VGPR_barplot.png", p_vgpr_bar, width = 8, height = 6, dpi = 200)

cat("Saved: cutoff_CRS_barplot.png\n")
cat("Saved: cutoff_VGPR_barplot.png\n")

#-------------------------------------------------------------------------------
# 10. Combined Cut-off Plot with Optimal Range Shading
#-------------------------------------------------------------------------------

p_combined_shaded <- ggplot(plot_data_long, aes(x = Cutoff, y = Rate, color = Outcome)) +
  annotate("rect", xmin = optimal_lower, xmax = optimal_upper,
           ymin = -Inf, ymax = Inf, fill = "#2ecc71", alpha = 0.2) +
  geom_line(size = 1.5) +
  geom_point(size = 2) +
  geom_vline(xintercept = optimal_lower, linetype = "dashed", color = "darkgreen", size = 0.8) +
  geom_vline(xintercept = optimal_upper, linetype = "dashed", color = "darkgreen", size = 0.8) +
  scale_color_manual(values = c("CRS Event" = "#e74c3c", "VGPR+ Response" = "#2ecc71")) +
  annotate("text", x = (optimal_lower + optimal_upper) / 2, y = 5,
           label = "Optimal\nRange", color = "darkgreen", fontface = "bold", size = 4) +
  labs(
    title = "Therapeutic Window Analysis: Cavg_120hr",
    subtitle = sprintf("Optimal Range: %.4f - %.4f mg/L (shaded green)", optimal_lower, optimal_upper),
    x = "Cavg_120hr (mg/L)",
    y = "Event Rate (%)",
    color = "Outcome"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave("cutoff_therapeutic_window.png", p_combined_shaded, width = 12, height = 8, dpi = 200)
cat("Saved: cutoff_therapeutic_window.png\n")

#-------------------------------------------------------------------------------
# 11. Print Final Summary
#-------------------------------------------------------------------------------

cat("\n")
cat("==========================================================\n")
cat("              CAVG_120HR CUT-OFF ANALYSIS SUMMARY\n")
cat("==========================================================\n\n")

cat(sprintf("Optimal Therapeutic Range: %.4f - %.4f mg/L\n\n", optimal_lower, optimal_upper))

cat("--- Event Rates by Category ---\n\n")
cat("Exposure Category | CRS Rate (%) | VGPR+ Rate (%)\n")
cat("------------------|--------------|---------------\n")
for (i in 1:nrow(summary_by_category)) {
  row <- summary_by_category[i, ]
  cat(sprintf("%-17s | %12.1f | %14.1f\n",
              as.character(row$Exposure_Category),
              row$CRS_rate,
              ifelse(is.na(row$VGPR_rate), 0, row$VGPR_rate)))
}

cat("\n==========================================================\n")
cat("                    OUTPUT FILES\n")
cat("==========================================================\n")
cat("  - pk_cutoff_category_summary.csv\n")
cat("  - cutoff_VGPR_rate.png\n")
cat("  - cutoff_CRS_rate.png\n")
cat("  - cutoff_combined_rates.png\n")
cat("  - cutoff_category_barplot.png\n")
cat("  - cutoff_CRS_barplot.png\n")
cat("  - cutoff_VGPR_barplot.png\n")
cat("  - cutoff_therapeutic_window.png\n")
cat("==========================================================\n")
cat("                 Analysis Complete!\n")
cat("==========================================================\n")
