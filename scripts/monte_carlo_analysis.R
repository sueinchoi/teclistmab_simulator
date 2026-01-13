#===============================================================================
# Teclistamab Monte Carlo PK Simulation Analysis
#
# 1000 Monte Carlo simulations per subject with IIV
# Calculate: Cavg, Cmax at 72hr, 120hr, 1st/3rd dosing intervals
#
# Model: 2-Compartment with Time-Dependent Clearance
# Uses individual dosing data from mrgsolve_dosing_full.csv
#===============================================================================

library(mrgsolve)
library(tidyverse)

# Set seed for reproducibility
set.seed(12345)

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

# Compile model
cat("Compiling mrgsolve model...\n")
mod <- mcode("teclistamab_mc", model_code, compile = TRUE)

#-------------------------------------------------------------------------------
# 2. IIV Parameters (CV to omega conversion)
#-------------------------------------------------------------------------------

cv_to_omega <- function(cv) sqrt(log(1 + cv^2))

IIV_OMEGA <- list(
  CL1 = cv_to_omega(0.536),   # 53.6%
  CL2 = cv_to_omega(1.07),    # 107%
  V1  = cv_to_omega(0.488),   # 48.8%
  KA  = cv_to_omega(0.452)    # 45.2%
)

#-------------------------------------------------------------------------------
# 3. Load Patient Parameters and Dosing Data
#-------------------------------------------------------------------------------

params <- read_csv("mrgsolve_params_full.csv", show_col_types = FALSE)
cat("\n=== Loaded", nrow(params), "patients ===\n")

# Load individual dosing data
dosing_all <- read_csv("mrgsolve_dosing_full.csv", show_col_types = FALSE) %>%
  filter(TIME >= 0) %>%  # Remove invalid negative times
  arrange(ID, TIME)

cat("=== Loaded dosing records for", n_distinct(dosing_all$ID), "patients ===\n")

# Show dosing summary
dosing_summary <- dosing_all %>%
  group_by(ID) %>%
  summarise(
    N_doses = n(),
    First_dose_time = min(TIME),
    Last_dose_time = max(TIME),
    Total_dose_mg = sum(AMT),
    .groups = "drop"
  )
print(dosing_summary)

#-------------------------------------------------------------------------------
# 4. Monte Carlo Simulation Function (Using Individual Dosing Data)
#-------------------------------------------------------------------------------

run_monte_carlo <- function(pt_id, params_df, dosing_df, n_sim = 1000, mod) {

  pt_params <- params_df %>% filter(ID == pt_id)
  if (nrow(pt_params) == 0) return(NULL)

  # Get patient-specific dosing data
  pt_dosing <- dosing_df %>%
    filter(ID == pt_id) %>%
    arrange(TIME) %>%
    select(time = TIME, amt = AMT, cmt = CMT, evid = EVID)

  if (nrow(pt_dosing) == 0) return(NULL)

  # Get dosing times for interval calculations
  dose_times <- pt_dosing$time

  # Simulation time: extend beyond last dose
  max_time <- max(pt_dosing$time) + 21 * 24  # Last dose + 21 days
  sim_times <- seq(0, max_time, by = 1)  # Hourly

  # Store results for all simulations
  all_sim_results <- vector("list", n_sim)

  for (sim_i in 1:n_sim) {
    # Generate individual parameters with IIV
    eta_CL1 <- rnorm(1, 0, IIV_OMEGA$CL1)
    eta_CL2 <- rnorm(1, 0, IIV_OMEGA$CL2)
    eta_V1  <- rnorm(1, 0, IIV_OMEGA$V1)
    eta_KA  <- rnorm(1, 0, IIV_OMEGA$KA)

    ind_CL1 <- pt_params$CL1 * exp(eta_CL1)
    ind_CL2 <- pt_params$CL2 * exp(eta_CL2)
    ind_V1  <- pt_params$V1 * exp(eta_V1)
    ind_KA  <- pt_params$KA * exp(eta_KA)

    # Prepare simulation data
    obs_data <- tibble(
      ID = 1,
      time = sim_times,
      amt = 0,
      cmt = 0,
      evid = 0
    )

    dosing_data <- pt_dosing %>%
      mutate(ID = 1)

    sim_data <- bind_rows(dosing_data, obs_data) %>%
      arrange(time, desc(evid))

    # Update model parameters
    mod_i <- mod %>%
      param(
        CL1  = ind_CL1,
        CL2  = ind_CL2,
        KDES = pt_params$KDES,
        V1   = ind_V1,
        V2   = pt_params$V2,
        Q    = pt_params$Q,
        KA   = ind_KA,
        F1   = pt_params$F1
      )

    # Run simulation
    out <- mod_i %>%
      data_set(sim_data) %>%
      mrgsim(carry_out = "amt,evid") %>%
      as_tibble() %>%
      filter(evid == 0) %>%
      select(time, DV)

    all_sim_results[[sim_i]] <- out %>%
      mutate(SIM = sim_i)
  }

  # Combine all simulations
  combined <- bind_rows(all_sim_results) %>%
    mutate(
      PT_ID = pt_id,
      PID = pt_params$PID,
      WT = pt_params$WT
    )

  list(
    simulations = combined,
    dose_times = dose_times,
    pt_info = pt_params
  )
}

#-------------------------------------------------------------------------------
# 5. PK Metrics Calculation Functions
#-------------------------------------------------------------------------------

# Calculate Cmax and Cavg up to a specific hour
calc_metrics_up_to_hour <- function(sim_data, end_hour) {
  sim_data %>%
    filter(time >= 0 & time <= end_hour) %>%
    group_by(SIM) %>%
    summarise(
      Cmax = max(DV, na.rm = TRUE),
      # Cavg = AUC / time (trapezoidal rule)
      AUC = {
        auc <- 0
        dv <- DV
        t <- time
        if (length(dv) > 1) {
          for (i in 2:length(dv)) {
            dt <- (t[i] - t[i-1]) / 24  # Convert to days
            avg_c <- (dv[i] + dv[i-1]) / 2
            auc <- auc + dt * avg_c
          }
        }
        auc
      },
      .groups = "drop"
    ) %>%
    mutate(Cavg = AUC / (end_hour / 24))
}

# Calculate Cmax and Cavg for a specific dosing interval
calc_metrics_dosing_interval <- function(sim_data, dose_times, interval_num) {
  if (length(dose_times) < interval_num) {
    return(tibble(SIM = unique(sim_data$SIM), Cmax = NA_real_, Cavg = NA_real_, AUC = NA_real_, interval_hours = NA_real_))
  }

  start_time <- dose_times[interval_num]

  # End time is next dose or +168hr (7 days) if last dose
  if (length(dose_times) > interval_num) {
    end_time <- dose_times[interval_num + 1]
  } else {
    end_time <- start_time + 168  # 7 days
  }

  sim_data %>%
    filter(time >= start_time & time <= end_time) %>%
    group_by(SIM) %>%
    summarise(
      Cmax = max(DV, na.rm = TRUE),
      AUC = {
        auc <- 0
        dv <- DV
        t <- time
        if (length(dv) > 1) {
          for (i in 2:length(dv)) {
            dt <- (t[i] - t[i-1]) / 24
            avg_c <- (dv[i] + dv[i-1]) / 2
            auc <- auc + dt * avg_c
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

cat("\n=== Running Monte Carlo Analysis (", n_simulations, " simulations per patient) ===\n")
cat("This may take a while...\n\n")

# Store all results
all_results <- list()

for (i in 1:nrow(params)) {
  pt_id <- params$ID[i]
  cat(sprintf("Processing Patient %d/%d (ID: %d, PID: %d)...\n",
              i, nrow(params), pt_id, params$PID[i]))

  # Run Monte Carlo simulation
  mc_result <- run_monte_carlo(pt_id, params, dosing_all, n_sim = n_simulations, mod = mod)

  if (is.null(mc_result)) {
    cat("  Skipped (no data)\n")
    next
  }

  sim_data <- mc_result$simulations
  dose_times <- mc_result$dose_times
  pt_info <- mc_result$pt_info

  # Calculate metrics at 72hr and 120hr
  metrics_72hr <- calc_metrics_up_to_hour(sim_data, 72)
  metrics_120hr <- calc_metrics_up_to_hour(sim_data, 120)

  # Calculate metrics for 1st and 3rd dosing intervals
  metrics_dose1 <- calc_metrics_dosing_interval(sim_data, dose_times, 1)
  metrics_dose3 <- calc_metrics_dosing_interval(sim_data, dose_times, 3)

  # Compile summary for this patient
  pt_summary <- tibble(
    ID = pt_id,
    PID = pt_info$PID,
    WT = pt_info$WT,
    ISS = pt_info$ISS,
    IGG = pt_info$IGG,
    N_doses = length(dose_times),

    # Dosing interval info
    Dose1_interval_hr = ifelse(length(dose_times) >= 2, dose_times[2] - dose_times[1], NA_real_),
    Dose3_interval_hr = ifelse(length(dose_times) >= 4, dose_times[4] - dose_times[3], NA_real_),

    # 72hr metrics
    Cmax_72hr_mean = mean(metrics_72hr$Cmax, na.rm = TRUE),
    Cmax_72hr_median = median(metrics_72hr$Cmax, na.rm = TRUE),
    Cmax_72hr_q5 = quantile(metrics_72hr$Cmax, 0.05, na.rm = TRUE),
    Cmax_72hr_q95 = quantile(metrics_72hr$Cmax, 0.95, na.rm = TRUE),
    Cavg_72hr_mean = mean(metrics_72hr$Cavg, na.rm = TRUE),
    Cavg_72hr_median = median(metrics_72hr$Cavg, na.rm = TRUE),
    Cavg_72hr_q5 = quantile(metrics_72hr$Cavg, 0.05, na.rm = TRUE),
    Cavg_72hr_q95 = quantile(metrics_72hr$Cavg, 0.95, na.rm = TRUE),

    # 120hr metrics
    Cmax_120hr_mean = mean(metrics_120hr$Cmax, na.rm = TRUE),
    Cmax_120hr_median = median(metrics_120hr$Cmax, na.rm = TRUE),
    Cmax_120hr_q5 = quantile(metrics_120hr$Cmax, 0.05, na.rm = TRUE),
    Cmax_120hr_q95 = quantile(metrics_120hr$Cmax, 0.95, na.rm = TRUE),
    Cavg_120hr_mean = mean(metrics_120hr$Cavg, na.rm = TRUE),
    Cavg_120hr_median = median(metrics_120hr$Cavg, na.rm = TRUE),
    Cavg_120hr_q5 = quantile(metrics_120hr$Cavg, 0.05, na.rm = TRUE),
    Cavg_120hr_q95 = quantile(metrics_120hr$Cavg, 0.95, na.rm = TRUE),

    # 1st dosing interval metrics
    Cmax_dose1_mean = mean(metrics_dose1$Cmax, na.rm = TRUE),
    Cmax_dose1_median = median(metrics_dose1$Cmax, na.rm = TRUE),
    Cmax_dose1_q5 = quantile(metrics_dose1$Cmax, 0.05, na.rm = TRUE),
    Cmax_dose1_q95 = quantile(metrics_dose1$Cmax, 0.95, na.rm = TRUE),
    Cavg_dose1_mean = mean(metrics_dose1$Cavg, na.rm = TRUE),
    Cavg_dose1_median = median(metrics_dose1$Cavg, na.rm = TRUE),
    Cavg_dose1_q5 = quantile(metrics_dose1$Cavg, 0.05, na.rm = TRUE),
    Cavg_dose1_q95 = quantile(metrics_dose1$Cavg, 0.95, na.rm = TRUE),

    # 3rd dosing interval metrics
    Cmax_dose3_mean = mean(metrics_dose3$Cmax, na.rm = TRUE),
    Cmax_dose3_median = median(metrics_dose3$Cmax, na.rm = TRUE),
    Cmax_dose3_q5 = quantile(metrics_dose3$Cmax, 0.05, na.rm = TRUE),
    Cmax_dose3_q95 = quantile(metrics_dose3$Cmax, 0.95, na.rm = TRUE),
    Cavg_dose3_mean = mean(metrics_dose3$Cavg, na.rm = TRUE),
    Cavg_dose3_median = median(metrics_dose3$Cavg, na.rm = TRUE),
    Cavg_dose3_q5 = quantile(metrics_dose3$Cavg, 0.05, na.rm = TRUE),
    Cavg_dose3_q95 = quantile(metrics_dose3$Cavg, 0.95, na.rm = TRUE)
  )

  all_results[[i]] <- pt_summary

  cat(sprintf("  Doses: %d, Cmax(72hr): %.4f [%.4f-%.4f], Cavg(72hr): %.4f [%.4f-%.4f]\n",
              length(dose_times),
              pt_summary$Cmax_72hr_median, pt_summary$Cmax_72hr_q5, pt_summary$Cmax_72hr_q95,
              pt_summary$Cavg_72hr_median, pt_summary$Cavg_72hr_q5, pt_summary$Cavg_72hr_q95))
}

#-------------------------------------------------------------------------------
# 7. Combine and Save Results
#-------------------------------------------------------------------------------

final_results <- bind_rows(all_results)

cat("\n=== Monte Carlo Analysis Complete ===\n")
cat("Patients analyzed:", nrow(final_results), "\n")
cat("Simulations per patient:", n_simulations, "\n\n")

# Save results
write_csv(final_results, "monte_carlo_pk_results.csv")
cat("Results saved to: monte_carlo_pk_results.csv\n")

# Print summary
cat("\n=== Summary Statistics (Median [90% PI]) ===\n\n")

summary_table <- final_results %>%
  summarise(
    `Cmax 72hr` = sprintf("%.4f [%.4f-%.4f]",
                          median(Cmax_72hr_median),
                          quantile(Cmax_72hr_median, 0.05),
                          quantile(Cmax_72hr_median, 0.95)),
    `Cavg 72hr` = sprintf("%.4f [%.4f-%.4f]",
                          median(Cavg_72hr_median),
                          quantile(Cavg_72hr_median, 0.05),
                          quantile(Cavg_72hr_median, 0.95)),
    `Cmax 120hr` = sprintf("%.4f [%.4f-%.4f]",
                           median(Cmax_120hr_median),
                           quantile(Cmax_120hr_median, 0.05),
                           quantile(Cmax_120hr_median, 0.95)),
    `Cavg 120hr` = sprintf("%.4f [%.4f-%.4f]",
                           median(Cavg_120hr_median),
                           quantile(Cavg_120hr_median, 0.05),
                           quantile(Cavg_120hr_median, 0.95)),
    `Cmax Dose1` = sprintf("%.4f [%.4f-%.4f]",
                           median(Cmax_dose1_median),
                           quantile(Cmax_dose1_median, 0.05),
                           quantile(Cmax_dose1_median, 0.95)),
    `Cavg Dose1` = sprintf("%.4f [%.4f-%.4f]",
                           median(Cavg_dose1_median),
                           quantile(Cavg_dose1_median, 0.05),
                           quantile(Cavg_dose1_median, 0.95)),
    `Cmax Dose3` = sprintf("%.4f [%.4f-%.4f]",
                           median(Cmax_dose3_median, na.rm = TRUE),
                           quantile(Cmax_dose3_median, 0.05, na.rm = TRUE),
                           quantile(Cmax_dose3_median, 0.95, na.rm = TRUE)),
    `Cavg Dose3` = sprintf("%.4f [%.4f-%.4f]",
                           median(Cavg_dose3_median, na.rm = TRUE),
                           quantile(Cavg_dose3_median, 0.05, na.rm = TRUE),
                           quantile(Cavg_dose3_median, 0.95, na.rm = TRUE))
  )

print(t(summary_table))

#-------------------------------------------------------------------------------
# 8. Detailed Results Table
#-------------------------------------------------------------------------------

cat("\n=== Per-Patient Results ===\n")
print(final_results %>%
        select(ID, PID, WT, N_doses,
               Cmax_72hr_median, Cavg_72hr_median,
               Cmax_120hr_median, Cavg_120hr_median,
               Cmax_dose1_median, Cavg_dose1_median,
               Cmax_dose3_median, Cavg_dose3_median))

cat("\n=== Done! ===\n")
