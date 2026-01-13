#===============================================================================
# Teclistamab PK Monte Carlo Simulation
#
# Individual PK simulations with Monte Carlo (1000 iterations)
# Outputs: pk_ae_merged_results.csv (PK metrics merged with AE data)
#
# Run this script ONCE to generate PK results, then use:
# - ae_pk_analysis.R for AE association analysis
# - pk_response_analysis.R for response association analysis
#===============================================================================

library(mrgsolve)
library(tidyverse)

set.seed(12345)

cat("==========================================================\n")
cat("    Teclistamab PK Monte Carlo Simulation\n")
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
mod <- mcode("teclistamab_sim", model_code, compile = TRUE)

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
params <- read_csv("output/tables/mrgsolve_params_full.csv", show_col_types = FALSE)
cat("  - Loaded", nrow(params), "patients from params file\n")

# Dosing data
dosing_all <- read_csv("output/tables/mrgsolve_dosing_full.csv", show_col_types = FALSE) %>%
  filter(TIME >= 0) %>%
  arrange(ID, TIME)
cat("  - Loaded", nrow(dosing_all), "dosing records\n")

# AE data (CRS, Neuro from existing file)
ae_crs_neuro <- read_csv("output/tables/ae_pk_metrics_final.csv", show_col_types = FALSE) %>%
  select(ID, PID, CRS_any, CRS_gr2, Neuro_any, Neuro_gr2,
         Neutropenia_gr3, Thrombocytopenia_gr3, Lymphopenia_gr3, Anemia_gr2, Other_any)
cat("  - Loaded CRS/Neuro AE data for", nrow(ae_crs_neuro), "patients\n")

# Detailed AE data (from Excel extraction)
ae_detailed <- read_csv("data/ae_detailed_data.csv", show_col_types = FALSE)
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
# 6. Run Monte Carlo Simulation for All Patients
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
# 7. Merge PK Results with AE Data and Save
#-------------------------------------------------------------------------------

cat("\nMerging PK results with AE data...\n")

analysis_data <- pk_results %>%
  left_join(ae_combined, by = c("ID", "PID"))

write_csv(analysis_data, "output/tables/pk_ae_merged_results.csv")
cat("Saved: output/tables/pk_ae_merged_results.csv\n")

# Also save PK-only results
write_csv(pk_results, "output/tables/pk_simulation_results.csv")
cat("Saved: output/tables/pk_simulation_results.csv\n")

#-------------------------------------------------------------------------------
# Summary
#-------------------------------------------------------------------------------

cat("\n==========================================================\n")
cat("                    SIMULATION COMPLETE\n")
cat("==========================================================\n\n")

cat("PK Metrics Summary:\n")
print(summary(pk_results %>% select(starts_with("Cmax"), starts_with("Cavg"))))

cat("\n==========================================================\n")
cat("                    OUTPUT FILES\n")
cat("==========================================================\n")
cat("  - pk_simulation_results.csv (PK metrics only)\n")
cat("  - pk_ae_merged_results.csv (PK + AE data merged)\n")
cat("\nNext steps:\n")
cat("  - Run ae_pk_analysis.R for AE association analysis\n")
cat("  - Run pk_response_analysis.R for response analysis\n")
cat("==========================================================\n")
