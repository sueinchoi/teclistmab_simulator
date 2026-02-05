
#===============================================================================
# Teclistamab Population PK Simulation with mrgsolve
# 
# Model: 2-Compartment with Time-Dependent Clearance
#        CL(t) = CL1 + CL2 * exp(-KDES * t)
#
# IMPORTANT: Time unit is HOURS, but rate constants are /DAY
#            → Must convert rate constants to /HOUR in ODE
#
# Reference: Miao et al. (2023) - Teclistamab Population PK Model
#===============================================================================

library(mrgsolve)
library(tidyverse)

#-------------------------------------------------------------------------------
# 1. Define the PK model
#-------------------------------------------------------------------------------

code <- '
$PARAM @annotated
// Individual parameters (all rate constants in /day units)
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
// CRITICAL: TIME is in HOURS, parameters are in /DAY
// Convert time to days for time-dependent CL calculation
double TDAY = SOLVERTIME / 24.0;

// Time-dependent clearance (L/day)
double CL = CL1 + CL2 * exp(-KDES * TDAY);

// CRITICAL: Convert rate constants from /day to /hour
// because SOLVERTIME is in hours
double KA_HR  = KA / 24.0;
double K10_HR = (CL / V1) / 24.0;
double K12_HR = (Q / V1) / 24.0;
double K21_HR = (Q / V2) / 24.0;

// ODEs (using /hour rate constants)
dxdt_DEPOT  = -KA_HR * DEPOT;
dxdt_CENT   = KA_HR * DEPOT * F1 - K10_HR * CENT - K12_HR * CENT + K21_HR * PERIPH;
dxdt_PERIPH = K12_HR * CENT - K21_HR * PERIPH;

$TABLE
double DV = CP;
double TDAY_OUT = TIME / 24.0;
double CL_T = CL1 + CL2 * exp(-KDES * TDAY_OUT);

$CAPTURE @annotated
DV      : Plasma concentration (mg/L)
TDAY_OUT: Time in days
CL_T    : Time-dependent clearance at this time (L/day)
'

# Compile the model
mod <- mcode("teclistamab_2cmt", code)

#-------------------------------------------------------------------------------
# 2. Load patient data
#-------------------------------------------------------------------------------

# Load individual PK parameters
params <- read_csv("mrgsolve_params_full.csv")
cat("\n=== Patient Parameters ===\n")
print(params %>% select(ID, WT, ISS, IGG, CL1, CL2, V1, V2, Q, KA, N_doses))

# Load dosing records (TIME in hours)
dosing <- read_csv("mrgsolve_dosing_full.csv")
cat("\n=== Dosing Records (first 20) ===\n")
print(head(dosing, 20))

#-------------------------------------------------------------------------------
# 3. Simulation function
#-------------------------------------------------------------------------------

simulate_patient <- function(pt_id, params_df, dosing_df) {
  
  pt_params <- params_df %>% filter(ID == pt_id)
  if (nrow(pt_params) == 0) return(NULL)
  
  pt_dosing <- dosing_df %>% 
    filter(ID == pt_id) %>%
    rename(time = TIME, amt = AMT, cmt = CMT, evid = EVID) %>%
    mutate(ID = pt_id)
  
  if (nrow(pt_dosing) == 0) return(NULL)
  
  # Simulation time (hourly, last dose + 21 days)
  max_time <- max(pt_dosing$time) + 21 * 24
  sim_times <- seq(0, max_time, by = 1)
  
  obs_data <- tibble(
    ID = pt_id,
    time = sim_times,
    amt = 0,
    cmt = 0,
    evid = 0
  )
  
  sim_data <- bind_rows(pt_dosing, obs_data) %>%
    arrange(time, desc(evid))
  
  # Update model parameters
  mod_i <- mod %>%
    param(
      CL1  = pt_params$CL1,
      CL2  = pt_params$CL2,
      KDES = pt_params$KDES,
      V1   = pt_params$V1,
      V2   = pt_params$V2,
      Q    = pt_params$Q,
      KA   = pt_params$KA,
      F1   = pt_params$F1
    )
  
  out <- mod_i %>%
    data_set(sim_data) %>%
    mrgsim(carry_out = "amt,evid") %>%
    as_tibble() %>%
    filter(evid == 0) %>%
    mutate(
      TIME_HOUR = time,
      TIME_DAY = time / 24,
      PID = pt_params$PID,
      WT = pt_params$WT,
      ISS = pt_params$ISS,
      IGG = pt_params$IGG,
      N_doses = pt_params$N_doses
    )
  
  return(out)
}

#-------------------------------------------------------------------------------
# 4. Run simulation
#-------------------------------------------------------------------------------

cat("\nRunning simulation...\n")

all_results <- map_dfr(
  params$ID, 
  ~simulate_patient(.x, params, dosing),
  .progress = TRUE
)

write_csv(all_results, "mrgsolve_simulation_correct.csv")
cat("\nSimulation saved!\n")

#-------------------------------------------------------------------------------
# 5. Create plots (5x6 layout)
#-------------------------------------------------------------------------------

plot_data <- all_results %>% filter(TIME_DAY <= 100)

# ISS colors
iss_colors <- c("1" = "#2ecc71", "2" = "#3498db", "3" = "#e74c3c")

p <- ggplot(plot_data, aes(x = TIME_DAY, y = DV)) +
  geom_line(aes(color = factor(ISS)), linewidth = 0.8) +
  facet_wrap(~paste0("ID ", ID, " (N=", N_doses, ")"), 
             scales = "free", ncol = 6) +
  scale_color_manual(values = iss_colors, name = "ISS Stage") +
  labs(
    x = "Time (days)",
    y = "Concentration (mg/L)",
    title = "Teclistamab Individual PK Profiles (n=30)",
    subtitle = "2-Compartment with Time-Dependent Clearance"
  ) +
  theme_bw() +
  theme(strip.text = element_text(size = 8), legend.position = "bottom")

ggsave("individual_pk_mrgsolve_5x6.png", p, width = 20, height = 18, dpi = 200)
ggsave("individual_pk_mrgsolve_5x6.pdf", p, width = 20, height = 18)

#-------------------------------------------------------------------------------
# 6. PK metrics
#-------------------------------------------------------------------------------

pk_metrics <- all_results %>%
  group_by(ID, PID, WT, ISS, IGG, N_doses) %>%
  summarise(
    Cmax = max(DV),
    Cmin = min(DV[DV > 0.01]),
    Cavg = mean(DV),
    .groups = "drop"
  )

write_csv(pk_metrics, "mrgsolve_pk_metrics_correct.csv")
print(pk_metrics)

cat("\n=== Complete! ===\n")
