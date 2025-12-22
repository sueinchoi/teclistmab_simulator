#===============================================================================
# Teclistamab PK Simulation Shiny App
#
# Individual PK parameter calculation based on patient characteristics
# 2-Compartment model with time-dependent clearance
# Monte Carlo simulation with inter-individual variability (IIV)
# CRS Risk Assessment based on Cavg
#
# Reference: Miao et al. (2023) - Teclistamab Population PK Model
#===============================================================================

library(shiny)
library(mrgsolve)
library(tidyverse)
library(DT)

#-------------------------------------------------------------------------------
# mrgsolve Model Definition
#-------------------------------------------------------------------------------

model_code <- '
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

# Compile model at startup
mod <- mcode("teclistamab_app", model_code, compile = TRUE)

#-------------------------------------------------------------------------------
# IIV Parameters (CV to omega conversion)
# omega = sqrt(log(1 + CV^2))
#-------------------------------------------------------------------------------

cv_to_omega <- function(cv) {
  sqrt(log(1 + cv^2))
}

# Inter-individual variability (CV)
IIV_CV <- list(
  CL1 = 0.536,   # 53.6%
  CL2 = 1.07,    # 107%
  V1  = 0.488,   # 48.8%
  KA  = 0.452    # 45.2%
)

# Convert to omega
IIV_OMEGA <- list(
  CL1 = cv_to_omega(IIV_CV$CL1),
  CL2 = cv_to_omega(IIV_CV$CL2),
  V1  = cv_to_omega(IIV_CV$V1),
  KA  = cv_to_omega(IIV_CV$KA)
)

#-------------------------------------------------------------------------------
# PK Parameter Calculation Functions (Typical Values)
#-------------------------------------------------------------------------------

calculate_pk_params_typical <- function(bw, iss, igg_type) {
  # ISS indicator variables
  iss_2 <- ifelse(iss == "II", 1, 0)
  iss_3 <- ifelse(iss == "III", 1, 0)

  # IgG indicator (Non-IgG = 1)
  non_igg <- ifelse(igg_type == "Non-IgG", 1, 0)

  # Calculate CL1 (L/day)
  CL1 <- 0.449 * (bw / 74)^0.704 * (1.31^iss_2) * (1.67^iss_3) * (0.689^non_igg)

  # Calculate CL2 (L/day)
  CL2 <- 0.547 * (0.295^non_igg)

  # Calculate V1 (L)
  V1 <- 4.13 * (bw / 74)^0.358

  # Calculate V2 (L)
  V2 <- 1.34 * (bw / 74)^1.40

  # Fixed parameters
  Q <- 0.039
  KA <- 0.133
  F1 <- 0.718
  KDES <- 0.0292

  list(
    CL1 = CL1, CL2 = CL2, V1 = V1, V2 = V2,
    Q = Q, KA = KA, F1 = F1, KDES = KDES
  )
}

#-------------------------------------------------------------------------------
# Generate Individual Parameters with IIV (Monte Carlo)
#-------------------------------------------------------------------------------

generate_individual_params <- function(typical_params, n_subjects, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Generate eta values (random effects)
  eta_CL1 <- rnorm(n_subjects, 0, IIV_OMEGA$CL1)
  eta_CL2 <- rnorm(n_subjects, 0, IIV_OMEGA$CL2)
  eta_V1  <- rnorm(n_subjects, 0, IIV_OMEGA$V1)
  eta_KA  <- rnorm(n_subjects, 0, IIV_OMEGA$KA)

  # Individual parameters = typical * exp(eta)
  tibble(
    ID = 1:n_subjects,
    CL1 = typical_params$CL1 * exp(eta_CL1),
    CL2 = typical_params$CL2 * exp(eta_CL2),
    V1  = typical_params$V1  * exp(eta_V1),
    V2  = typical_params$V2,  # No IIV on V2
    Q   = typical_params$Q,
    KA  = typical_params$KA  * exp(eta_KA),
    F1  = typical_params$F1,
    KDES = typical_params$KDES
  )
}

#-------------------------------------------------------------------------------
# NCA Calculation Function (First Dose Only)
#-------------------------------------------------------------------------------

calculate_nca_first_dose <- function(sim_data, day_stepup2) {
  # Filter data for first dose interval (Day 1 to step-up 2 day)
  # TIME_DAY starts from 1 (Day 1 = time 0)
  conc_data <- sim_data %>%
    filter(TIME_DAY >= 1 & TIME_DAY <= day_stepup2) %>%
    filter(DV > 0) %>%
    select(TIME_DAY, DV)

  if (nrow(conc_data) < 2) {
    return(NULL)
  }

  # Basic NCA parameters
  Cmax <- max(conc_data$DV)
  Tmax <- conc_data$TIME_DAY[which.max(conc_data$DV)]
  Clast <- tail(conc_data$DV, 1)
  Tlast <- tail(conc_data$TIME_DAY, 1)

  # AUC using trapezoidal rule
  AUC_total <- 0
  for (i in 2:nrow(conc_data)) {
    dt <- conc_data$TIME_DAY[i] - conc_data$TIME_DAY[i-1]
    avg_conc <- (conc_data$DV[i] + conc_data$DV[i-1]) / 2
    AUC_total <- AUC_total + dt * avg_conc
  }

  # Calculate terminal half-life (using last 50% of data for first dose)
  n_points <- nrow(conc_data)
  terminal_start <- max(1, floor(n_points * 0.5))
  terminal_data <- conc_data[terminal_start:n_points, ]

  if (nrow(terminal_data) >= 3 && all(terminal_data$DV > 0)) {
    log_conc <- log(terminal_data$DV)
    time_vals <- terminal_data$TIME_DAY

    fit <- lm(log_conc ~ time_vals)
    lambda_z <- -coef(fit)[2]

    if (lambda_z > 0) {
      t_half <- log(2) / lambda_z
      AUC_inf <- AUC_total + Clast / lambda_z
    } else {
      t_half <- NA
      AUC_inf <- NA
      lambda_z <- NA
    }
  } else {
    t_half <- NA
    AUC_inf <- NA
    lambda_z <- NA
  }

  tibble(
    Parameter = c("Cmax", "Tmax", "AUC(first dose)", "AUC(0-inf)", "t1/2", "Lambda_z"),
    Value = c(Cmax, Tmax, AUC_total, AUC_inf, t_half, lambda_z),
    Unit = c("mg/L", "day", "mg·day/L", "mg·day/L", "day", "1/day")
  )
}

#-------------------------------------------------------------------------------
# Calculate NCA Summary Statistics for Monte Carlo
#-------------------------------------------------------------------------------

calculate_nca_summary <- function(all_sim_data, day_stepup2) {
  # Calculate NCA for each subject
  nca_results <- all_sim_data %>%
    group_by(ID) %>%
    group_modify(~ {
      nca <- calculate_nca_first_dose(.x, day_stepup2)
      if (is.null(nca)) {
        tibble(Parameter = character(), Value = numeric(), Unit = character())
      } else {
        nca
      }
    }) %>%
    ungroup()

  if (nrow(nca_results) == 0) return(NULL)

  # Calculate summary statistics
  nca_summary <- nca_results %>%
    group_by(Parameter, Unit) %>%
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      Median = median(Value, na.rm = TRUE),
      Q5 = quantile(Value, 0.05, na.rm = TRUE),
      Q95 = quantile(Value, 0.95, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      `CV%` = (SD / Mean) * 100,
      `90% PI` = paste0("[", signif(Q5, 3), " - ", signif(Q95, 3), "]")
    ) %>%
    select(Parameter, Mean, SD, `CV%`, Median, `90% PI`, Unit)

  nca_summary
}

#-------------------------------------------------------------------------------
# Calculate Cavg up to specific time point (AUC / time)
#-------------------------------------------------------------------------------

calculate_cavg_up_to_time <- function(sim_data, end_hour) {
  # Calculate Cavg from time 0 up to end_hour for each subject
  # Cavg = AUC(0-t) / t
  # TIME_HOUR is the internal simulation time (0-based)

  cavg_data <- sim_data %>%
    filter(TIME_HOUR >= 0 & TIME_HOUR <= end_hour) %>%
    group_by(ID) %>%
    arrange(TIME_HOUR) %>%
    summarise(
      # Calculate AUC using trapezoidal rule
      AUC = {
        auc <- 0
        dv <- DV
        th <- TIME_HOUR
        for (i in 2:length(dv)) {
          dt <- (th[i] - th[i-1]) / 24  # Convert to days
          avg_c <- (dv[i] + dv[i-1]) / 2
          auc <- auc + dt * avg_c
        }
        auc
      },
      Cmax = max(DV),
      .groups = "drop"
    ) %>%
    mutate(
      Cavg = AUC / (end_hour / 24)  # Cavg = AUC / time in days
    )

  if (nrow(cavg_data) == 0) return(list(mean = NA, median = NA, q5 = NA, q95 = NA, cmax_mean = NA, cmax_median = NA))

  list(
    mean = mean(cavg_data$Cavg, na.rm = TRUE),
    median = median(cavg_data$Cavg, na.rm = TRUE),
    q5 = quantile(cavg_data$Cavg, 0.05, na.rm = TRUE),
    q95 = quantile(cavg_data$Cavg, 0.95, na.rm = TRUE),
    cmax_mean = mean(cavg_data$Cmax, na.rm = TRUE),
    cmax_median = median(cavg_data$Cmax, na.rm = TRUE)
  )
}

#-------------------------------------------------------------------------------
# UI Definition
#-------------------------------------------------------------------------------

ui <- fluidPage(

  # Custom CSS
  tags$head(
    tags$style(HTML("
      .main-header {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        padding: 20px;
        margin-bottom: 20px;
        border-radius: 10px;
      }
      .param-box {
        background-color: #f8f9fa;
        border-radius: 10px;
        padding: 15px;
        margin-bottom: 15px;
        border: 1px solid #dee2e6;
      }
      .result-box {
        background-color: #e8f4f8;
        border-radius: 10px;
        padding: 15px;
        border: 1px solid #b8daff;
        margin-bottom: 20px;
      }
      .warning-box {
        background-color: #fff3cd;
        border-radius: 10px;
        padding: 15px;
        border: 2px solid #ffc107;
        margin-bottom: 20px;
      }
      .danger-box {
        background-color: #f8d7da;
        border-radius: 10px;
        padding: 15px;
        border: 2px solid #dc3545;
        margin-bottom: 20px;
      }
      .success-box {
        background-color: #d4edda;
        border-radius: 10px;
        padding: 15px;
        border: 2px solid #28a745;
        margin-bottom: 20px;
      }
      .info-box {
        background-color: #cce5ff;
        border-radius: 10px;
        padding: 15px;
        border: 2px solid #004085;
        margin-bottom: 20px;
      }
      .btn-simulate {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        border: none;
        font-size: 18px;
        padding: 12px 30px;
      }
      .btn-simulate:hover {
        background: linear-gradient(135deg, #764ba2 0%, #667eea 100%);
      }
      .crs-metric {
        font-size: 16px;
        padding: 8px 12px;
        margin: 5px 0;
        border-radius: 5px;
      }
      .crs-warning {
        background-color: #fff3cd;
        border-left: 4px solid #ffc107;
      }
      .crs-danger {
        background-color: #f8d7da;
        border-left: 4px solid #dc3545;
      }
      .crs-safe {
        background-color: #d4edda;
        border-left: 4px solid #28a745;
      }
      .crs-info {
        background-color: #cce5ff;
        border-left: 4px solid #004085;
      }
    "))
  ),

  # Header
  div(class = "main-header",
      h1("Teclistamab PK Simulator", style = "margin: 0;"),
      p("Monte Carlo Simulation with CRS Risk Assessment (Cavg-based)", style = "margin: 5px 0 0 0; opacity: 0.9;")
  ),

  # Main Layout
  fluidRow(
    # Input Panel
    column(4,
           div(class = "param-box",
               h4(icon("user"), " Patient Characteristics"),
               hr(),
               numericInput("bw", "Body Weight (kg):", value = 70, min = 30, max = 150, step = 0.1),
               selectInput("iss", "ISS Stage:", choices = c("I", "II", "III"), selected = "II"),
               selectInput("igg_type", "Immunoglobulin Type:", choices = c("IgG", "Non-IgG"), selected = "IgG")
           ),

           div(class = "param-box",
               h4(icon("syringe"), " Dosing Schedule"),
               hr(),
               p(strong("Step-up Dosing Days:"), style = "margin-bottom: 10px;"),
               fluidRow(
                 column(4, numericInput("day_stepup1", "Step-up 1 (Day):", value = 1, min = 1, max = 10, step = 1)),
                 column(4, numericInput("day_stepup2", "Step-up 2 (Day):", value = 4, min = 2, max = 14, step = 1)),
                 column(4, numericInput("day_treatment", "Treatment (Day):", value = 7, min = 3, max = 21, step = 1))
               ),
               hr(),
               p(strong("Dose Amounts:"), style = "margin-bottom: 10px;"),
               numericInput("dose1", "Step-up 1st Dose (mg/kg):", value = 0.06, min = 0.01, max = 1, step = 0.01),
               numericInput("dose2", "Step-up 2nd Dose (mg/kg):", value = 0.3, min = 0.01, max = 1, step = 0.01),
               numericInput("dose_treat", "Treatment Dose (mg/kg):", value = 1.5, min = 0.1, max = 10, step = 0.1),
               numericInput("n_treatment_doses", "Number of Treatment Doses:", value = 6, min = 1, max = 50, step = 1),
               numericInput("sim_duration", "Simulation Duration after Last Dose (days):", value = 21, min = 7, max = 90, step = 1)
           ),

           div(class = "param-box",
               h4(icon("random"), " Monte Carlo Settings"),
               hr(),
               numericInput("n_subjects", "Number of Virtual Subjects:", value = 100, min = 10, max = 1000, step = 10),
               numericInput("seed", "Random Seed (optional):", value = 12345, min = 1, step = 1),
               p(tags$small("IIV: CL1=53.6%, CL2=107%, V1=48.8%, KA=45.2%"), style = "color: #6c757d;")
           ),

           div(style = "text-align: center; margin-top: 20px;",
               actionButton("simulate", "Run Monte Carlo Simulation",
                            class = "btn btn-primary btn-simulate btn-lg",
                            icon = icon("play"))
           )
    ),

    # Output Panel
    column(8,
           # CRS Risk Assessment
           uiOutput("crs_warning_ui"),

           # Calculated PK Parameters (Typical)
           div(class = "result-box",
               h4(icon("calculator"), " Typical PK Parameters"),
               hr(),
               tableOutput("pk_params_table")
           ),

           # PK Curve
           div(class = "result-box",
               h4(icon("chart-line"), " Time-Concentration Profile (Monte Carlo)"),
               hr(),
               plotOutput("pk_plot", height = "450px")
           ),

           # Cavg at specific timepoints
           div(class = "result-box",
               h4(icon("crosshairs"), " Cavg at Key Timepoints (CRS Risk Metrics)"),
               hr(),
               tableOutput("cavg_table")
           ),

           # NCA Parameters (First Dose)
           div(class = "result-box",
               h4(icon("table"), " NCA Parameters (First Dose)"),
               hr(),
               DTOutput("nca_table")
           )
    )
  ),

  # Footer
  div(style = "text-align: center; margin-top: 30px; padding: 20px; color: #6c757d;",
      p("Model: 2-Compartment with Time-Dependent Clearance"),
      p("Reference: Miao et al. (2023) - Teclistamab Population PK Model")
  )
)

#-------------------------------------------------------------------------------
# Server Logic
#-------------------------------------------------------------------------------

server <- function(input, output, session) {

  # Reactive: Calculate typical PK parameters
  typical_params <- reactive({
    calculate_pk_params_typical(input$bw, input$iss, input$igg_type)
  })

  # Display typical PK parameters
  output$pk_params_table <- renderTable({
    params <- typical_params()

    data.frame(
      Parameter = c("CL1 (Linear CL)", "CL2 (Time-dep CL)", "V1 (Central)",
                    "V2 (Peripheral)", "Q", "KA", "F1", "KDES"),
      `Typical Value` = c(
        sprintf("%.4f", params$CL1),
        sprintf("%.4f", params$CL2),
        sprintf("%.3f", params$V1),
        sprintf("%.4f", params$V2),
        sprintf("%.4f", params$Q),
        sprintf("%.3f", params$KA),
        sprintf("%.3f", params$F1),
        sprintf("%.4f", params$KDES)
      ),
      Unit = c("L/day", "L/day", "L", "L", "L/day", "1/day", "-", "1/day"),
      `IIV (CV%)` = c("53.6%", "107%", "48.8%", "-", "-", "45.2%", "-", "-"),
      check.names = FALSE
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE, width = "100%")

  # Reactive: Run Monte Carlo simulation
  sim_result <- eventReactive(input$simulate, {

    withProgress(message = 'Running Monte Carlo simulation...', value = 0, {

      params_typical <- typical_params()
      bw <- input$bw

      # Generate individual parameters for all subjects
      incProgress(0.1, detail = "Generating individual parameters...")
      ind_params <- generate_individual_params(
        params_typical,
        input$n_subjects,
        seed = input$seed
      )

      # Create dosing schedule based on user input
      # Convert Day to hours (Day 1 = time 0)
      time_stepup1 <- (input$day_stepup1 - 1) * 24  # Day 1 = 0 hours
      time_stepup2 <- (input$day_stepup2 - 1) * 24
      time_treatment_start <- (input$day_treatment - 1) * 24

      dose1_mg <- input$dose1 * bw
      dose2_mg <- input$dose2 * bw
      dose_treat_mg <- input$dose_treat * bw

      # Treatment dose times
      treatment_times <- time_treatment_start + seq(0, (input$n_treatment_doses - 1) * 7 * 24, by = 7 * 24)

      # All dose times and amounts
      dose_times <- c(time_stepup1, time_stepup2, treatment_times)
      dose_amounts <- c(dose1_mg, dose2_mg, rep(dose_treat_mg, input$n_treatment_doses))

      incProgress(0.1, detail = "Creating dosing schedule...")

      # Simulation end time
      max_time <- max(dose_times) + input$sim_duration * 24
      sim_times <- seq(0, max_time, by = 1)  # hourly

      # Run simulation for each subject
      incProgress(0.1, detail = "Running simulations...")

      all_results <- map_dfr(1:input$n_subjects, function(i) {

        if (i %% 10 == 0) {
          incProgress(0.5 / input$n_subjects * 10,
                      detail = paste0("Simulating subject ", i, "/", input$n_subjects))
        }

        # Get individual parameters
        ind_p <- ind_params %>% filter(ID == i)

        # Dosing data for this subject
        dosing_data <- tibble(
          ID = i,
          time = dose_times,
          amt = dose_amounts,
          cmt = 1,
          evid = 1
        )

        # Observation data
        obs_data <- tibble(
          ID = i,
          time = sim_times,
          amt = 0,
          cmt = 0,
          evid = 0
        )

        sim_data <- bind_rows(dosing_data, obs_data) %>%
          arrange(time, desc(evid))

        # Update model parameters
        mod_i <- mod %>%
          param(
            CL1 = ind_p$CL1,
            CL2 = ind_p$CL2,
            KDES = ind_p$KDES,
            V1 = ind_p$V1,
            V2 = ind_p$V2,
            Q = ind_p$Q,
            KA = ind_p$KA,
            F1 = ind_p$F1
          )

        # Run simulation
        out <- mod_i %>%
          data_set(sim_data) %>%
          mrgsim(carry_out = "amt,evid") %>%
          as_tibble() %>%
          filter(evid == 0) %>%
          mutate(
            TIME_HOUR = time,
            TIME_DAY = time / 24 + 1  # Day 1 starts at time 0
          )

        out
      })

      incProgress(0.2, detail = "Complete!")

      # Dosing info for plot
      dosing_info <- tibble(
        time_day = dose_times / 24 + 1,
        dose_mg = dose_amounts,
        dose_type = c("Step-up 1", "Step-up 2", rep("Treatment", input$n_treatment_doses))
      )

      # Calculate Cavg up to Day 3 (72hr) and Day 5 (120hr)
      # Day 3 Cavg = AUC(0-72hr) / 3 days
      # Day 5 Cavg = AUC(0-120hr) / 5 days
      cavg_day3 <- calculate_cavg_up_to_time(all_results, 72)   # Cavg from 0 to 72hr
      cavg_day5 <- calculate_cavg_up_to_time(all_results, 120)  # Cavg from 0 to 120hr

      list(
        simulation = all_results,
        dosing_info = dosing_info,
        ind_params = ind_params,
        typical_params = params_typical,
        cavg_day3 = cavg_day3,
        cavg_day5 = cavg_day5,
        bw = bw
      )
    })
  })

  # CRS Warning UI
  output$crs_warning_ui <- renderUI({
    if (is.null(input$simulate) || input$simulate == 0) return(NULL)
    result <- sim_result()
    if (is.null(result)) return(NULL)

    cavg_day3 <- result$cavg_day3$median
    cavg_day5 <- result$cavg_day5$median

    warnings <- list()
    infos <- list()

    # Check Day 3 Cavg >= 0.18 (CRS Grade II risk)
    if (!is.na(cavg_day3) && cavg_day3 >= 0.18) {
      warnings <- c(warnings, list(
        div(class = "crs-metric crs-danger",
            icon("exclamation-triangle"),
            strong(" CRS Grade II 위험 높음! "),
            sprintf("Day 3 Cavg (%.4f µg/mL) ≥ 0.18 µg/mL", cavg_day3)
        )
      ))
    } else if (!is.na(cavg_day3)) {
      infos <- c(infos, list(
        div(class = "crs-metric crs-safe",
            icon("check"),
            sprintf(" Day 3 Cavg (%.4f µg/mL) < 0.18 µg/mL - 정상 범위", cavg_day3)
        )
      ))
    }

    # Check Day 5 Cavg
    if (!is.na(cavg_day5)) {
      if (cavg_day5 < 0.32) {
        # Low efficacy warning
        warnings <- c(warnings, list(
          div(class = "crs-metric crs-warning",
              icon("exclamation-circle"),
              strong(" Efficacy 저하 우려 "),
              sprintf("Day 5 Cavg (%.4f µg/mL) < 0.32 µg/mL", cavg_day5)
          )
        ))
      } else if (cavg_day5 > 0.45) {
        # CRS risk warning
        warnings <- c(warnings, list(
          div(class = "crs-metric crs-danger",
              icon("exclamation-triangle"),
              strong(" CRS 위험 높음! "),
              sprintf("Day 5 Cavg (%.4f µg/mL) > 0.45 µg/mL", cavg_day5)
          )
        ))
      } else {
        # Optimal range
        infos <- c(infos, list(
          div(class = "crs-metric crs-safe",
              icon("check-circle"),
              strong(" Optimal 범위! "),
              sprintf("Day 5 Cavg (%.4f µg/mL) ∈ [0.32 - 0.45] µg/mL", cavg_day5)
          )
        ))
      }
    }

    if (length(warnings) == 0) {
      div(class = "success-box",
          h4(icon("check-circle"), " CRS Risk Assessment"),
          hr(),
          infos,
          div(class = "crs-metric crs-safe",
              icon("thumbs-up"),
              " 모든 CRS 위험 지표 정상 범위"
          )
      )
    } else {
      div(class = "danger-box",
          h4(icon("exclamation-triangle"), " CRS Risk Assessment"),
          hr(),
          warnings,
          if (length(infos) > 0) infos
      )
    }
  })

  # Cavg at specific timepoints table
  output$cavg_table <- renderTable({
    if (is.null(input$simulate) || input$simulate == 0) return(NULL)
    result <- sim_result()
    if (is.null(result)) return(NULL)

    cavg_day3 <- result$cavg_day3
    cavg_day5 <- result$cavg_day5

    data.frame(
      Parameter = c("Day 3 Cavg (0-72 hr)", "Day 5 Cavg (0-120 hr)"),
      `Median (µg/mL)` = c(
        sprintf("%.4f", cavg_day3$median),
        sprintf("%.4f", cavg_day5$median)
      ),
      `Mean (µg/mL)` = c(
        sprintf("%.4f", cavg_day3$mean),
        sprintf("%.4f", cavg_day5$mean)
      ),
      `90% PI` = c(
        sprintf("[%.4f - %.4f]", cavg_day3$q5, cavg_day3$q95),
        sprintf("[%.4f - %.4f]", cavg_day5$q5, cavg_day5$q95)
      ),
      `Optimal Range` = c("< 0.18 µg/mL", "0.32 - 0.45 µg/mL"),
      `Risk` = c("≥0.18: CRS Gr.II↑", "<0.32: Efficacy↓, >0.45: CRS↑"),
      check.names = FALSE
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE, width = "100%")

  # Plot PK curves
  output$pk_plot <- renderPlot({
    if (is.null(input$simulate) || input$simulate == 0) return(NULL)
    result <- sim_result()
    if (is.null(result)) return(NULL)

    sim_data <- result$simulation
    dosing_info <- result$dosing_info

    # Calculate summary statistics (convert TIME_DAY back to 0-based for plotting)
    summary_data <- sim_data %>%
      mutate(TIME_PLOT = TIME_DAY - 1) %>%
      group_by(TIME_PLOT) %>%
      summarise(
        median = median(DV),
        q5 = quantile(DV, 0.05),
        q25 = quantile(DV, 0.25),
        q75 = quantile(DV, 0.75),
        q95 = quantile(DV, 0.95),
        .groups = "drop"
      )

    # Convert dosing times to 0-based for plotting
    dosing_plot <- dosing_info %>%
      mutate(time_plot = time_day - 1)

    # Plot
    ggplot() +
      # 90% prediction interval
      geom_ribbon(data = summary_data,
                  aes(x = TIME_PLOT, ymin = q5, ymax = q95),
                  fill = "#667eea", alpha = 0.2) +
      # 50% prediction interval
      geom_ribbon(data = summary_data,
                  aes(x = TIME_PLOT, ymin = q25, ymax = q75),
                  fill = "#667eea", alpha = 0.3) +
      # Median line
      geom_line(data = summary_data,
                aes(x = TIME_PLOT, y = median),
                color = "#667eea", linewidth = 1.2) +
      # Dose markers
      geom_vline(data = dosing_plot,
                 aes(xintercept = time_plot, color = dose_type),
                 linetype = "dashed", alpha = 0.7) +
      # Day 3 and Day 5 markers
      geom_vline(xintercept = 3, linetype = "dotted", color = "orange", linewidth = 0.8) +
      geom_vline(xintercept = 5, linetype = "dotted", color = "red", linewidth = 0.8) +
      annotate("text", x = 3.2, y = 0.002, label = "72hr", color = "orange", hjust = 0, size = 3) +
      annotate("text", x = 5.2, y = 0.002, label = "120hr", color = "red", hjust = 0, size = 3) +
      # Threshold lines for Cavg reference
      geom_hline(yintercept = 0.18, linetype = "dashed", color = "orange", alpha = 0.7) +
      geom_hline(yintercept = 0.32, linetype = "dashed", color = "blue", alpha = 0.5) +
      geom_hline(yintercept = 0.45, linetype = "dashed", color = "red", alpha = 0.7) +
      annotate("text", x = max(summary_data$TIME_PLOT) - 5, y = 0.20,
               label = "0.18 (Day3 Cavg)", color = "orange", size = 3) +
      annotate("text", x = max(summary_data$TIME_PLOT) - 5, y = 0.35,
               label = "0.32", color = "blue", size = 3) +
      annotate("text", x = max(summary_data$TIME_PLOT) - 5, y = 0.50,
               label = "0.45 (Day5 Cavg)", color = "red", size = 3) +
      scale_color_manual(values = c("Step-up 1" = "#e74c3c",
                                    "Step-up 2" = "#f39c12",
                                    "Treatment" = "#27ae60"),
                         name = "Dose Type") +
      scale_y_log10(
        breaks = c(0.001, 0.01, 0.1, 0.18, 0.32, 0.45, 1, 10, 100),
        labels = c("0.001", "0.01", "0.1", "0.18", "0.32", "0.45", "1", "10", "100")
      ) +
      annotation_logticks(sides = "l") +
      labs(
        x = "Time (days)",
        y = "Concentration (µg/mL)",
        title = paste0("Teclistamab PK Profile (n=", input$n_subjects, " subjects)"),
        subtitle = paste0("BW: ", input$bw, " kg | ISS: ", input$iss,
                          " | ", input$igg_type,
                          "\nShaded: 90% PI (light) and 50% PI (dark), Line: Median")
      ) +
      theme_bw(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 11, color = "gray40"),
        panel.grid.minor = element_line(color = "gray90"),
        legend.position = "bottom"
      ) +
      coord_cartesian(ylim = c(0.001, NA), xlim = c(0, NA))
  })

  # NCA table (First Dose)
  output$nca_table <- renderDT({
    if (is.null(input$simulate) || input$simulate == 0) return(NULL)
    result <- sim_result()
    if (is.null(result)) return(NULL)

    nca_summary <- calculate_nca_summary(result$simulation, input$day_stepup2)

    if (is.null(nca_summary)) {
      return(NULL)
    }

    nca_summary %>%
      mutate(
        Mean = signif(Mean, 4),
        SD = signif(SD, 3),
        `CV%` = round(`CV%`, 1),
        Median = signif(Median, 4)
      ) %>%
      datatable(
        options = list(
          dom = 't',
          pageLength = 10,
          ordering = FALSE
        ),
        rownames = FALSE,
        class = 'cell-border stripe',
        caption = htmltools::tags$caption(
          style = 'caption-side: top; text-align: left; color: gray;',
          paste0('NCA calculated for first dose interval (Day 1 to Day ', input$day_stepup2, ') across ',
                 input$n_subjects, ' virtual subjects')
        )
      )
  })
}

#-------------------------------------------------------------------------------
# Run Application
#-------------------------------------------------------------------------------

shinyApp(ui = ui, server = server)
