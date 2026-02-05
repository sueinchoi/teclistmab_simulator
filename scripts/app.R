#===============================================================================
# Teclistamab PK Simulation Shiny App with AI Dosing Optimizer
#
# Features:
# 1. Individual PK parameter calculation based on patient characteristics
# 2. 2-Compartment model with time-dependent clearance
# 3. Monte Carlo simulation with inter-individual variability (IIV)
# 4. Risk Assessment based on Cavg (Response & CRS)
# 5. AI Dosing Optimizer Agent - suggests optimal dosing regimens
#
# Reference: Miao et al. (2023) - Teclistamab Population PK Model
# Risk Cut-offs: ROC analysis (Response: Youden, CRS: Accuracy)
#===============================================================================

library(shiny)
library(mrgsolve)
library(tidyverse)
library(DT)
library(shinycssloaders)

#-------------------------------------------------------------------------------
# Optimal Cut-off Values (from ROC analysis)
#-------------------------------------------------------------------------------

# Response (>=VGPR): Youden's Index based
RESPONSE_LOWER_CUTOFF <- 0.296  # Below this: low efficacy risk

# CRS Grade 2+: Accuracy based
CRS_UPPER_CUTOFF <- 0.350  # Above this: high CRS risk

# Optimal therapeutic window
OPTIMAL_LOWER <- RESPONSE_LOWER_CUTOFF
OPTIMAL_UPPER <- CRS_UPPER_CUTOFF

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
#-------------------------------------------------------------------------------

cv_to_omega <- function(cv) {
  sqrt(log(1 + cv^2))
}

IIV_CV <- list(
  CL1 = 0.536,
  CL2 = 1.07,
  V1  = 0.488,
  KA  = 0.452
)

IIV_OMEGA <- list(
  CL1 = cv_to_omega(IIV_CV$CL1),
  CL2 = cv_to_omega(IIV_CV$CL2),
  V1  = cv_to_omega(IIV_CV$V1),
  KA  = cv_to_omega(IIV_CV$KA)
)

#-------------------------------------------------------------------------------
# PK Parameter Calculation Functions
#-------------------------------------------------------------------------------

calculate_pk_params_typical <- function(bw, iss, igg_type) {
  iss_2 <- ifelse(iss == "II", 1, 0)
  iss_3 <- ifelse(iss == "III", 1, 0)
  non_igg <- ifelse(igg_type == "Non-IgG", 1, 0)

  CL1 <- 0.449 * (bw / 74)^0.704 * (1.31^iss_2) * (1.67^iss_3) * (0.689^non_igg)
  CL2 <- 0.547 * (0.295^non_igg)
  V1 <- 4.13 * (bw / 74)^0.358
  V2 <- 1.34 * (bw / 74)^1.40
  Q <- 0.039
  KA <- 0.133
  F1 <- 0.718
  KDES <- 0.0292

  list(CL1 = CL1, CL2 = CL2, V1 = V1, V2 = V2, Q = Q, KA = KA, F1 = F1, KDES = KDES)
}

generate_individual_params <- function(typical_params, n_subjects, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  eta_CL1 <- rnorm(n_subjects, 0, IIV_OMEGA$CL1)
  eta_CL2 <- rnorm(n_subjects, 0, IIV_OMEGA$CL2)
  eta_V1  <- rnorm(n_subjects, 0, IIV_OMEGA$V1)
  eta_KA  <- rnorm(n_subjects, 0, IIV_OMEGA$KA)

  tibble(
    ID = 1:n_subjects,
    CL1 = typical_params$CL1 * exp(eta_CL1),
    CL2 = typical_params$CL2 * exp(eta_CL2),
    V1  = typical_params$V1  * exp(eta_V1),
    V2  = typical_params$V2,
    Q   = typical_params$Q,
    KA  = typical_params$KA  * exp(eta_KA),
    F1  = typical_params$F1,
    KDES = typical_params$KDES
  )
}

#-------------------------------------------------------------------------------
# Simulation Helper Functions
#-------------------------------------------------------------------------------

calculate_cavg_up_to_time <- function(sim_data, end_hour) {
  cavg_data <- sim_data %>%
    filter(TIME_HOUR >= 0 & TIME_HOUR <= end_hour) %>%
    group_by(ID) %>%
    arrange(TIME_HOUR) %>%
    summarise(
      AUC = {
        auc <- 0
        dv <- DV
        th <- TIME_HOUR
        for (i in 2:length(dv)) {
          dt <- (th[i] - th[i-1]) / 24
          avg_c <- (dv[i] + dv[i-1]) / 2
          auc <- auc + dt * avg_c
        }
        auc
      },
      Cmax = max(DV),
      .groups = "drop"
    ) %>%
    mutate(Cavg = AUC / (end_hour / 24))

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
# Core Simulation Function (for single scenario)
#-------------------------------------------------------------------------------

run_single_simulation <- function(mod, typical_params, bw,
                                   dose1, dose2, dose_treat,
                                   day_stepup1, day_stepup2, day_treatment,
                                   n_treatment_doses, sim_duration,
                                   n_subjects, seed = NULL) {

  ind_params <- generate_individual_params(typical_params, n_subjects, seed = seed)

  time_stepup1 <- (day_stepup1 - 1) * 24
  time_stepup2 <- (day_stepup2 - 1) * 24
  time_treatment_start <- (day_treatment - 1) * 24

  dose1_mg <- dose1 * bw
  dose2_mg <- dose2 * bw
  dose_treat_mg <- dose_treat * bw

  treatment_times <- time_treatment_start + seq(0, (n_treatment_doses - 1) * 7 * 24, by = 7 * 24)
  dose_times <- c(time_stepup1, time_stepup2, treatment_times)
  dose_amounts <- c(dose1_mg, dose2_mg, rep(dose_treat_mg, n_treatment_doses))

  max_time <- max(dose_times) + sim_duration * 24
  sim_times <- seq(0, max_time, by = 2)  # every 2 hours for speed

  all_results <- map_dfr(1:n_subjects, function(i) {
    ind_p <- ind_params %>% filter(ID == i)

    dosing_data <- tibble(
      ID = i, time = dose_times, amt = dose_amounts, cmt = 1, evid = 1
    )

    obs_data <- tibble(
      ID = i, time = sim_times, amt = 0, cmt = 0, evid = 0
    )

    sim_data <- bind_rows(dosing_data, obs_data) %>%
      arrange(time, desc(evid))

    mod_i <- mod %>%
      param(CL1 = ind_p$CL1, CL2 = ind_p$CL2, KDES = ind_p$KDES,
            V1 = ind_p$V1, V2 = ind_p$V2, Q = ind_p$Q,
            KA = ind_p$KA, F1 = ind_p$F1)

    out <- mod_i %>%
      data_set(sim_data) %>%
      mrgsim(carry_out = "amt,evid") %>%
      as_tibble() %>%
      filter(evid == 0) %>%
      mutate(TIME_HOUR = time, TIME_DAY = time / 24 + 1)

    out
  })

  cavg_day3 <- calculate_cavg_up_to_time(all_results, 72)
  cavg_day5 <- calculate_cavg_up_to_time(all_results, 120)

  list(
    simulation = all_results,
    cavg_day3 = cavg_day3,
    cavg_day5 = cavg_day5,
    params = list(
      dose1 = dose1, dose2 = dose2, dose_treat = dose_treat,
      day_stepup1 = day_stepup1, day_stepup2 = day_stepup2,
      day_treatment = day_treatment, n_treatment_doses = n_treatment_doses
    )
  )
}

#-------------------------------------------------------------------------------
# AI Dosing Optimizer - Parse Request and Generate Scenarios
#-------------------------------------------------------------------------------

parse_ai_request <- function(request_text) {
  # Parse natural language request and determine optimization goal
  request_lower <- tolower(request_text)

  goals <- list(
    reduce_crs = FALSE,
    improve_efficacy = FALSE,
    find_optimal = FALSE,
    reduce_dose = FALSE,
    extend_interval = FALSE
  )

  # Korean keywords
  if (grepl("crs|부작용|안전|위험.*줄|낮", request_lower)) {
    goals$reduce_crs <- TRUE
  }
  if (grepl("효과|반응|vgpr|efficacy|response", request_lower)) {
    goals$improve_efficacy <- TRUE
  }
  if (grepl("최적|optimal|best|추천|권장", request_lower)) {
    goals$find_optimal <- TRUE
  }
  if (grepl("용량.*줄|낮은.*용량|dose.*reduc|lower.*dose", request_lower)) {
    goals$reduce_dose <- TRUE
  }
  if (grepl("간격|interval|늘리|연장|extend", request_lower)) {
    goals$extend_interval <- TRUE
  }

  # Default: find optimal if no specific goal
  if (!any(unlist(goals))) {
    goals$find_optimal <- TRUE
  }

  goals
}

generate_dosing_scenarios <- function(base_params, goals) {
  # Generate alternative dosing scenarios based on optimization goals

  scenarios <- list()

  # Current (baseline)
  scenarios[["Current"]] <- base_params

  if (goals$reduce_crs || goals$find_optimal) {
    # Lower step-up doses
    scenarios[["Step-up 50%"]] <- modifyList(base_params, list(
      dose1 = base_params$dose1 * 0.5,
      dose2 = base_params$dose2 * 0.5
    ))

    # Extend step-up interval
    scenarios[["Extended Step-up (D1,5,8)"]] <- modifyList(base_params, list(
      day_stepup2 = 5,
      day_treatment = 8
    ))

    # Lower treatment dose
    scenarios[["Treatment 80%"]] <- modifyList(base_params, list(
      dose_treat = base_params$dose_treat * 0.8
    ))
  }

  if (goals$improve_efficacy || goals$find_optimal) {
    # Higher step-up doses
    scenarios[["Step-up 150%"]] <- modifyList(base_params, list(
      dose1 = base_params$dose1 * 1.5,
      dose2 = base_params$dose2 * 1.5
    ))

    # Higher treatment dose
    scenarios[["Treatment 120%"]] <- modifyList(base_params, list(
      dose_treat = base_params$dose_treat * 1.2
    ))

    # Shorter step-up interval
    scenarios[["Shorter Step-up (D1,3,5)"]] <- modifyList(base_params, list(
      day_stepup2 = 3,
      day_treatment = 5
    ))
  }

  if (goals$reduce_dose) {
    # Various dose reductions
    scenarios[["All Doses 70%"]] <- modifyList(base_params, list(
      dose1 = base_params$dose1 * 0.7,
      dose2 = base_params$dose2 * 0.7,
      dose_treat = base_params$dose_treat * 0.7
    ))

    scenarios[["Step-up Only 50%"]] <- modifyList(base_params, list(
      dose1 = base_params$dose1 * 0.5,
      dose2 = base_params$dose2 * 0.5
    ))
  }

  if (goals$extend_interval) {
    # Extended intervals
    scenarios[["Step-up D1,5,9"]] <- modifyList(base_params, list(
      day_stepup2 = 5,
      day_treatment = 9
    ))

    scenarios[["Step-up D1,6,11"]] <- modifyList(base_params, list(
      day_stepup2 = 6,
      day_treatment = 11
    ))
  }

  # Additional balanced scenarios
  scenarios[["Balanced 1"]] <- modifyList(base_params, list(
    dose1 = base_params$dose1 * 0.75,
    dose2 = base_params$dose2 * 0.75,
    day_stepup2 = 5
  ))

  scenarios[["Balanced 2"]] <- modifyList(base_params, list(
    dose2 = base_params$dose2 * 0.8,
    day_stepup2 = 5,
    day_treatment = 8
  ))

  scenarios
}

score_scenario <- function(cavg_day5, goals) {
  # Score a scenario based on Cavg Day 5 and goals
  # Higher score = better

  if (is.na(cavg_day5)) return(-Inf)

  score <- 0

  # In optimal range: +100
  if (cavg_day5 >= OPTIMAL_LOWER && cavg_day5 <= OPTIMAL_UPPER) {
    score <- score + 100
  }

  # Distance penalty from optimal range
  if (cavg_day5 < OPTIMAL_LOWER) {
    score <- score - (OPTIMAL_LOWER - cavg_day5) * 200  # Efficacy penalty
  }
  if (cavg_day5 > OPTIMAL_UPPER) {
    score <- score - (cavg_day5 - OPTIMAL_UPPER) * 300  # CRS penalty (higher weight)
  }

  # Goal-specific adjustments
  if (goals$reduce_crs) {
    # Prefer lower Cavg
    score <- score - cavg_day5 * 50
  }
  if (goals$improve_efficacy) {
    # Prefer higher Cavg (up to upper limit)
    if (cavg_day5 <= OPTIMAL_UPPER) {
      score <- score + cavg_day5 * 30
    }
  }

  score
}

classify_risk <- function(cavg_day5) {
  if (is.na(cavg_day5)) return("Unknown")

  if (cavg_day5 < OPTIMAL_LOWER) {
    return("Low Efficacy Risk")
  } else if (cavg_day5 > OPTIMAL_UPPER) {
    return("High CRS Risk")
  } else {
    return("Optimal")
  }
}

#-------------------------------------------------------------------------------
# UI Definition
#-------------------------------------------------------------------------------

ui <- fluidPage(

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
      .ai-box {
        background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
        border-radius: 10px;
        padding: 15px;
        border: 2px solid #667eea;
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
      .btn-simulate {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        border: none;
        font-size: 16px;
        padding: 10px 25px;
      }
      .btn-simulate:hover {
        background: linear-gradient(135deg, #764ba2 0%, #667eea 100%);
      }
      .btn-ai {
        background: linear-gradient(135deg, #11998e 0%, #38ef7d 100%);
        border: none;
        font-size: 16px;
        padding: 10px 25px;
        color: white;
      }
      .btn-ai:hover {
        background: linear-gradient(135deg, #38ef7d 0%, #11998e 100%);
        color: white;
      }
      .risk-optimal {
        background-color: #d4edda;
        border-left: 5px solid #28a745;
        padding: 15px;
        margin: 10px 0;
        border-radius: 5px;
      }
      .risk-low-efficacy {
        background-color: #fff3cd;
        border-left: 5px solid #ffc107;
        padding: 15px;
        margin: 10px 0;
        border-radius: 5px;
      }
      .risk-high-crs {
        background-color: #f8d7da;
        border-left: 5px solid #dc3545;
        padding: 15px;
        margin: 10px 0;
        border-radius: 5px;
      }
      .scenario-table .optimal-row {
        background-color: #d4edda !important;
      }
      .recommendation-box {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        padding: 20px;
        border-radius: 10px;
        margin: 15px 0;
      }
    "))
  ),

  # Header
  div(class = "main-header",
      h1("Teclistamab PK Simulator", style = "margin: 0;"),
      p("Monte Carlo Simulation with AI Dosing Optimizer", style = "margin: 5px 0 0 0; opacity: 0.9;"),
      p(sprintf("Optimal Cavg Day 5: %.3f - %.3f mg/L", OPTIMAL_LOWER, OPTIMAL_UPPER),
        style = "margin: 0; font-size: 12px; opacity: 0.8;")
  ),

  # Main Layout with Tabs
  tabsetPanel(
    #---------------------------------------------------------------------------
    # Tab 1: Single Simulation
    #---------------------------------------------------------------------------
    tabPanel("Single Simulation",
             fluidRow(
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
                          p(strong("Step-up Dosing Days:")),
                          fluidRow(
                            column(4, numericInput("day_stepup1", "Step-up 1:", value = 1, min = 1, max = 10)),
                            column(4, numericInput("day_stepup2", "Step-up 2:", value = 4, min = 2, max = 14)),
                            column(4, numericInput("day_treatment", "Treatment:", value = 7, min = 3, max = 21))
                          ),
                          hr(),
                          p(strong("Dose Amounts (mg/kg):")),
                          numericInput("dose1", "Step-up 1st:", value = 0.06, min = 0.01, max = 1, step = 0.01),
                          numericInput("dose2", "Step-up 2nd:", value = 0.3, min = 0.01, max = 1, step = 0.01),
                          numericInput("dose_treat", "Treatment:", value = 1.5, min = 0.1, max = 10, step = 0.1),
                          numericInput("n_treatment_doses", "# Treatment Doses:", value = 6, min = 1, max = 50),
                          numericInput("sim_duration", "Simulation Duration (days after last dose):", value = 21, min = 7, max = 90)
                      ),

                      div(class = "param-box",
                          h4(icon("random"), " Monte Carlo"),
                          hr(),
                          numericInput("n_subjects", "Virtual Subjects:", value = 100, min = 10, max = 1000, step = 10),
                          numericInput("seed", "Random Seed:", value = 12345, min = 1)
                      ),

                      div(style = "text-align: center;",
                          actionButton("simulate", "Run Simulation", class = "btn btn-primary btn-simulate", icon = icon("play"))
                      )
               ),

               column(8,
                      uiOutput("risk_assessment_ui"),
                      div(class = "result-box",
                          h4(icon("chart-line"), " PK Profile"),
                          plotOutput("pk_plot", height = "400px") %>% withSpinner()
                      ),
                      div(class = "result-box",
                          h4(icon("table"), " Cavg Summary"),
                          tableOutput("cavg_table")
                      )
               )
             )
    ),

    #---------------------------------------------------------------------------
    # Tab 2: AI Dosing Optimizer
    #---------------------------------------------------------------------------
    tabPanel("AI Dosing Optimizer",
             fluidRow(
               column(4,
                      div(class = "ai-box",
                          h4(icon("robot"), " AI Dosing Optimizer"),
                          hr(),
                          p("Describe what you want to achieve:"),
                          textAreaInput("ai_request",
                                        label = NULL,
                                        placeholder = "Examples:\n- CRS 위험을 줄이고 싶어요\n- 효과를 유지하면서 안전한 용량을 찾아줘\n- Find optimal dosing for this patient\n- Step-up 간격을 늘려서 시뮬레이션 해줘",
                                        rows = 5,
                                        width = "100%"),
                          hr(),
                          p(strong("Current Patient:")),
                          verbatimTextOutput("ai_patient_summary"),
                          hr(),
                          numericInput("ai_n_subjects", "Subjects per scenario:", value = 50, min = 20, max = 200),
                          div(style = "text-align: center; margin-top: 15px;",
                              actionButton("run_ai", "Optimize Dosing", class = "btn btn-ai btn-lg", icon = icon("magic"))
                          )
                      ),

                      div(class = "param-box",
                          h4(icon("info-circle"), " Reference Cut-offs"),
                          hr(),
                          p(strong("Optimal Range (Cavg Day 5):")),
                          p(sprintf("%.3f - %.3f mg/L", OPTIMAL_LOWER, OPTIMAL_UPPER)),
                          hr(),
                          p(tags$small("Response (>=VGPR): Youden's Index")),
                          p(tags$small("CRS Grade 2+: Accuracy-based")),
                          hr(),
                          tags$ul(
                            tags$li(tags$span(style = "color: #ffc107;", "Below optimal: Low efficacy risk")),
                            tags$li(tags$span(style = "color: #28a745;", "In range: Optimal")),
                            tags$li(tags$span(style = "color: #dc3545;", "Above optimal: High CRS risk"))
                          )
                      )
               ),

               column(8,
                      uiOutput("ai_recommendation_ui"),
                      div(class = "result-box",
                          h4(icon("table"), " Scenario Comparison"),
                          DTOutput("ai_scenario_table") %>% withSpinner()
                      ),
                      div(class = "result-box",
                          h4(icon("chart-bar"), " Cavg Day 5 by Scenario"),
                          plotOutput("ai_comparison_plot", height = "350px") %>% withSpinner()
                      )
               )
             )
    )
  ),

  # Footer
  div(style = "text-align: center; margin-top: 30px; padding: 20px; color: #6c757d;",
      p("Model: 2-Compartment with Time-Dependent Clearance | Reference: Miao et al. (2023)")
  )
)

#-------------------------------------------------------------------------------
# Server Logic
#-------------------------------------------------------------------------------

server <- function(input, output, session) {

  # Reactive: Typical PK parameters
  typical_params <- reactive({
    calculate_pk_params_typical(input$bw, input$iss, input$igg_type)
  })

  # Reactive: Run single simulation
  sim_result <- eventReactive(input$simulate, {
    withProgress(message = 'Running simulation...', value = 0, {

      params_typical <- typical_params()
      incProgress(0.2)

      result <- run_single_simulation(
        mod = mod,
        typical_params = params_typical,
        bw = input$bw,
        dose1 = input$dose1,
        dose2 = input$dose2,
        dose_treat = input$dose_treat,
        day_stepup1 = input$day_stepup1,
        day_stepup2 = input$day_stepup2,
        day_treatment = input$day_treatment,
        n_treatment_doses = input$n_treatment_doses,
        sim_duration = input$sim_duration,
        n_subjects = input$n_subjects,
        seed = input$seed
      )

      incProgress(0.8)
      result
    })
  })

  # Risk Assessment UI
  output$risk_assessment_ui <- renderUI({
    if (is.null(input$simulate) || input$simulate == 0) return(NULL)
    result <- sim_result()
    if (is.null(result)) return(NULL)

    cavg_day5 <- result$cavg_day5$median
    risk_class <- classify_risk(cavg_day5)

    if (risk_class == "Optimal") {
      div(class = "risk-optimal",
          h4(icon("check-circle"), " Risk Assessment: OPTIMAL"),
          p(sprintf("Cavg Day 5: %.4f mg/L", cavg_day5)),
          p(sprintf("Within optimal range (%.3f - %.3f mg/L)", OPTIMAL_LOWER, OPTIMAL_UPPER)),
          p(icon("thumbs-up"), " Good balance between efficacy and safety")
      )
    } else if (risk_class == "Low Efficacy Risk") {
      div(class = "risk-low-efficacy",
          h4(icon("exclamation-triangle"), " Risk Assessment: LOW EFFICACY RISK"),
          p(sprintf("Cavg Day 5: %.4f mg/L", cavg_day5)),
          p(sprintf("Below optimal threshold (%.3f mg/L)", OPTIMAL_LOWER)),
          p(icon("lightbulb"), " Consider: Higher doses or shorter step-up interval")
      )
    } else {
      div(class = "risk-high-crs",
          h4(icon("exclamation-circle"), " Risk Assessment: HIGH CRS RISK"),
          p(sprintf("Cavg Day 5: %.4f mg/L", cavg_day5)),
          p(sprintf("Above safety threshold (%.3f mg/L)", OPTIMAL_UPPER)),
          p(icon("lightbulb"), " Consider: Lower doses or extended step-up interval")
      )
    }
  })

  # Cavg Table
  output$cavg_table <- renderTable({
    if (is.null(input$simulate) || input$simulate == 0) return(NULL)
    result <- sim_result()
    if (is.null(result)) return(NULL)

    cavg_day3 <- result$cavg_day3
    cavg_day5 <- result$cavg_day5

    data.frame(
      Parameter = c("Cavg Day 3 (0-72 hr)", "Cavg Day 5 (0-120 hr)"),
      `Median` = c(sprintf("%.4f", cavg_day3$median), sprintf("%.4f", cavg_day5$median)),
      `Mean` = c(sprintf("%.4f", cavg_day3$mean), sprintf("%.4f", cavg_day5$mean)),
      `90% PI` = c(
        sprintf("[%.4f - %.4f]", cavg_day3$q5, cavg_day3$q95),
        sprintf("[%.4f - %.4f]", cavg_day5$q5, cavg_day5$q95)
      ),
      `Optimal Range` = c("-", sprintf("%.3f - %.3f", OPTIMAL_LOWER, OPTIMAL_UPPER)),
      check.names = FALSE
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE, width = "100%")

  # PK Plot
  output$pk_plot <- renderPlot({
    if (is.null(input$simulate) || input$simulate == 0) return(NULL)
    result <- sim_result()
    if (is.null(result)) return(NULL)

    sim_data <- result$simulation

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

    ggplot(summary_data) +
      geom_ribbon(aes(x = TIME_PLOT, ymin = q5, ymax = q95), fill = "#667eea", alpha = 0.2) +
      geom_ribbon(aes(x = TIME_PLOT, ymin = q25, ymax = q75), fill = "#667eea", alpha = 0.3) +
      geom_line(aes(x = TIME_PLOT, y = median), color = "#667eea", linewidth = 1.2) +
      geom_vline(xintercept = 5, linetype = "dotted", color = "red", linewidth = 0.8) +
      geom_hline(yintercept = OPTIMAL_LOWER, linetype = "dashed", color = "#ffc107", alpha = 0.8) +
      geom_hline(yintercept = OPTIMAL_UPPER, linetype = "dashed", color = "#dc3545", alpha = 0.8) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = OPTIMAL_LOWER, ymax = OPTIMAL_UPPER,
               fill = "#28a745", alpha = 0.1) +
      annotate("text", x = 5.5, y = 0.002, label = "Day 5", color = "red", hjust = 0, size = 3) +
      annotate("text", x = max(summary_data$TIME_PLOT) - 5, y = OPTIMAL_LOWER * 0.8,
               label = sprintf("%.3f", OPTIMAL_LOWER), color = "#ffc107", size = 3) +
      annotate("text", x = max(summary_data$TIME_PLOT) - 5, y = OPTIMAL_UPPER * 1.2,
               label = sprintf("%.3f", OPTIMAL_UPPER), color = "#dc3545", size = 3) +
      scale_y_log10(breaks = c(0.001, 0.01, 0.1, 0.3, 0.5, 1, 10, 100),
                    labels = c("0.001", "0.01", "0.1", "0.3", "0.5", "1", "10", "100")) +
      annotation_logticks(sides = "l") +
      labs(x = "Time (days)", y = "Concentration (mg/L)",
           title = paste0("Teclistamab PK Profile (n=", input$n_subjects, ")"),
           subtitle = paste0("BW: ", input$bw, "kg | ISS: ", input$iss, " | ", input$igg_type,
                             " | Green zone: Optimal range")) +
      theme_bw(base_size = 14) +
      theme(plot.title = element_text(face = "bold"),
            plot.subtitle = element_text(size = 11, color = "gray40")) +
      coord_cartesian(ylim = c(0.001, NA), xlim = c(0, NA))
  })

  #---------------------------------------------------------------------------
  # AI Dosing Optimizer
  #---------------------------------------------------------------------------

  # Patient summary for AI tab
  output$ai_patient_summary <- renderText({
    paste0("BW: ", input$bw, " kg\n",
           "ISS: ", input$iss, "\n",
           "IgG: ", input$igg_type, "\n",
           "Current: ", input$dose1, "/", input$dose2, "/", input$dose_treat, " mg/kg")
  })

  # Run AI optimization
  ai_result <- eventReactive(input$run_ai, {
    # Validate input - use shiny::validate instead of req
    ai_request_text <- input$ai_request
    if (is.null(ai_request_text) || nchar(trimws(ai_request_text)) == 0) {
      return(NULL)
    }

    withProgress(message = 'AI Optimizer running...', value = 0, {

      # Parse request
      goals <- parse_ai_request(ai_request_text)
      incProgress(0.1, detail = "Analyzing request...")

      # Current base parameters
      base_params <- list(
        dose1 = input$dose1,
        dose2 = input$dose2,
        dose_treat = input$dose_treat,
        day_stepup1 = input$day_stepup1,
        day_stepup2 = input$day_stepup2,
        day_treatment = input$day_treatment,
        n_treatment_doses = input$n_treatment_doses
      )

      # Generate scenarios
      scenarios <- generate_dosing_scenarios(base_params, goals)
      incProgress(0.2, detail = "Generated scenarios...")

      params_typical <- typical_params()

      # Run simulations for each scenario
      results <- list()
      n_scenarios <- length(scenarios)

      for (i in seq_along(scenarios)) {
        scenario_name <- names(scenarios)[i]
        params <- scenarios[[scenario_name]]

        incProgress(0.6 / n_scenarios, detail = paste0("Simulating: ", scenario_name))

        sim_result <- tryCatch({
          run_single_simulation(
            mod = mod,
            typical_params = params_typical,
            bw = input$bw,
            dose1 = params$dose1,
            dose2 = params$dose2,
            dose_treat = params$dose_treat,
            day_stepup1 = params$day_stepup1,
            day_stepup2 = params$day_stepup2,
            day_treatment = params$day_treatment,
            n_treatment_doses = params$n_treatment_doses,
            sim_duration = 14,
            n_subjects = input$ai_n_subjects,
            seed = input$seed + i
          )
        }, error = function(e) NULL)

        if (!is.null(sim_result)) {
          results[[scenario_name]] <- list(
            params = params,
            cavg_day5 = sim_result$cavg_day5,
            cavg_day3 = sim_result$cavg_day3,
            score = score_scenario(sim_result$cavg_day5$median, goals),
            risk = classify_risk(sim_result$cavg_day5$median)
          )
        }
      }

      incProgress(0.1, detail = "Complete!")

      list(
        goals = goals,
        results = results,
        request = ai_request_text
      )
    })
  })

  # AI Recommendation UI
  output$ai_recommendation_ui <- renderUI({
    result <- ai_result()
    if (is.null(result) || length(result$results) == 0) return(NULL)

    # Find best scenario
    scores <- sapply(result$results, function(x) x$score)
    if (length(scores) == 0) return(NULL)

    best_name <- names(which.max(scores))
    best_result <- result$results[[best_name]]

    if (is.null(best_result)) return(NULL)

    params <- best_result$params

    div(class = "recommendation-box",
        h4(icon("star"), " AI Recommendation"),
        hr(style = "border-color: rgba(255,255,255,0.3);"),
        p(strong("Best Scenario: "), best_name),
        p(sprintf("Predicted Cavg Day 5: %.4f mg/L (%s)", best_result$cavg_day5$median, best_result$risk)),
        hr(style = "border-color: rgba(255,255,255,0.3);"),
        p(strong("Recommended Dosing:")),
        tags$ul(
          tags$li(sprintf("Step-up 1: %.3f mg/kg (Day %d)", params$dose1, params$day_stepup1)),
          tags$li(sprintf("Step-up 2: %.3f mg/kg (Day %d)", params$dose2, params$day_stepup2)),
          tags$li(sprintf("Treatment: %.2f mg/kg (Day %d+)", params$dose_treat, params$day_treatment))
        )
    )
  })

  # Scenario comparison table
  output$ai_scenario_table <- renderDT({
    result <- ai_result()
    if (is.null(result) || length(result$results) == 0) return(NULL)

    # Build comparison table
    comparison <- map_dfr(names(result$results), function(name) {
      r <- result$results[[name]]
      tibble(
        Scenario = name,
        `Step-up 1` = sprintf("%.3f (D%d)", r$params$dose1, r$params$day_stepup1),
        `Step-up 2` = sprintf("%.3f (D%d)", r$params$dose2, r$params$day_stepup2),
        Treatment = sprintf("%.2f (D%d)", r$params$dose_treat, r$params$day_treatment),
        `Cavg Day 5` = sprintf("%.4f", r$cavg_day5$median),
        `90% PI` = sprintf("[%.4f - %.4f]", r$cavg_day5$q5, r$cavg_day5$q95),
        Risk = r$risk,
        Score = round(r$score, 1)
      )
    }) %>%
      arrange(desc(Score))

    datatable(
      comparison,
      options = list(
        pageLength = 15,
        dom = 't',
        ordering = TRUE
      ),
      rownames = FALSE,
      class = 'cell-border stripe'
    ) %>%
      formatStyle(
        'Risk',
        backgroundColor = styleEqual(
          c("Optimal", "Low Efficacy Risk", "High CRS Risk"),
          c("#d4edda", "#fff3cd", "#f8d7da")
        )
      ) %>%
      formatStyle(
        'Score',
        background = styleColorBar(range(comparison$Score), '#667eea'),
        backgroundSize = '90% 70%',
        backgroundRepeat = 'no-repeat',
        backgroundPosition = 'left'
      )
  })

  # Comparison plot
  output$ai_comparison_plot <- renderPlot({
    result <- ai_result()
    if (is.null(result) || length(result$results) == 0) return(NULL)

    plot_data <- map_dfr(names(result$results), function(name) {
      r <- result$results[[name]]
      tibble(
        Scenario = name,
        Cavg = r$cavg_day5$median,
        Lower = r$cavg_day5$q5,
        Upper = r$cavg_day5$q95,
        Risk = r$risk
      )
    }) %>%
      mutate(Scenario = factor(Scenario, levels = Scenario[order(Cavg)]))

    ggplot(plot_data, aes(x = Scenario, y = Cavg, fill = Risk)) +
      geom_hline(yintercept = OPTIMAL_LOWER, linetype = "dashed", color = "#ffc107", linewidth = 1) +
      geom_hline(yintercept = OPTIMAL_UPPER, linetype = "dashed", color = "#dc3545", linewidth = 1) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = OPTIMAL_LOWER, ymax = OPTIMAL_UPPER,
               fill = "#28a745", alpha = 0.15) +
      geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.3, color = "gray50") +
      geom_point(size = 4, shape = 21, color = "black") +
      scale_fill_manual(values = c(
        "Optimal" = "#28a745",
        "Low Efficacy Risk" = "#ffc107",
        "High CRS Risk" = "#dc3545"
      )) +
      labs(
        x = "Dosing Scenario",
        y = "Cavg Day 5 (mg/L)",
        title = "Cavg Day 5 Comparison Across Scenarios",
        subtitle = "Error bars: 90% prediction interval | Green zone: Optimal range"
      ) +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"),
        legend.position = "bottom"
      )
  })
}

#-------------------------------------------------------------------------------
# Run Application
#-------------------------------------------------------------------------------

shinyApp(ui = ui, server = server)
