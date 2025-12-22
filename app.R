#===============================================================================
# Teclistamab PK Simulation Shiny App
#
# Individual PK parameter calculation based on patient characteristics
# 2-Compartment model with time-dependent clearance
#
# Reference: Miao et al. (2023) - Teclistamab Population PK Model
#===============================================================================

library(shiny)
library(mrgsolve)
library(tidyverse)
library(PKNCA)
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
# PK Parameter Calculation Functions
#-------------------------------------------------------------------------------

calculate_pk_params <- function(bw, iss, igg_type) {
  # ISS indicator variables
  iss_2 <- ifelse(iss == "II", 1, 0)
  iss_3 <- ifelse(iss == "III", 1, 0)

  # IgG indicator (Non-IgG = 1)
  non_igg <- ifelse(igg_type == "Non-IgG", 1, 0)

  # Calculate CL1 (L/day)
  # CL1 = 0.449 × (BWT/74)^0.704 × 1.31^(ISS=II) × 1.67^(ISS=III) × 0.689^(Non-IgG)
  CL1 <- 0.449 * (bw / 74)^0.704 * (1.31^iss_2) * (1.67^iss_3) * (0.689^non_igg)

  # Calculate CL2 (L/day)
  # CL2 = 0.547 × 0.295^(Non-IgG)
  CL2 <- 0.547 * (0.295^non_igg)

  # Calculate V1 (L)
  # V1 = 4.13 × (BWT/74)^0.358
  V1 <- 4.13 * (bw / 74)^0.358

  # Calculate V2 (L)
  # V2 = 1.34 × (BWT/74)^1.40
  V2 <- 1.34 * (bw / 74)^1.40

  # Fixed parameters
  Q <- 0.039      # Intercompartmental clearance (L/day)
  KA <- 0.133     # Absorption rate constant (1/day)
  F1 <- 0.718     # Bioavailability
  KDES <- 0.0292  # Clearance decay rate (1/day)

  list(
    CL1 = CL1,
    CL2 = CL2,
    V1 = V1,
    V2 = V2,
    Q = Q,
    KA = KA,
    F1 = F1,
    KDES = KDES
  )
}

#-------------------------------------------------------------------------------
# NCA Calculation Function
#-------------------------------------------------------------------------------

calculate_nca <- function(sim_data) {
  # Use data after last dose for steady-state NCA
  conc_data <- sim_data %>%
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

  # Calculate terminal half-life (using last 20% of data points)
  n_points <- nrow(conc_data)
  terminal_start <- max(1, floor(n_points * 0.8))
  terminal_data <- conc_data[terminal_start:n_points, ]

  if (nrow(terminal_data) >= 3 && all(terminal_data$DV > 0)) {
    log_conc <- log(terminal_data$DV)
    time_vals <- terminal_data$TIME_DAY

    fit <- lm(log_conc ~ time_vals)
    lambda_z <- -coef(fit)[2]

    if (lambda_z > 0) {
      t_half <- log(2) / lambda_z
      # AUC extrapolated to infinity
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

  # Cavg (average concentration)
  Cavg <- AUC_total / (Tlast - conc_data$TIME_DAY[1])

  # Cmin (trough concentration - minimum after Tmax)
  post_tmax <- conc_data %>% filter(TIME_DAY > Tmax)
  Cmin <- if(nrow(post_tmax) > 0) min(post_tmax$DV) else NA

  tibble(
    Parameter = c("Cmax", "Tmax", "Cmin", "Cavg", "AUC(0-last)", "AUC(0-inf)", "t1/2", "Lambda_z"),
    Value = c(Cmax, Tmax, Cmin, Cavg, AUC_total, AUC_inf, t_half, lambda_z),
    Unit = c("mg/L", "day", "mg/L", "mg/L", "mg·day/L", "mg·day/L", "day", "1/day")
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
    "))
  ),

  # Header
  div(class = "main-header",
      h1("Teclistamab PK Simulator", style = "margin: 0;"),
      p("Individual PK Simulation based on Patient Characteristics", style = "margin: 5px 0 0 0; opacity: 0.9;")
  ),

  # Main Layout
  fluidRow(
    # Input Panel
    column(4,
           div(class = "param-box",
               h4(icon("user"), " Patient Characteristics"),
               hr(),
               numericInput("bw",
                            "Body Weight (kg):",
                            value = 70,
                            min = 30,
                            max = 150,
                            step = 0.1),
               selectInput("iss",
                           "ISS Stage:",
                           choices = c("I", "II", "III"),
                           selected = "II"),
               selectInput("igg_type",
                           "Immunoglobulin Type:",
                           choices = c("IgG", "Non-IgG"),
                           selected = "IgG")
           ),

           div(class = "param-box",
               h4(icon("syringe"), " Dosing Information"),
               hr(),
               numericInput("dose",
                            "Dose (mg/kg):",
                            value = 1.5,
                            min = 0.1,
                            max = 10,
                            step = 0.1),
               numericInput("n_doses",
                            "Number of Doses:",
                            value = 8,
                            min = 1,
                            max = 50,
                            step = 1),
               selectInput("dosing_interval",
                           "Dosing Interval:",
                           choices = c("Weekly (QW)" = 7,
                                       "Every 2 Weeks (Q2W)" = 14),
                           selected = 7),
               numericInput("sim_duration",
                            "Simulation Duration after Last Dose (days):",
                            value = 21,
                            min = 7,
                            max = 90,
                            step = 1)
           ),

           div(style = "text-align: center; margin-top: 20px;",
               actionButton("simulate", "Run Simulation",
                            class = "btn btn-primary btn-simulate btn-lg",
                            icon = icon("play"))
           )
    ),

    # Output Panel
    column(8,
           # Calculated PK Parameters
           div(class = "result-box", style = "margin-bottom: 20px;",
               h4(icon("calculator"), " Calculated PK Parameters"),
               hr(),
               tableOutput("pk_params_table")
           ),

           # PK Curve
           div(class = "result-box", style = "margin-bottom: 20px;",
               h4(icon("chart-line"), " Time-Concentration Profile"),
               hr(),
               plotOutput("pk_plot", height = "400px")
           ),

           # NCA Parameters
           div(class = "result-box",
               h4(icon("table"), " NCA Parameters"),
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

  # Reactive: Calculate PK parameters
  pk_params <- reactive({
    calculate_pk_params(input$bw, input$iss, input$igg_type)
  })

  # Display calculated PK parameters
  output$pk_params_table <- renderTable({
    params <- pk_params()

    data.frame(
      Parameter = c("CL1 (Linear CL)", "CL2 (Time-dep CL)", "V1 (Central)",
                    "V2 (Peripheral)", "Q", "KA", "F1", "KDES"),
      Value = c(
        sprintf("%.4f", params$CL1),
        sprintf("%.4f", params$CL2),
        sprintf("%.3f", params$V1),
        sprintf("%.4f", params$V2),
        sprintf("%.4f", params$Q),
        sprintf("%.3f", params$KA),
        sprintf("%.3f", params$F1),
        sprintf("%.4f", params$KDES)
      ),
      Unit = c("L/day", "L/day", "L", "L", "L/day", "1/day", "-", "1/day")
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE, width = "100%")

  # Reactive: Run simulation
  sim_result <- eventReactive(input$simulate, {

    withProgress(message = 'Running simulation...', value = 0, {

      params <- pk_params()

      # Calculate actual dose in mg
      dose_mg <- input$dose * input$bw

      # Create dosing schedule (TIME in hours)
      interval_hours <- as.numeric(input$dosing_interval) * 24
      dose_times <- seq(0, (input$n_doses - 1) * interval_hours, by = interval_hours)

      dosing_data <- tibble(
        ID = 1,
        time = dose_times,
        amt = dose_mg,
        cmt = 1,
        evid = 1
      )

      incProgress(0.3, detail = "Creating dosing schedule...")

      # Simulation time grid (hourly)
      max_time <- max(dose_times) + input$sim_duration * 24
      sim_times <- seq(0, max_time, by = 1)

      obs_data <- tibble(
        ID = 1,
        time = sim_times,
        amt = 0,
        cmt = 0,
        evid = 0
      )

      sim_data <- bind_rows(dosing_data, obs_data) %>%
        arrange(time, desc(evid))

      incProgress(0.3, detail = "Running mrgsolve...")

      # Update model with individual parameters
      mod_i <- mod %>%
        param(
          CL1 = params$CL1,
          CL2 = params$CL2,
          KDES = params$KDES,
          V1 = params$V1,
          V2 = params$V2,
          Q = params$Q,
          KA = params$KA,
          F1 = params$F1
        )

      # Run simulation
      out <- mod_i %>%
        data_set(sim_data) %>%
        mrgsim(carry_out = "amt,evid") %>%
        as_tibble() %>%
        filter(evid == 0) %>%
        mutate(
          TIME_HOUR = time,
          TIME_DAY = time / 24
        )

      incProgress(0.4, detail = "Complete!")

      list(
        simulation = out,
        dosing = dosing_data,
        params = params
      )
    })
  })

  # Plot PK curve
  output$pk_plot <- renderPlot({
    req(sim_result())

    result <- sim_result()
    sim_data <- result$simulation
    dosing <- result$dosing

    # Dose times for vertical lines
    dose_days <- dosing$time / 24

    ggplot(sim_data, aes(x = TIME_DAY, y = DV)) +
      geom_line(color = "#667eea", linewidth = 1.2) +
      geom_vline(xintercept = dose_days, linetype = "dashed",
                 color = "#e74c3c", alpha = 0.5) +
      scale_y_log10(
        breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
        labels = c("0.001", "0.01", "0.1", "1", "10", "100")
      ) +
      annotation_logticks(sides = "l") +
      labs(
        x = "Time (days)",
        y = "Concentration (mg/L)",
        title = paste0("Teclistamab PK Profile (", input$dose, " mg/kg, ",
                       input$iss, " stage, ", input$igg_type, ")"),
        subtitle = paste0("BW: ", input$bw, " kg | Dose: ",
                          round(input$dose * input$bw, 1), " mg | ",
                          input$n_doses, " doses")
      ) +
      theme_bw(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_line(color = "gray90"),
        legend.position = "bottom"
      )
  })

  # Calculate and display NCA parameters
  output$nca_table <- renderDT({
    req(sim_result())

    result <- sim_result()
    nca_params <- calculate_nca(result$simulation)

    if (is.null(nca_params)) {
      return(NULL)
    }

    nca_params %>%
      mutate(Value = ifelse(is.na(Value), "N/A", sprintf("%.4f", Value))) %>%
      datatable(
        options = list(
          dom = 't',
          pageLength = 10,
          ordering = FALSE
        ),
        rownames = FALSE,
        class = 'cell-border stripe'
      )
  })
}

#-------------------------------------------------------------------------------
# Run Application
#-------------------------------------------------------------------------------

shinyApp(ui = ui, server = server)
