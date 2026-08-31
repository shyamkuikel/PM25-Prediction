# ==========================================================
# DLNM ANALYSIS: TEMPERATURE AND PM2.5
# ==========================================================
# Purpose:
# Evaluate nonlinear and delayed associations between
# temperature and PM2.5 concentrations across monitoring
# stations using Distributed Lag Nonlinear Models (DLNM).
#
# Covariates:
# Relative humidity, wind speed, atmospheric pressure,
# day of week, and temporal trend.
#
# Author: Shyam Prasad Kuikel
# ==========================================================


# ==========================================================
# 1. INSTALL AND LOAD REQUIRED PACKAGES
# ==========================================================

required_packages <- c("readxl", "dlnm", "splines")

installed_packages <- rownames(installed.packages())

for (pkg in required_packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg)
  }
}

library(readxl)
library(dlnm)
library(splines)


# ==========================================================
# 2. INPUT DATA
# ==========================================================

file_name <- "Meterological data ALAbama.xlsx"

all_sheets <- excel_sheets(file_name)

cat("Available stations/sheets:\n")
print(all_sheets)


# ==========================================================
# 3. FUNCTION FOR STATION-SPECIFIC DLNM ANALYSIS
# ==========================================================

get_station_result <- function(sheet_name,
                               file_name,
                               lag_days = 15) {

  cat("\nProcessing:", sheet_name, "\n")

  # Read station data
  data <- tryCatch(
    read_excel(file_name, sheet = sheet_name),
    error = function(e) NULL
  )

  if (is.null(data)) {
    cat("  -> Could not read sheet\n")
    return(NULL)
  }


  # --------------------------------------------------------
  # Required variables
  # --------------------------------------------------------

  needed_cols <- c(
    "Only_Date",
    "Temperature",
    "Relative Humidity",
    "Wind Speed",
    "Pressure",
    "Daily Mean PM2.5 Concentration"
  )

  if (!all(needed_cols %in% names(data))) {
    cat("  -> Missing required columns\n")
    return(NULL)
  }


  # Keep required columns only
  data <- data[, needed_cols]

  names(data) <- c(
    "date",
    "temp",
    "rh",
    "ws",
    "pressure",
    "pm25"
  )


  # --------------------------------------------------------
  # Data type conversion
  # --------------------------------------------------------

  data$date <- as.Date(data$date)
  data$temp <- as.numeric(data$temp)
  data$rh <- as.numeric(data$rh)
  data$ws <- as.numeric(data$ws)
  data$pressure <- as.numeric(data$pressure)
  data$pm25 <- as.numeric(data$pm25)


  # Remove incomplete observations
  data <- na.omit(data)


  # --------------------------------------------------------
  # Basic quality checks
  # --------------------------------------------------------

  if (nrow(data) < 30) {
    cat("  -> Too few observations\n")
    return(NULL)
  }

  if (sd(data$temp, na.rm = TRUE) < 1) {
    cat("  -> Low temperature variation\n")
    return(NULL)
  }

  if (sd(data$pm25, na.rm = TRUE) < 0.2) {
    cat("  -> Low PM2.5 variation\n")
    return(NULL)
  }


  # ========================================================
  # 4. TIME VARIABLES
  # ========================================================

  data$dow <- factor(
    weekdays(data$date),
    levels = c(
      "Monday",
      "Tuesday",
      "Wednesday",
      "Thursday",
      "Friday",
      "Saturday",
      "Sunday"
    )
  )

  data$time <- as.numeric(data$date)


  # ========================================================
  # 5. TEMPERATURE CROSS-BASIS
  # ========================================================

  cb.temp <- tryCatch(
    crossbasis(
      data$temp,
      lag = lag_days,
      argvar = list(
        fun = "ns",
        df = 4
      ),
      arglag = list(
        fun = "ns",
        df = 4
      )
    ),
    error = function(e) NULL
  )

  if (is.null(cb.temp)) {
    cat("  -> Cross-basis construction failed\n")
    return(NULL)
  }


  # ========================================================
  # 6. FIT DLNM
  # ========================================================

  model <- tryCatch(
    glm(
      pm25 ~ cb.temp +
        ns(rh, df = 2) +
        ns(ws, df = 2) +
        ns(pressure, df = 2) +
        dow +
        ns(time, df = 4),
      family = gaussian(),
      data = data
    ),
    error = function(e) NULL
  )

  if (is.null(model)) {
    cat("  -> Model fitting failed\n")
    return(NULL)
  }


  # ========================================================
  # 7. DLNM PREDICTION
  # ========================================================

  temp_seq <- seq(
    floor(min(data$temp, na.rm = TRUE)),
    ceiling(max(data$temp, na.rm = TRUE)),
    by = 1
  )

  reference_temp <- median(
    data$temp,
    na.rm = TRUE
  )

  pred <- tryCatch(
    crosspred(
      cb.temp,
      model,
      at = temp_seq,
      cen = reference_temp,
      bylag = 1
    ),
    error = function(e) NULL
  )

  if (is.null(pred) || is.null(pred$matfit)) {
    cat("  -> Prediction failed\n")
    return(NULL)
  }

  z <- pred$matfit

  if (all(is.na(z))) {
    cat("  -> Prediction surface contains only NA values\n")
    return(NULL)
  }

  lag_seq <- 0:lag_days


  # ========================================================
  # 8. CHECK PREDICTION MATRIX ORIENTATION
  # ========================================================

  if (
    nrow(z) == length(temp_seq) &&
    ncol(z) == length(lag_seq)
  ) {

    # Matrix orientation is already correct

  } else if (
    nrow(z) == length(lag_seq) &&
    ncol(z) == length(temp_seq)
  ) {

    z <- t(z)

  } else {

    cat(
      "  -> Prediction matrix dimension mismatch:",
      "nrow =", nrow(z),
      "ncol =", ncol(z),
      "\n"
    )

    return(NULL)
  }


  # Replace isolated missing prediction values
  if (anyNA(z)) {
    z[is.na(z)] <- mean(z, na.rm = TRUE)
  }


  # Skip nearly flat surfaces
  if (sd(as.vector(z), na.rm = TRUE) < 0.05) {
    cat("  -> Very low variation in prediction surface\n")
    return(NULL)
  }


  local_min <- min(z, na.rm = TRUE)
  local_max <- max(z, na.rm = TRUE)

  if (
    !is.finite(local_min) ||
    !is.finite(local_max) ||
    local_min == local_max
  ) {

    cat("  -> Invalid prediction range\n")
    return(NULL)
  }


  # Return station results
  list(
    name = sheet_name,
    temp = temp_seq,
    lag = lag_seq,
    z = z,
    reference_temp = reference_temp,
    local_min = local_min,
    local_max = local_max
  )
}


# ==========================================================
# 9. RUN ANALYSIS FOR ALL STATIONS
# ==========================================================

results <- lapply(
  all_sheets,
  get_station_result,
  file_name = file_name,
  lag_days = 15
)

names(results) <- all_sheets

results <- results[
  !sapply(results, is.null)
]

if (length(results) == 0) {
  stop("No valid stations found.")
}

cat("\nValid stations:\n")
print(names(results))


# ==========================================================
# 10. COLOR PALETTE
# ==========================================================

cols <- colorRampPalette(
  c(
    "darkblue",
    "blue",
    "white",
    "red",
    "darkred"
  )
)(100)


# ==========================================================
# 11. GROUP STATIONS
# ==========================================================

plot_names <- names(results)

groups <- split(
  plot_names,
  ceiling(seq_along(plot_names) / 4)
)


# ==========================================================
# 12. GENERATE DLNM CONTOUR FIGURES
# ==========================================================

for (g in seq_along(groups)) {

  current_names <- groups[[g]]

  out_name <- paste0(
    "Temperature_SeparateCbar_Figure_",
    g,
    ".png"
  )


  png(
    out_name,
    width = 4200,
    height = 2600,
    res = 300
  )


  layout(
    matrix(
      c(
        1, 2, 3, 4,
        5, 6, 7, 8
      ),
      nrow = 2,
      byrow = TRUE
    ),
    widths = c(
      1, 0.16,
      1, 0.16
    ),
    heights = c(1, 1)
  )


  for (i in 1:4) {

    if (i <= length(current_names)) {

      res <- results[[current_names[i]]]

      z <- res$z

      breaks_local <- seq(
        res$local_min,
        res$local_max,
        length.out = 101
      )


      # ----------------------------------------------------
      # MAIN DLNM PLOT
      # ----------------------------------------------------

      par(
        mar = c(5, 5, 3, 1.5)
      )

      if (!(
        nrow(z) == length(res$temp) &&
        ncol(z) == length(res$lag)
      )) {

        plot.new()

        title(
          main = paste(
            res$name,
            "\nInvalid matrix shape"
          )
        )

      } else {

        image(
          x = res$temp,
          y = res$lag,
          z = z,
          col = cols,
          breaks = breaks_local,
          xlab = "Temperature (°C)",
          ylab = "Lag (Days)",
          main = res$name,
          cex.main = 1.2,
          cex.lab = 1.1,
          cex.axis = 1
        )


        contour(
          x = res$temp,
          y = res$lag,
          z = z,
          add = TRUE,
          drawlabels = FALSE,
          col = "black",
          lwd = 1
        )

        box()
      }


      # ----------------------------------------------------
      # STATION-SPECIFIC COLOR BAR
      # ----------------------------------------------------

      par(
        mar = c(5, 1, 3, 4)
      )

      y_vals <- seq(
        res$local_min,
        res$local_max,
        length.out = 100
      )

      x_vals <- c(0, 1)

      z_bar <- matrix(
        rep(
          y_vals,
          each = 2
        ),
        nrow = 2,
        byrow = TRUE
      )


      image(
        x = x_vals,
        y = y_vals,
        z = z_bar,
        col = cols,
        breaks = breaks_local,
        xaxt = "n",
        yaxt = "n",
        xlab = "",
        ylab = ""
      )

      axis(
        4,
        las = 1,
        cex.axis = 0.9
      )

      mtext(
        "Effect",
        side = 4,
        line = 2,
        cex = 1
      )

      box()

    } else {

      par(
        mar = c(5, 5, 3, 1.5)
      )

      plot.new()

      par(
        mar = c(5, 1, 3, 4)
      )

      plot.new()
    }
  }


  mtext(
    "DLNM Contour: Temperature and PM2.5",
    side = 3,
    outer = TRUE,
    line = -1.5,
    cex = 1.5,
    font = 2
  )


  dev.off()

  cat(
    "Saved:",
    out_name,
    "\n"
  )
}
