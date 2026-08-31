# PM₂.₅ Prediction and Environmental Modeling

## Overview

This repository contains research code for analyzing PM₂.₅ concentrations and their relationships with meteorological variables using statistical modeling, machine learning, and deep learning approaches.

The current implementation includes a Distributed Lag Nonlinear Model (DLNM) workflow for investigating nonlinear and delayed associations between temperature and PM₂.₅.

## Current Analysis

### Distributed Lag Nonlinear Model (DLNM)

The R script `R/dlnm_temperature_pm25.R` performs station-specific DLNM analysis of temperature and PM₂.₅.

The workflow includes:

- Multi-station data processing
- Data quality checks
- Temperature cross-basis construction
- Nonlinear exposure-response modeling
- Lagged-effect analysis
- Adjustment for relative humidity, wind speed, and atmospheric pressure
- Adjustment for day-of-week and temporal trends
- Station-specific prediction
- DLNM contour visualization

The current analysis considers lag periods from 0 to 15 days.

## Required Variables

The input dataset should contain the following variables:

- Date
- Temperature
- Relative Humidity
- Wind Speed
- Atmospheric Pressure
- Daily Mean PM₂.₅ Concentration

## Statistical Framework

Temperature is modeled using a DLNM cross-basis with natural cubic splines for both the exposure and lag dimensions.

PM₂.₅ concentration is treated as a continuous response using a Gaussian generalized linear model. Relative humidity, wind speed, atmospheric pressure, day of week, and temporal trend are included as covariates.

## Repository Structure

```text
PM25-Prediction/
├── R/
│   └── dlnm_temperature_pm25.R
└── README.md
