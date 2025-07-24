####-------------------------- Beginning of code for application to spectra data --------------------------####
# ============================================================================
# File: Application_spectra.R
# Purpose: Load FNMR data, estimate derivatives, and generate plots for analysis
#
# This script demonstrates how to use the FNMR data processing and derivative estimation
# functions defined in FNMR.R. It loads the data, preprocesses it, estimates derivatives
# for all time points, and visualizes the results for inspection.
# ============================================================================
rm(list = ls())  # Clear the workspace
# Load required libraries
library(here)      # For robust file path handling
library(ggplot2)   # For advanced plotting
library(dplyr)     # For data manipulation
library(tidyr)     # For data tidying

# Source the FNMR processing functions
source(here::here("R", "FNMR.R"))

# ------------------------------------------------------------------------------
# 1. Load and preprocess the FNMR data
# ------------------------------------------------------------------------------
# Load the FNMR data for all time points (returns a named list of data frames)
fnmr_data <- load_data()
# Normalize the chemical shift for each time point to [0, 1]
fnmr_data <- lapply(fnmr_data, preprocess_shift)

head(fnmr_data$data.0)  # Preview the first few rows of the first time point
# ------------------------------------------------------------------------------
# 2. Estimate derivatives for all time points
# ------------------------------------------------------------------------------
# This will add fitted values and both corrected/uncorrected derivatives to each data frame
# processed_results <- process_all_fnmr(
#     fnmr_data,
#     nseg = 6,  # Number of segments for smoothing
#     bdeg = 4,  # B-spline degree
#     pord = 2,  # Polynomial order for smoothing
#     x.grid = NULL
# )

processed_results <- process_all_fnmr(
    fnmr_data,
    nseg = 6,  # Number of segments for smoothing
    # c(5, 6, 5, 5, 5, 5, 5, 5),  # or named vector
    bdeg = 4,
    pord = 2
)
# ------------------------------------------------------------------------------
# 3. Inspect the processed data for a specific time point (e.g., 0 hours)
# ------------------------------------------------------------------------------
cat("\nPreview of processed data for 0 hours:\n")
print(head(processed_results$data.0))
cat("\nRange of chem_orig for 0 hours: ", range(processed_results$data.0$chem_orig), "\n")

processed_results1 <- processed_results  # Keep a copy for later use
# ------------------------------------------------------------------------------
# Plot the spectrum for each time point and mark the peak (maximum intensity)
# ------------------------------------------------------------------------------
# This will show the spectrum and mark the peak for each time point, allowing manual verification
plot_fnmr_peaks(processed_results1, xvar = 'chem_orig')


# ------------------------------------------------------------------------------
# 8x3 grid panel plot: mean, 1st, and 2nd derivatives for all 8 time points (base R)
# ------------------------------------------------------------------------------
# This will show the fitted mean, first and second derivatives (naive & resub) for all 8 time points
plot_fnmr_panels_base(processed_results1, xvar = 'chem_orig')


# ------------------------------------------------------------------------------
# Zoom in on the derivatives around the true peak (from last time point)
# ------------------------------------------------------------------------------
plot_fnmr_deriv_zoom(processed_results1, xvar = 'chem_orig', window = 1)

# ------------------------------------------------------------------------------
# Table comparing the first and second derivative values at the true peak for each time point
# ------------------------------------------------------------------------------
deriv_at_peak_table <- get_deriv_at_true_peak(processed_results1, xvar = 'chem_orig')
cat('\nDerivative values at the true peak (from last time point):\n')
print(deriv_at_peak_table)



## ploting maxima detection with  no interpolation####
plot_fnmr_panels_maxima_detection_interp(processed_results1, xvar = 'chem_orig', boundary_n = 500, save_pdf = TRUE, true_max_shift = -37.43,
 use_interpolation = FALSE)


## ploting maxima detection with interpolation####
plot_fnmr_panels_maxima_detection_interp(processed_results1, xvar = 'chem_orig', boundary_n = 500, save_pdf = TRUE, true_max_shift = -38,
 use_interpolation = TRUE)
