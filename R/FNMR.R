
## ------------------------------------------------------------------
# Example usage: Process all FNMR data and inspect results
# Uncomment and run the following lines in your R session:
#
# # Load and preprocess the data
# fnmr_data <- load_data()
# fnmr_data <- lapply(fnmr_data, preprocess_shift)
#
# # Estimate derivatives for all time points
# processed_results <- process_all_fnmr(fnmr_data)
#
# # Inspect the processed data for one time point (e.g., 0 hours)
# head(processed_results$data.0)
#
# # Plot the fitted values and derivatives if desired
# plot(processed_results$data.0$chem, processed_results$data.0$fitted, type='l', main='Fitted (Naive) - 0 Hours')
# plot(processed_results$data.0$chem, processed_results$data.0$naive_deriv, type='l', main='Naive Derivative - 0 Hours')
# plot(processed_results$data.0$chem, processed_results$data.0$resub_deriv, type='l', main='Resub Derivative - 0 Hours')

#===================================================================#
# FNMR Data Processing and Visualization Script
# ---------------------------------------------------------------
# This script loads, processes, and visualizes FNMR (Fluorescence Nuclear Magnetic Resonance)
# data collected at different time points (hours). It generates plots for each time point,
# preprocesses the data for further analysis, and provides utility functions for derivative correction.
#
# Author: Christopher Odoom
# Date: July 14, 2025
#===================================================================#


# Load the 'here' package for robust file path handling
library(here)

# Source the necessary functions for data loading and plotting
source(here::here("R", "MainFunctions.R"))

#' Load FNMR Data for Specified Hours
#'
#' Reads the FNMR CSV file and extracts chemical shift and intensity columns
#' for each specified hour. Returns a named list of data frames for each hour.
#'
#' @param hour Numeric vector of hours to extract data for (default: all available)
#' @return Named list of data frames, each with columns 'chem' and 'Inten'
load_data <- function(hour = c(0, 6, 9, 15, 19, 22, 24, 48)) {
  # Load the FNMR data, skipping the first row (header info)
  data <- read.csv(here::here("data", "raw", "FNMRData.csv"), skip=1)
  data_list <- list()
  for (h in hour) {
    # Determine the column indices for chemical shift and intensity for each hour
    idx <- which(hour == h)
    col_chem <- 2 * idx - 1
    col_inten <- 2 * idx
    dat <- data[, c(col_chem, col_inten)]
    # Remove rows with missing chemical shift values
    dat <- dat[!is.na(dat[, 1]), ]
    colnames(dat) <- c("chem", "Inten")
    # Store in the list with a name indicating the hour
    data_list[[paste0("data.", h)]] <- dat
  }
   return(data_list)
}

 

# data.0 <- fnmr_data$data.0
# data.6 <- fnmr_data$data.6
# data.9 <- fnmr_data$data.9
# data.15 <- fnmr_data$data.15
# data.19 <- fnmr_data$data.19
# data.22 <- fnmr_data$data.22
# data.24 <- fnmr_data$data.24
# data.48 <- fnmr_data$data.48

# Load all FNMR data for the specified hours
#fnmr_data <- load_data()

# # Extract data for each hour into separate variables for convenience
# data.0  <- fnmr_data$data.0
# data.6  <- fnmr_data$data.6
# data.9  <- fnmr_data$data.9
# data.15 <- fnmr_data$data.15
# data.19 <- fnmr_data$data.19
# data.22 <- fnmr_data$data.22
# data.24 <- fnmr_data$data.24
# data.48 <- fnmr_data$data.48


#' Plot FNMR Data for a Single Time Point
#'
#' Plots chemical shift vs intensity for a given data frame.
#' @param data Data frame with columns 'chem' and 'Inten'
#' @param title Plot title (typically the hour)
plot_data <- function(data, title) {
  plot(data$chem, data$Inten, main=title,
       type="l", axes=FALSE, xlab="Shift", ylab="Intensity")
  points(data$chem, data$Inten, pch=16, cex=0.5)
  axis(1)
  axis(2)
}

# # Generate and save plots for all time points in a single PNG file
# png(filename = here::here("output", "plots", "FNMR_plots", "FNMR_plots.png"), width=1600, height=800, res=150)
# par(mfrow=c(2,4)) # Arrange plots in a 2x4 grid
# plot_data(data.0,  "0 Hours")
# plot_data(data.6,  "6 Hours")
# plot_data(data.9,  "9 Hours")
# plot_data(data.15, "15 Hours")
# plot_data(data.19, "19 Hours")
# plot_data(data.22, "22 Hours")
# plot_data(data.24, "24 Hours")
# plot_data(data.48, "48 Hours")
# dev.off()



#' Preprocess Chemical Shift to [0, 1] Range
#'
#' Normalizes the 'chem' (chemical shift) column of a data frame to the [0, 1] interval.
#' Useful for standardizing the x-axis before fitting models or comparing across samples.
#' @param data Data frame with a 'chem' column
#' @return Data frame with normalized 'chem' column
preprocess_shift <- function(data) {
  data$chem_orig <- data$chem # Store original chemical shift
  data$chem <- (data$chem - min(data$chem)) / (max(data$chem) - min(data$chem))
  return(data)
}

# # Preprocess the chemical shift for all time points
# data.0  <- preprocess_shift(data.0)
# head(data.0) # Show first few rows for inspection
# data.6  <- preprocess_shift(data.6)
# data.9  <- preprocess_shift(data.9)
# data.15 <- preprocess_shift(data.15)
# data.19 <- preprocess_shift(data.19)
# data.22 <- preprocess_shift(data.22)
# data.24 <- preprocess_shift(data.24)
# data.48 <- preprocess_shift(data.48)


#' Correct Derivative Scale After Normalization
#'
#' Adjusts the first or second derivative values back to the original chemical shift scale
#' after fitting on normalized [0, 1] data. This is necessary because differentiation
#' with respect to a scaled variable changes the magnitude of the derivative.
#'
#' @param derivative Numeric vector of derivative values (from normalized scale)
#' @param data Data frame with original 'chem' column (before normalization)
#' @param order Integer, 1 for first derivative, 2 for second derivative
#' @return Numeric vector of corrected derivative values
get_corrected_derivative <- function(derivative, data, order) {
  rng <- max(data$chem) - min(data$chem)
  if (order == 1) {
    # First derivative: divide by the range of the original x-axis
    return(derivative / rng)
  } else if (order == 2) {
    # Second derivative: divide by the square of the range
    return(derivative / (rng^2))
  } else {
    stop("Order must be either 1 or 2")
  }
}



#' Estimate Derivatives Using Naive and Resubstitution Methods
#'
#' This function estimates the mean function and its derivatives using both the naive and
#' resubstitution (resub) methods, as defined in MainFunctions.R. It returns a list containing
#' the results from both methods for further analysis or comparison.
#' @param data Data frame with 'chem' and 'Inten' columns (normalized or original scale)
#' @param r Integer, order of derivative to estimate (default: 1)
#' @param nseg Number of B-spline segments (default: 35)
#' @param bdeg B-spline degree (default: 4)
#' @param pord Penalty order (default: 2)
#' @param x.grid Grid of x values for prediction (default: data$chem)
#' @param compute_CI Logical, if TRUE computes pointwise confidence intervals for derivatives
#' @return List with elements 'naive' and 'resub', each containing the estimated values
estimate_derivative <- function(data, r = 1, nseg = 35, bdeg = 4, pord = 2, x.grid = NULL, compute_CI = TRUE) {
  # Use the same grid as the data if not specified
  if (is.null(x.grid)) x.grid <- data$chem
  # Call the correct functions from MainFunctions.R
  if(r==1){
    pord <- pord-1
    naive_result <- naive.est.opt(x = data$chem, y = data$Inten, r = r, nseg = nseg, bdeg = bdeg, pord = pord, x.grid = x.grid)
    resub_result <- resub.est(x = data$chem, y = data$Inten, r = r, nseg = nseg, bdeg = bdeg, pord = pord, x.grid = x.grid)
  } else if (r == 2) {
    naive_result <- naive.est.opt(x = data$chem, y = data$Inten, r = r, nseg = nseg, bdeg = (bdeg+1), pord = pord, x.grid = x.grid)
    resub_result <- resub.est(x = data$chem, y = data$Inten, r = r, nseg = nseg, bdeg = (bdeg+1), pord = pord, x.grid = x.grid)
  } else {
    stop("Order of derivative must be either 1 or 2")
  }
  # naive_result <- naive.est.opt(x = data$chem, y = data$Inten, r = r, nseg = nseg, bdeg = bdeg, pord = pord, x.grid = x.grid)
  # resub_result <- resub.est(x = data$chem, y = data$Inten, r = r, nseg = nseg, bdeg = bdeg, pord = pord, x.grid = x.grid)

  # Use original chemical shift for correction if available
  chem_for_correction <- if ("chem_orig" %in% names(data)) {
    data.frame(chem = data$chem_orig)
  } else {
    data
  }
  # Correct the derivatives to the original scale
  naive_deriv_corrected <- get_corrected_derivative(naive_result$fr.hat, chem_for_correction, r)
  resub_deriv_corrected <- get_corrected_derivative(resub_result$fr.hat, chem_for_correction, r)

  # Compute CIs if requested
  if (compute_CI) {
    # Naive
    print("The naive lambda is:")
    print(naive_result$lambda)
    print("********************")
    naive_CI <- pointwiseConfInt(
      deriv.est = naive_result$fr.hat,
      r = r,
      BS = naive_result$BS,
      pord = pord,
      sig.est = naive_result$sig.hat,
      lambda.est = naive_result$lambda,
      x = x.grid
    )
    # Resub
    print("The resub lambda is:")
    print(resub_result$lambda)
    print("********************")
    resub_CI <- pointwiseConfInt(
      deriv.est = resub_result$fr.hat,
      r = r,
      BS = naive_result$BS,
      pord = pord,
      sig.est = naive_result$sig.hat,
      lambda.est = resub_result$lambda,
      x = x.grid
    )
    # Correct CIs to original scale
    naive_CI_lower <- get_corrected_derivative(naive_CI$CI.lower, chem_for_correction, r)
    naive_CI_upper <- get_corrected_derivative(naive_CI$CI.upper, chem_for_correction, r)
    resub_CI_lower <- get_corrected_derivative(resub_CI$CI.lower, chem_for_correction, r)
    resub_CI_upper <- get_corrected_derivative(resub_CI$CI.upper, chem_for_correction, r)
  } else {
    naive_CI_lower <- naive_CI_upper <- resub_CI_lower <- resub_CI_upper <- rep(NA, length(x.grid))
  }

  # Add both corrected and uncorrected columns to the data frame
  data$fitted <- as.numeric(naive_result$f.hat)
  data$naive_deriv_uncorrected <- as.numeric(naive_result$fr.hat)
  data$resub_deriv_uncorrected <- as.numeric(resub_result$fr.hat)
  data$naive_deriv <- as.numeric(naive_deriv_corrected)
  data$resub_deriv <- as.numeric(resub_deriv_corrected)
  # Add CIs
  data$naive_deriv_CI.lower <- as.numeric(naive_CI_lower)
  data$naive_deriv_CI.upper <- as.numeric(naive_CI_upper)
  data$resub_deriv_CI.lower <- as.numeric(resub_CI_lower)
  data$resub_deriv_CI.upper <- as.numeric(resub_CI_upper)
  # Return the updated data frame
  return(data)
}

#' Process All FNMR Data for Derivative Estimation
##' Process All FNMR Data for Derivative Estimation (1st and 2nd Derivatives)
##'
##' This function estimates both the first and second derivatives for all time points in the data list.
##' It returns a named list of data frames, each containing fitted values, first and second derivatives (corrected and uncorrected).
##'
##' @param data_list Named list of data frames (output of load_data and preprocess_shift)
##' @param nseg Number of B-spline segments (default: 35)
##' @param bdeg B-spline degree (default: 4)
##' @param pord Penalty order (default: 2)
##' @param x.grid Grid of x values for prediction (default: data$chem)
##' @return Named list of data frames with fitted values and derivatives

process_all_fnmr <- function(data_list, nseg = 35, bdeg = 4, pord = 2, x.grid = NULL) {
  processed_data <- list()
  tp_names <- names(data_list)
  # Helper to get per-dataset value (by name or index)
  get_param <- function(param, idx, nm) {
    if (length(param) == 1) {
      return(param)
    } else if (!is.null(names(param)) && nm %in% names(param)) {
      return(param[[nm]])
    } else if (length(param) >= idx) {
      return(param[[idx]])
    } else {
      stop(sprintf("Parameter vector too short for dataset %s", nm))
    }
  }
  for (i in seq_along(tp_names)) {
    name <- tp_names[i]
    cat("\nProcessing:", name, "...\n")
    data <- data_list[[name]]
    cat("  Input data: ", nrow(data), "rows\n")
    nseg_i <- get_param(nseg, i, name)
    bdeg_i <- get_param(bdeg, i, name)
    pord_i <- get_param(pord, i, name)
    # Estimate first derivative (r = 1)
    result1 <- estimate_derivative(data, r = 1, nseg = nseg_i, bdeg = bdeg_i, pord = pord_i, x.grid = x.grid, compute_CI = TRUE)
    # Estimate second derivative (r = 2)
    result2 <- estimate_derivative(data, r = 2, nseg = nseg_i, bdeg = bdeg_i, pord = pord_i, x.grid = x.grid, compute_CI = TRUE)
    # Combine results: add second derivative columns to the first derivative data frame
    result1$naive_2nd_deriv_uncorrected <- as.numeric(result2$naive_deriv_uncorrected)
    result1$resub_2nd_deriv_uncorrected <- as.numeric(result2$resub_deriv_uncorrected)
    result1$naive_2nd_deriv <- as.numeric(result2$naive_deriv)
    result1$resub_2nd_deriv <- as.numeric(result2$resub_deriv)
    # Add CIs for 2nd derivative
    result1$naive_2nd_deriv_CI.lower <- as.numeric(result2$naive_deriv_CI.lower)
    result1$naive_2nd_deriv_CI.upper <- as.numeric(result2$naive_deriv_CI.upper)
    result1$resub_2nd_deriv_CI.lower <- as.numeric(result2$resub_deriv_CI.lower)
    result1$resub_2nd_deriv_CI.upper <- as.numeric(result2$resub_deriv_CI.upper)
    cat("  Output columns:", paste(colnames(result1), collapse=", "), "\n")
    cat("  First few fitted values:", paste(head(result1$fitted, 3), collapse=", "), "\n")
    cat("  First few 2nd derivatives (naive, corrected):", paste(head(result1$naive_2nd_deriv, 3), collapse=", "), "\n")
    processed_data[[name]] <- result1
  }
  cat("\nAll time points processed.\n")
  return(processed_data)
}


#' Plot FNMR Mean and Derivatives in Panel Layout (ggplot2, Biometrika style, each panel y-axis)
#'
#' Creates a multi-panel plot for each time point, showing:
#'   - Mean function (fitted, thick black) with original data as grey points
#'   - First derivative (naive: dashed, resub: dotdash, both thin)
#'   - Second derivative (naive: dashed, resub: dotdash, both thin)
#' Panels are ordered: Mean | First Derivative | Second Derivative
#' Rows are ordered: data.0, data.6, ..., data.48
#' Each panel has its own y-axis (ticks and labels).
#' Saves the plot as a PDF to output/plots/FNMR_plots/fnmr_panel_plot.pdf
#'
#' @param processed_results List of processed data frames (from process_all_fnmr)
#' @param time_points Optional vector of time points to plot (default: all)
#' @param xvar Name of x variable (default: 'chem_orig')
#' @param save_pdf Logical, if TRUE saves the plot as PDF
#' @return ggplot2 object (prints the panel plot)
plot_fnmr_panels_gg <- function(processed_results, time_points = NULL, xvar = 'chem_orig', save_pdf = TRUE) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)

  # Set the desired order for time points
  all_time_points <- names(processed_results)
  if (!is.null(time_points)) {
    all_time_points <- intersect(all_time_points, time_points)
  }

  # Gather all data into a tidy data frame
  all_data <- bind_rows(
    lapply(names(processed_results), function(nm) {
      df <- processed_results[[nm]]
      df$time_point <- nm
      df
    }),
    .id = 'id')

  # Filter and set factor order for time_point
  all_data <- all_data %>% filter(time_point %in% all_time_points)
  all_data$time_point <- factor(all_data$time_point, levels = all_time_points)

  # Prepare data for plotting
  plot_data <- all_data %>%
    select(time_point, {{xvar}}, Inten, fitted,
           naive_deriv, resub_deriv,
           naive_2nd_deriv, resub_2nd_deriv) %>%
    pivot_longer(
      cols = c(fitted, naive_deriv, resub_deriv, naive_2nd_deriv, resub_2nd_deriv),
      names_to = 'feature', values_to = 'value') %>%
    mutate(
      panel = factor(case_when(
        feature == 'fitted' ~ 'Mean',
        feature %in% c('naive_deriv', 'resub_deriv') ~ 'First Derivative',
        feature %in% c('naive_2nd_deriv', 'resub_2nd_deriv') ~ 'Second Derivative',
        TRUE ~ feature
      ), levels = c('Mean', 'First Derivative', 'Second Derivative')),
      method = case_when(
        feature == 'fitted' ~ 'Fitted',
        feature == 'naive_deriv' ~ 'Naive',
        feature == 'resub_deriv' ~ 'Resub',
        feature == 'naive_2nd_deriv' ~ 'Naive',
        feature == 'resub_2nd_deriv' ~ 'Resub',
        TRUE ~ feature
      )
    )

  # For mean, only plot one line (no method distinction)
  plot_data$method[plot_data$panel == 'Mean'] <- 'Fitted'

  # For original data points (only for mean panel)
  points_data <- all_data %>%
    select(time_point, {{xvar}}, Inten) %>%
    mutate(panel = factor('Mean', levels = c('Mean', 'First Derivative', 'Second Derivative')))

  # Split data for plotting: mean vs derivatives
  mean_data <- plot_data %>% filter(panel == 'Mean')
  deriv_data <- plot_data %>% filter(panel != 'Mean')

  # Fix: sort by xvar for each group to ensure lines are not broken
  mean_data <- mean_data %>% arrange(time_point, !!sym(xvar))
  deriv_data <- deriv_data %>% arrange(time_point, panel, method, !!sym(xvar))
  points_data <- points_data %>% arrange(time_point, !!sym(xvar))

  # Plot
  gg <- ggplot() +
    # Original data as grey points in mean panel
    geom_point(data = points_data, aes_string(x = xvar, y = 'Inten'), color = 'grey60', size = 0.5, alpha = 0.7) +
    # Fitted mean as thick black line
    geom_line(data = mean_data, aes_string(x = xvar, y = 'value', group = 'time_point'), color = 'black', size = 1.1) +
    # Derivatives as thin black lines, linetype by method
    geom_line(data = deriv_data, aes_string(x = xvar, y = 'value', linetype = 'method', group = 'interaction(time_point, method, panel)'), color = 'black', size = 0.6) +
    facet_grid(time_point ~ panel, scales = 'free_y', switch = 'both') +
    scale_linetype_manual(values = c('Fitted' = 'solid', 'Naive' = 'dashed', 'Resub' = 'dotdash')) +
    labs(x = 'Chemical Shift', y = NULL, linetype = 'Method') +
    theme_bw(base_size = 11) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = 'bold', size = 11),
      legend.position = 'top',
      panel.grid = element_blank(),
      axis.title.y = element_text(angle = 90, vjust = 0.5),
      panel.spacing = unit(0.7, 'lines'),
      axis.text = element_text(color = 'black'),
      axis.ticks = element_line(color = 'black'),
      legend.key = element_blank(),
      panel.border = element_rect(color = 'black', fill = NA)
    )

  if (save_pdf) {
    pdf_dir <- 'output/plots/FNMR_plots'
    if (!dir.exists(pdf_dir)) dir.create(pdf_dir, recursive = TRUE)
    pdf_file <- file.path(pdf_dir, 'fnmr_panel_plot.pdf')
    ggsave(pdf_file, gg, width = 12, height = 12, units = 'in')
    message('Panel plot saved to ', pdf_file)
  }

  print(gg)
  invisible(gg)
}

#' Plot FNMR Mean and Derivatives in a 6x3 Panel Layout (ggplot2)
#'
#' Plots the mean, 1st, and 2nd derivatives for all time points in a 6-column by 3-row layout.
#' Each column is a time point, each row is a feature (mean, 1st, 2nd derivative).
#' @param processed_results List of processed data frames (from process_all_fnmr)
#' @param time_points Optional vector of time points to plot (default: all, in order)
#' @param xvar Name of x variable (default: 'chem_orig')
#' @param save_pdf Logical, if TRUE saves the plot as PDF
#' @return ggplot2 object (prints the panel plot)
plot_fnmr_panels_grid <- function(processed_results, time_points = NULL, xvar = 'chem_orig', save_pdf = TRUE) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)

  # Set the desired order for time points (default: all, in order)
  all_time_points <- names(processed_results)
  if (!is.null(time_points)) {
    all_time_points <- intersect(all_time_points, time_points)
  }
  # For 6x3, select first 6 time points (or all if <6)
  all_time_points <- all_time_points[1:min(6, length(all_time_points))]

  # Gather all data into a tidy data frame
  all_data <- bind_rows(
    lapply(names(processed_results), function(nm) {
      df <- processed_results[[nm]]
      df$time_point <- nm
      df
    }),
    .id = 'id')

  # Filter and set factor order for time_point
  all_data <- all_data %>% filter(time_point %in% all_time_points)
  all_data$time_point <- factor(all_data$time_point, levels = all_time_points)

  # Prepare data for plotting
  plot_data <- all_data %>%
    select(time_point, {{xvar}}, Inten, fitted,
           naive_deriv, resub_deriv,
           naive_2nd_deriv, resub_2nd_deriv) %>%
    pivot_longer(
      cols = c(fitted, naive_deriv, resub_deriv, naive_2nd_deriv, resub_2nd_deriv),
      names_to = 'feature', values_to = 'value') %>%
    mutate(
      feature_row = factor(case_when(
        feature == 'fitted' ~ 'Mean',
        feature %in% c('naive_deriv', 'resub_deriv') ~ 'First Derivative',
        feature %in% c('naive_2nd_deriv', 'resub_2nd_deriv') ~ 'Second Derivative',
        TRUE ~ feature
      ), levels = c('Mean', 'First Derivative', 'Second Derivative')),
      method = case_when(
        feature == 'fitted' ~ 'Fitted',
        feature == 'naive_deriv' ~ 'Naive',
        feature == 'resub_deriv' ~ 'Resub',
        feature == 'naive_2nd_deriv' ~ 'Naive',
        feature == 'resub_2nd_deriv' ~ 'Resub',
        TRUE ~ feature
      )
    )

  # For mean, only plot one line (no method distinction)
  plot_data$method[plot_data$feature_row == 'Mean'] <- 'Fitted'

  # For original data points (only for mean row)
  points_data <- all_data %>%
    select(time_point, {{xvar}}, Inten) %>%
    mutate(feature_row = factor('Mean', levels = c('Mean', 'First Derivative', 'Second Derivative')))

  # Split data for plotting: mean vs derivatives
  mean_data <- plot_data %>% filter(feature_row == 'Mean')
  deriv_data <- plot_data %>% filter(feature_row != 'Mean')

  # Sort for correct line drawing
  mean_data <- mean_data %>% arrange(time_point, !!sym(xvar))
  deriv_data <- deriv_data %>% arrange(time_point, feature_row, method, !!sym(xvar))
  points_data <- points_data %>% arrange(time_point, !!sym(xvar))

  # Plot
  gg <- ggplot() +
    geom_point(data = points_data, aes_string(x = xvar, y = 'Inten'), color = 'grey60', size = 0.5, alpha = 0.7) +
    geom_line(data = mean_data, aes_string(x = xvar, y = 'value', group = 'time_point'), color = 'black', size = 1.1) +
    geom_line(data = deriv_data, aes_string(x = xvar, y = 'value', linetype = 'method', group = 'interaction(time_point, method, feature_row)'), color = 'black', size = 0.6) +
    facet_grid(feature_row ~ time_point, scales = 'free_y', switch = 'both') +
    scale_linetype_manual(values = c('Fitted' = 'solid', 'Naive' = 'dashed', 'Resub' = 'dotdash')) +
    labs(x = 'Chemical Shift', y = NULL, linetype = 'Method') +
    theme_bw(base_size = 11) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = 'bold', size = 11),
      legend.position = 'top',
      panel.grid = element_blank(),
      axis.title.y = element_text(angle = 90, vjust = 0.5),
      panel.spacing = unit(0.7, 'lines'),
      axis.text = element_text(color = 'black'),
      axis.ticks = element_line(color = 'black'),
      legend.key = element_blank(),
      panel.border = element_rect(color = 'black', fill = NA)
    )

  if (save_pdf) {
    pdf_dir <- 'output/plots/FNMR_plots'
    if (!dir.exists(pdf_dir)) dir.create(pdf_dir, recursive = TRUE)
    pdf_file <- file.path(pdf_dir, 'fnmr_panel_grid_plot.pdf')
    ggsave(pdf_file, gg, width = 16, height = 8, units = 'in')
    message('Panel grid plot saved to ', pdf_file)
  }

  print(gg)
  invisible(gg)
}

#' Plot FNMR Mean and Derivatives in an 8x3 Panel Layout (base R, improved layout)
#'
#' Plots the mean, 1st, and 2nd derivatives for all 8 time points in a base R 8-column by 3-row layout.
#' Each column is a time point, each row is a feature (mean, 1st, 2nd derivative).
#' Each plot has an L-shaped box (left and bottom axes only).
#' The first row is the mean function and fitted plot for all 8 hours.
#' Under each, the respective first and second derivatives.
#' For the first derivative, a grey horizontal line at zero is added.
#' @param processed_results List of processed data frames (from process_all_fnmr)
#' @param xvar Name of x variable (default: 'chem_orig')
#' @param save_pdf Logical, if TRUE saves the plot as PDF
#' @return NULL
plot_fnmr_panels_base <- function(processed_results, xvar = 'chem_orig', save_pdf = TRUE) {
  # Ensure time points are sorted by hour (extract numeric from names like 'data.0', 'data.6', ...)
  get_hour <- function(nm) as.numeric(sub('.*\\.', '', nm))
  time_points <- names(processed_results)
  hours <- sapply(time_points, get_hour)
  ord <- order(hours)
  time_points <- time_points[ord]
  hours <- hours[ord]
  n_col <- length(time_points)
  n_row <- 3
  if (save_pdf) {
    pdf_dir <- 'output/plots/FNMR_plots'
    if (!dir.exists(pdf_dir)) dir.create(pdf_dir, recursive = TRUE)
    pdf_file <- file.path(pdf_dir, 'fnmr_panel_base_plot.pdf')
    pdf(pdf_file, width = 2.2*n_col, height = 7)
  }
  old_par <- par(no.readonly = TRUE)
  # Increase left margin for y-axis labels, remove top margin for overall title
  par(mfrow = c(n_row, n_col), mar = c(6, 4.5, 2, 1), oma = c(6, 3.5, 0, 0))
  # Loop over rows (features) and columns (time points)
  for (row in 1:n_row) {
    for (col in 1:n_col) {
      tp <- time_points[col]
      dat <- processed_results[[tp]]
      x <- dat[[xvar]]
      # Only first column gets y-axis label
      ylab <- ''
      if (col == 1) {
        ylab <- switch(row, 'Intensity', '1st Derivative', '2nd Derivative')
      }
      # Only bottom row gets x-axis label
      xlab <- if (row == n_row) 'Chemical Shift' else ''
      # Only top row gets column title
      main <- ''
      if (row == 1) {
        main <- paste0(hours[col], ' Hours')
      }
      if (row == 1) {
        # Mean function
        plot(x, dat$Inten, pch = 16, cex = 0.3, col = 'grey60', main = main, xlab = xlab, ylab = ylab, axes = FALSE, type = 'n')
        points(x, dat$Inten, pch = 16, cex = 0.3, col = 'grey60')
        lines(x, dat$fitted, col = 'black', lwd = 2)
      } else if (row == 2) {
        # 1st Derivative
        plot(x, dat$naive_deriv, type = 'l', lty = 2, col = 'black', lwd = 2.2, main = main, xlab = xlab, ylab = ylab, axes = FALSE)
        lines(x, dat$resub_deriv, lty = 1, col = 'black', lwd = 1)
        abline(h = 0, col = 'grey70', lwd = 1.2)
      } else if (row == 3) {
        # 2nd Derivative
        plot(x, dat$naive_2nd_deriv, type = 'l', lty = 2, col = 'black', lwd = 2.2, main = main, xlab = xlab, ylab = ylab, axes = FALSE)
        lines(x, dat$resub_2nd_deriv, lty = 1, col = 'black', lwd = 1)
      }
      # Axes
      if (row == n_row) {
        axis(1)
        mtext('Chemical Shifts(ppm)', side = 1, line = 4, cex = 1.2)
      } else {
        axis(1, labels = FALSE)
      }
      axis(2)
      box(bty = 'l')
    }
  }
  # Add legend centered below the plot (outer margin)
  par(xpd = NA)
  usr <- par('usr')
  plt <- par('plt')
  # Center legend horizontally in the figure region
  legend_x <- mean(grconvertX(c(plt[1], plt[2]), from = 'nfc', to = 'user'))
  legend_y <- grconvertY(-0.18, from = 'nfc', to = 'user')
  legend(
    x = legend_x, y = legend_y,
    legend = c('Naive', 'Resub'),
    lty = c(2, 1), lwd = c(2.2, 1), col = 'black', horiz = TRUE, bty = 'n',
    xjust = 0.5, yjust = 0,
    seg.len = 2, cex = 1.1, text.col = 'black',
    title = 'Derivative Method'
  )
  par(xpd = FALSE)
  if (save_pdf) dev.off()
  par(old_par)
  message('Base R panel plot saved to ', ifelse(save_pdf, pdf_file, 'screen'))
}

#' Plot FNMR Mean and Derivatives in an 8x3 Panel Layout (base R, with boundary removal option)
#'
#' Like plot_fnmr_panels_base, but allows removal of N boundary points from each end for derivatives.
#' @param processed_results List of processed data frames (from process_all_fnmr)
#' @param xvar Name of x variable (default: 'chem_orig')
#' @param boundary_n Integer, number of points to remove from each end for derivatives (default: 0)
#' @param save_pdf Logical, if TRUE saves the plot as PDF
#' @return NULL
plot_fnmr_panels_base_noboundary <- function(processed_results, xvar = 'chem_orig', boundary_n = 0, save_pdf = TRUE) {
  # Ensure time points are sorted by hour (extract numeric from names like 'data.0', 'data.6', ...)
  get_hour <- function(nm) as.numeric(sub('.*\\.', '', nm))
  time_points <- names(processed_results)
  hours <- sapply(time_points, get_hour)
  ord <- order(hours)
  time_points <- time_points[ord]
  hours <- hours[ord]
  n_col <- length(time_points)
  n_row <- 3
  if (save_pdf) {
    pdf_dir <- 'output/plots/FNMR_plots'
    if (!dir.exists(pdf_dir)) dir.create(pdf_dir, recursive = TRUE)
    pdf_file <- file.path(pdf_dir, 'fnmr_panel_base_plot_noboundary.pdf')
    pdf(pdf_file, width = 2.2*n_col, height = 7)
  }
  old_par <- par(no.readonly = TRUE)
  par(mfrow = c(n_row, n_col), mar = c(3, 4.5, 2, 1), oma = c(4, 3.5, 0, 0))
  # Loop over rows (features) and columns (time points)
  for (row in 1:n_row) {
    for (col in 1:n_col) {
      tp <- time_points[col]
      dat <- processed_results[[tp]]
      x <- dat[[xvar]]
      # Only first column gets y-axis label
      ylab <- ''
      if (col == 1) {
        ylab <- switch(row, 'Intensity', '1st Derivative', '2nd Derivative')
      }
      # Only bottom row gets x-axis label
      xlab <- if (row == n_row) 'Chemical Shift' else ''
      # Only top row gets column title
      main <- ''
      if (row == 1) {
        main <- paste0(hours[col], ' Hours')
      }
      if (row == 1) {
        # Mean function
        plot(x, dat$Inten, pch = 16, cex = 0.3, col = 'grey60', main = main, xlab = xlab, ylab = ylab, axes = FALSE, type = 'n')
        points(x, dat$Inten, pch = 16, cex = 0.3, col = 'grey60')
        lines(x, dat$fitted, col = 'black', lwd = 2)
      } else if (row == 2) {
        # 1st Derivative
        idx <- seq_along(x)
        if (boundary_n > 0 && length(idx) > 2*boundary_n) {
          keep <- idx > boundary_n & idx <= (length(idx) - boundary_n)
        } else {
          keep <- rep(TRUE, length(idx))
        }
        plot(x[keep], dat$naive_deriv[keep], type = 'l', lty = 2, col = 'black', lwd = 2.2, main = main, xlab = xlab, ylab = ylab, axes = FALSE)
        lines(x[keep], dat$resub_deriv[keep], lty = 1, col = 'black', lwd = 1)
        abline(h = 0, col = 'grey70', lwd = 1.2)
      } else if (row == 3) {
        # 2nd Derivative
        idx <- seq_along(x)
        if (boundary_n > 0 && length(idx) > 2*boundary_n) {
          keep <- idx > boundary_n & idx <= (length(idx) - boundary_n)
        } else {
          keep <- rep(TRUE, length(idx))
        }
        plot(x[keep], dat$naive_2nd_deriv[keep], type = 'l', lty = 2, col = 'black', lwd = 2.2, main = main, xlab = xlab, ylab = ylab, axes = FALSE)
        lines(x[keep], dat$resub_2nd_deriv[keep], lty = 1, col = 'black', lwd = 1)
      }
      # Axes
      if (row == n_row) {
        axis(1)
      } else {
        axis(1, labels = FALSE)
      }
      axis(2)
      box(bty = 'l')
    }
  }
  # Add legend centered below the plot (outer margin)
  par(xpd = NA)
  plt <- par('plt')
  legend_x <- mean(grconvertX(c(plt[1], plt[2]), from = 'nfc', to = 'user'))
  legend_y <- grconvertY(-0.18, from = 'nfc', to = 'user')
  legend(
    x = legend_x, y = legend_y,
    legend = c('Naive', 'Resub'),
    lty = c(2, 1), lwd = c(2.2, 1), col = 'black', horiz = TRUE, bty = 'n',
    xjust = 0.5, yjust = 0,
    seg.len = 2, cex = 1.1, text.col = 'black',
    title = 'Derivative Method'
  )
  par(xpd = FALSE)
  if (save_pdf) dev.off()
  par(old_par)
  message('Base R panel plot (noboundary) saved to ', ifelse(save_pdf, pdf_file, 'screen'))
}

#' Plot FNMR Mean and Derivatives in an 8x3 Panel Layout (base R, resub derivatives only, with boundary removal option)
#'
#' Like plot_fnmr_panels_base_noboundary, but only plots the resub derivative estimates (not naive).
#' @param processed_results List of processed data frames (from process_all_fnmr)
#' @param xvar Name of x variable (default: 'chem_orig')
#' @param boundary_n Integer, number of points to remove from each end for derivatives (default: 0)
#' @param save_pdf Logical, if TRUE saves the plot as PDF
plot_fnmr_panels_base_resubonly <- function(processed_results, xvar = 'chem_orig', boundary_n = 0, save_pdf = TRUE) {
  get_hour <- function(nm) as.numeric(sub('.*\\.', '', nm))
  time_points <- names(processed_results)
  hours <- sapply(time_points, get_hour)
  ord <- order(hours)
  time_points <- time_points[ord]
  hours <- hours[ord]
  n_col <- length(time_points)
  n_row <- 3
  if (save_pdf) {
    pdf_dir <- 'output/plots/FNMR_plots'
    if (!dir.exists(pdf_dir)) dir.create(pdf_dir, recursive = TRUE)
    pdf_file <- file.path(pdf_dir, 'fnmr_panel_base_plot_resubonly.pdf')
    pdf(pdf_file, width = 2.2*n_col, height = 7)
  }
  old_par <- par(no.readonly = TRUE)
  par(mfrow = c(n_row, n_col), mar = c(3, 4.5, 2, 1), oma = c(4, 3.5, 0, 0))
  for (row in 1:n_row) {
    for (col in 1:n_col) {
      tp <- time_points[col]
      dat <- processed_results[[tp]]
      x <- dat[[xvar]]
      # Only first column gets y-axis label
      ylab <- ''
      if (col == 1) {
        ylab <- switch(row, 'Intensity', '1st Derivative', '2nd Derivative')
      }
      # Only bottom row gets x-axis label
      xlab <- if (row == n_row) 'Chemical Shift' else ''
      # Only top row gets column title
      main <- ''
      if (row == 1) {
        main <- paste0(hours[col], ' Hours')
      }
      if (row == 1) {
        # Mean function (always plot all points)
        plot(x, dat$Inten, pch = 16, cex = 0.3, col = 'grey60', main = main, xlab = xlab, ylab = ylab, axes = FALSE, type = 'n')
        points(x, dat$Inten, pch = 16, cex = 0.3, col = 'grey60')
        lines(x, dat$fitted, col = 'black', lwd = 2)
      } else if (row == 2) {
        # 1st Derivative (resub only)
        idx <- seq_along(x)
        if (boundary_n > 0 && length(idx) > 2*boundary_n) {
          keep <- idx > boundary_n & idx <= (length(idx) - boundary_n)
        } else {
          keep <- rep(TRUE, length(idx))
        }
        plot(x[keep], dat$resub_deriv[keep], type = 'l', lty = 2, col = 'black', lwd = 2.5, main = main, xlab = xlab, ylab = ylab, axes = FALSE)
        abline(h = 0, col = 'grey70', lwd = 1.5)
      } else if (row == 3) {
        # 2nd Derivative (resub only)
        idx <- seq_along(x)
        if (boundary_n > 0 && length(idx) > 2*boundary_n) {
          keep <- idx > boundary_n & idx <= (length(idx) - boundary_n)
        } else {
          keep <- rep(TRUE, length(idx))
        }
        plot(x[keep], dat$resub_2nd_deriv[keep], type = 'l', lty = 3, col = 'black', lwd = 3, main = main, xlab = xlab, ylab = ylab, axes = FALSE)
      }
      # Axes
      if (row == n_row) {
        axis(1)
      } else {
        axis(1, labels = FALSE)
      }
      axis(2)
      box(bty = 'l')
    }
  }
  # Add legend centered below the plot (outer margin)
  par(xpd = NA)
  plt <- par('plt')
  legend_x <- mean(grconvertX(c(plt[1], plt[2]), from = 'nfc', to = 'user'))
  legend_y <- grconvertY(-0.18, from = 'nfc', to = 'user')
  legend(
    x = legend_x, y = legend_y,
    legend = c('Resub'),
    lty = 3, lwd = 2.2, col = 'black', horiz = TRUE, bty = 'n',
    xjust = 0.5, yjust = 0,
    seg.len = 2, cex = 1.1, text.col = 'black',
    title = 'Derivative Method'
  )
  par(xpd = FALSE)
  if (save_pdf) dev.off()
  par(old_par)
  message('Base R panel plot (resub only) saved to ', ifelse(save_pdf, pdf_file, 'screen'))
}

#' Plot FNMR spectrum with peak marked for each time point
#'
#' For each time point, plot the spectrum (chem_orig vs Inten) and mark the peak (maximum intensity)
#' with a vertical line and annotation. Optionally, return a data.frame of peak positions.
#' @param processed_results List of processed data frames (from process_all_fnmr)
#' @param xvar Name of x variable (default: 'chem_orig')
#' @param save_pdf Logical, if TRUE saves the plot as PDF
#' @return Data frame of peak positions (time point, chem_orig, Inten)
plot_fnmr_peaks <- function(processed_results, xvar = 'chem_orig', save_pdf = TRUE) {
  time_points <- names(processed_results)
  n_col <- 4
  n_row <- ceiling(length(time_points) / n_col)
  if (save_pdf) {
    pdf_dir <- 'output/plots/FNMR_plots'
    if (!dir.exists(pdf_dir)) dir.create(pdf_dir, recursive = TRUE)
    pdf_file <- file.path(pdf_dir, 'fnmr_peaks_plot.pdf')
    pdf(pdf_file, width = 4.5*n_col, height = 3.5*n_row)
  }
  old_par <- par(no.readonly = TRUE)
  par(mfrow = c(n_row, n_col), mar = c(4, 4, 2, 1))
  peak_df <- data.frame(time_point = character(), chem_orig = numeric(), Inten = numeric(), stringsAsFactors = FALSE)
  for (tp in time_points) {
    dat <- processed_results[[tp]]
    x <- dat[[xvar]]
    y <- dat$Inten
    peak_idx <- which.max(y)
    peak_x <- x[peak_idx]
    peak_y <- y[peak_idx]
    plot(x, y, type = 'l', main = tp, xlab = 'Chemical Shift', ylab = 'Intensity')
    abline(v = peak_x, col = 'red', lty = 2, lwd = 2)
    points(peak_x, peak_y, col = 'red', pch = 19)
    text(peak_x, peak_y, labels = sprintf('%.2f', peak_x), pos = 4, col = 'red', cex = 0.9, offset = 0.5)
    peak_df <- rbind(peak_df, data.frame(time_point = tp, chem_orig = peak_x, Inten = peak_y))
  }
  if (save_pdf) dev.off()
  par(old_par)
  message('Peak plots saved to ', ifelse(save_pdf, pdf_file, 'screen'))
  invisible(peak_df)
}

#' Plot zoomed-in derivatives around the true peak (from last time point)
#'
#' Plots the first and second derivatives for all time points, zoomed in around the true peak
#' (peak of the last time point). The x-axis is set to a window around the true peak (e.g., +/- 1 unit).
#' @param processed_results List of processed data frames (from process_all_fnmr)
#' @param xvar Name of x variable (default: 'chem_orig')
#' @param window Numeric, half-width of the window around the true peak (default: 1)
#' @param save_pdf Logical, if TRUE saves the plot as PDF
plot_fnmr_deriv_zoom <- function(processed_results, xvar = 'chem_orig', window = 1, save_pdf = TRUE) {
  # Get the true peak from the last time point
  last_tp <- tail(names(processed_results), 1)
  dat_last <- processed_results[[last_tp]]
  x_last <- dat_last[[xvar]]
  y_last <- dat_last$Inten
  peak_idx <- which.max(y_last)
  true_peak <- x_last[peak_idx]
  xlim <- c(true_peak - window, true_peak + window)
  time_points <- names(processed_results)
  n_row <- 2
  n_col <- length(time_points)
  if (save_pdf) {
    pdf_dir <- 'output/plots/FNMR_plots'
    if (!dir.exists(pdf_dir)) dir.create(pdf_dir, recursive = TRUE)
    pdf_file <- file.path(pdf_dir, 'fnmr_deriv_zoom_plot.pdf')
    pdf(pdf_file, width = 2.2*n_col, height = 5)
  }
  old_par <- par(no.readonly = TRUE)
  par(mfrow = c(n_row, n_col), mar = c(3, 4, 2, 1), oma = c(2, 3, 0, 0))
  for (col in seq_along(time_points)) {
    tp <- time_points[col]
    dat <- processed_results[[tp]]
    x <- dat[[xvar]]
    # First Derivative
    plot(x, dat$naive_deriv, type = 'l', lty = 2, col = 'black', lwd = 2.2, main = paste0(tp, '\n1st Derivative'), xlab = '', ylab = ifelse(col == 1, '1st Derivative', ''), xlim = xlim, axes = FALSE)
    lines(x, dat$resub_deriv, lty = 3, col = 'black', lwd = 2.2)
    abline(h = 0, col = 'grey70', lwd = 1.2)
    axis(1)
    axis(2)
    box(bty = 'l')
    # Second Derivative
    plot(x, dat$naive_2nd_deriv, type = 'l', lty = 2, col = 'black', lwd = 2.2, main = paste0(tp, '\n2nd Derivative'), xlab = ifelse(n_row == 2, 'Chemical Shift', ''), ylab = ifelse(col == 1, '2nd Derivative', ''), xlim = xlim, axes = FALSE)
    lines(x, dat$resub_2nd_deriv, lty = 3, col = 'black', lwd = 2.2)
    abline(v = true_peak, col = 'red', lty = 2, lwd = 1.5)
    axis(1)
    axis(2)
    box(bty = 'l')
  }
  if (save_pdf) dev.off()
  par(old_par)
  message('Zoomed-in derivative plots saved to ', ifelse(save_pdf, pdf_file, 'screen'))
}

#' Create a table comparing the first and second derivative values at the true peak for each time point
#'
#' Returns a data.frame with time point, chem_orig (true peak), first and second derivatives (naive & resub)
#' @param processed_results List of processed data frames (from process_all_fnmr)
#' @param xvar Name of x variable (default: 'chem_orig')
#' @return Data frame with derivative values at the true peak for each time point
get_deriv_at_true_peak <- function(processed_results, xvar = 'chem_orig') {
  # Get the true peak from the last time point
  last_tp <- tail(names(processed_results), 1)
  dat_last <- processed_results[[last_tp]]
  x_last <- dat_last[[xvar]]
  y_last <- dat_last$Inten
  peak_idx <- which.max(y_last)
  true_peak <- x_last[peak_idx]
  # For each time point, find the closest x to the true peak
  res <- data.frame(time_point = character(), chem_orig = numeric(),
                   naive_deriv = numeric(), resub_deriv = numeric(),
                   naive_2nd_deriv = numeric(), resub_2nd_deriv = numeric(),
                   stringsAsFactors = FALSE)
  for (tp in names(processed_results)) {
    dat <- processed_results[[tp]]
    x <- dat[[xvar]]
    idx <- which.min(abs(x - true_peak))
    res <- rbind(res, data.frame(
      time_point = tp,
      chem_orig = x[idx],
      naive_deriv = dat$naive_deriv[idx],
      resub_deriv = dat$resub_deriv[idx],
      naive_2nd_deriv = dat$naive_2nd_deriv[idx],
      resub_2nd_deriv = dat$resub_2nd_deriv[idx]
    ))
  }
  return(res)
}

#' Plot FNMR Mean and Derivatives in an 8x3 Panel Layout (base R, resub derivatives only, with boundary removal option)
#'
#' Like plot_fnmr_panels_base_noboundary, but only plots the resub derivative estimates (not naive).
#' @param processed_results List of processed data frames (from process_all_fnmr)
#' @param xvar Name of x variable (default: 'chem_orig')
#' @param boundary_n Integer, number of points to remove from each end for derivatives (default: 0)
#' @param save_pdf Logical, if TRUE saves the plot as PDF
plot_fnmr_panels_base_resubonly_CI <- function(processed_results, xvar = 'chem_orig', boundary_n = 0, save_pdf = TRUE) {
  get_hour <- function(nm) as.numeric(sub('.*\\.', '', nm))
  time_points <- names(processed_results)
  hours <- sapply(time_points, get_hour)
  ord <- order(hours)
  time_points <- time_points[ord]
  hours <- hours[ord]
  n_col <- length(time_points)
  n_row <- 3
  if (save_pdf) {
    pdf_dir <- 'output/plots/FNMR_plots'
    if (!dir.exists(pdf_dir)) dir.create(pdf_dir, recursive = TRUE)
    pdf_file <- file.path(pdf_dir, 'fnmr_panel_base_plot_resubonly_CI.pdf')
    pdf(pdf_file, width = 2.2*n_col, height = 7)
  }
  old_par <- par(no.readonly = TRUE)
  par(mfrow = c(n_row, n_col), mar = c(3, 4.5, 2, 1), oma = c(4, 3.5, 0, 0))
  row_labels <- c('(a)', '(b)', '(c)')
  for (row in 1:n_row) {
    for (col in 1:n_col) {
      tp <- time_points[col]
      dat <- processed_results[[tp]]
      x <- dat[[xvar]]
      # Only first column gets y-axis label
      ylab <- ''
      if (col == 1) {
        ylab <- switch(row, 'Intensity', '1st Derivative', '2nd Derivative')
      }
      # Only bottom row gets x-axis label
      xlab <- if (row == n_row) 'Chemical Shift' else ''
      # Only top row gets column title
      main <- ''
      if (row == 1) {
        main <- paste0(hours[col], ' Hours')
      }

      # Add (a), (b), (c) to left margin of each row, only for first column, after plot is created
      add_row_label <- function() {
        # Place label in the left outer margin, top of each row, only for first column
        # Increase line for further up, decrease for closer to plot
        mtext(row_labels[row], side = 2, line = 5, adj = 1, las = 1, font = 2, cex = 1.5)
      }
      if (row == 1) {
        # Mean function (always plot all points)
        plot(x, dat$Inten, pch = 5, cex = 0.2, col = 'grey60', main = main, xlab = xlab, ylab = ylab, axes = FALSE, type = 'n')
        points(x, dat$Inten, pch = 5, cex = 0.2, col = 'grey60')
        lines(x, dat$fitted, col = 'black', lwd = 3)
        # Add confidence interval polygon
        if (col == 1) add_row_label()
      } else if (row == 2) {
        # 1st Derivative (resub only, with CI)
        idx <- seq_along(x)
        if (boundary_n > 0 && length(idx) > 2*boundary_n) {
          keep <- idx > boundary_n & idx <= (length(idx) - boundary_n)
        } else {
          keep <- rep(TRUE, length(idx))
        }
        # Remove NA/Inf and sort by x
        xk <- x[keep]
        yk <- dat$resub_deriv[keep]
        ci_low <- dat$resub_deriv_CI.lower[keep]
        ci_up  <- dat$resub_deriv_CI.upper[keep]
        valid <- is.finite(xk) & is.finite(yk) & is.finite(ci_low) & is.finite(ci_up)
        xk <- xk[valid]
        yk <- yk[valid]
        ci_low <- ci_low[valid]
        ci_up  <- ci_up[valid]
        ordx <- order(xk)
        xk <- xk[ordx]
        yk <- yk[ordx]
        ci_low <- ci_low[ordx]
        ci_up <- ci_up[ordx]
        plot(xk, yk, type = 'n', main = main, xlab = xlab, ylab = ylab, axes = FALSE)
        if (length(xk) > 1) {
          polygon(c(xk, rev(xk)), c(ci_low, rev(ci_up)), col = rgb(0.7,0.7,0.7,0.5), border = NA)
        }
        lines(xk, yk, lty = 2, col = 'black', lwd = 3)
        abline(h = 0, col = 'grey70', lwd = 2)
        if (col == 1) add_row_label()
      } else if (row == 3) {
        # 2nd Derivative (resub only, with CI)
        idx <- seq_along(x)
        if (boundary_n > 0 && length(idx) > 2*boundary_n) {
          keep <- idx > boundary_n & idx <= (length(idx) - boundary_n)
        } else {
          keep <- rep(TRUE, length(idx))
        }
        xk <- x[keep]
        yk <- dat$resub_2nd_deriv[keep]
        ci_low <- dat$resub_2nd_deriv_CI.lower[keep]
        ci_up  <- dat$resub_2nd_deriv_CI.upper[keep]
        valid <- is.finite(xk) & is.finite(yk) & is.finite(ci_low) & is.finite(ci_up)
        xk <- xk[valid]
        yk <- yk[valid]
        ci_low <- ci_low[valid]
        ci_up  <- ci_up[valid]
        ordx <- order(xk)
        xk <- xk[ordx]
        yk <- yk[ordx]
        ci_low <- ci_low[ordx]
        ci_up <- ci_up[ordx]
        plot(xk, yk, type = 'n', main = main, xlab = xlab, ylab = ylab, axes = FALSE)
        if (length(xk) > 1) {
          polygon(c(xk, rev(xk)), c(ci_low, rev(ci_up)), col = rgb(0.7,0.7,0.7,0.5), border = NA)
        }
        lines(xk, yk, lty = 5, col = 'black', lwd = 3)
        if (col == 1) add_row_label()
      }
      # Axes
      if (row == n_row) {
        axis(1)
      } else {
        axis(1, labels = FALSE)
      }
      axis(2)
      box(bty = 'l')
    }
  }
  # Add legend centered below the plot (outer margin)
  par(xpd = NA)
  plt <- par('plt')
  legend_x <- mean(grconvertX(c(plt[1], plt[2]), from = 'nfc', to = 'user'))
  legend_y <- grconvertY(-0.18, from = 'nfc', to = 'user')
  legend(
    x = legend_x, y = legend_y,
    legend = c('Resub', '95% CI'),
    lty = c(3, 1), lwd = c(2.2, 6), col = c('black', rgb(0.7,0.7,0.7,0.5)), horiz = TRUE, bty = 'n',
    xjust = 0.5, yjust = 0,
    seg.len = 2, cex = 1.1, text.col = c('black', 'grey40'),
    title = 'Derivative Method'
  )
  par(xpd = FALSE)
  if (save_pdf) dev.off()
  par(old_par)
  message('Base R panel plot (resub only, CI) saved to ', ifelse(save_pdf, pdf_file, 'screen'))
}



# # Debug/test block for pointwiseConfInt (run interactively or as a test)
# if (interactive()) {
#   cat("\n--- Debugging pointwiseConfInt ---\n")
#   # Load and preprocess a single time point
#   fnmr_data <- load_data()
#   fnmr_data <- lapply(fnmr_data, preprocess_shift)
#   data0 <- fnmr_data[[7]]
#   # Estimate resub derivative for r=1
#   naive_result <- naive.est.opt(x = data0$chem, y = data0$Inten, r = 2, nseg = 7, bdeg = 4, pord = 2, x.grid = data0$chem)
#   # Try to compute CI
#   try({
#     ci_out <- pointwiseConfInt(
#       deriv.est = naive_result$fr.hat,
#       r = 2,
#       BS = naive_result$BS,
#       pord = 2,
#       sig.est = naive_result$sig.hat,
#       lambda.est = naive_result$lambda,
#       x = data0$chem
#     )
#     print(str(ci_out))
#     # Add shaded CI background with white plot background
#     op <- par(bg = 'white')
#     plot(data0$chem, naive_result$fr.hat, type = 'n', main = 'Resub Derivative with CI', ylab = 'Derivative', xlab = 'chem')
#     polygon(
#       c(data0$chem, rev(data0$chem)),
#       c(ci_out$CI.lower, rev(ci_out$CI.upper)),
#       col = rgb(0.7, 0.7, 0.9, 0.5), border = NA
#     )
#     lines(data0$chem, naive_result$fr.hat, col = 'black', lwd = 2)
#     lines(data0$chem, ci_out$CI.lower, col = 'red', lty = 2)
#     lines(data0$chem, ci_out$CI.upper, col = 'blue', lty = 2)
#     par(op)
#   })
# }


#' Panel Plot: Resub and Naive Derivative Estimates with CIs (Base R, No Color)
#'
#' Plots intensity, 1st, and 2nd derivatives for all time points, overlaying both resubstitution and naive estimates with their confidence intervals.
#' All features and style are preserved from plot_fnmr_panels_base_resubonly_CI, but now includes naive estimates for direct comparison.
#' @param processed_results List of processed results for each time point
#' @param xvar Name of x variable (default: 'chem_orig')
#' @param boundary_n Integer, number of points to remove from each end for derivatives (default: 0)
#' @param save_pdf Logical, if TRUE saves the plot as PDF
plot_fnmr_panels_base_resubonly_CI_with_naive <- function(processed_results, xvar = 'chem_orig', boundary_n = 0, save_pdf = TRUE) {
  get_hour <- function(nm) as.numeric(sub('.*\\.', '', nm))
  time_points <- names(processed_results)
  hours <- sapply(time_points, get_hour)
  ord <- order(hours)
  time_points <- time_points[ord]
  hours <- hours[ord]
  n_col <- length(time_points)
  n_row <- 3
  if (save_pdf) {
    pdf_dir <- 'output/plots/FNMR_plots'
    if (!dir.exists(pdf_dir)) dir.create(pdf_dir, recursive = TRUE)
    pdf_file <- file.path(pdf_dir, 'fnmr_panel_base_plot_resubonly_CI_with_naive.pdf')
    pdf(pdf_file, width = 2.2*n_col, height = 7)
  }
  old_par <- par(no.readonly = TRUE)
  par(mfrow = c(n_row, n_col), mar = c(3, 4.5, 2, 1), oma = c(4, 3.5, 0, 0))
  row_labels <- c('(a)', '(b)', '(c)')
  for (row in 1:n_row) {
    for (col in 1:n_col) {
      tp <- time_points[col]
      dat <- processed_results[[tp]]
      x <- dat[[xvar]]
      ylab <- ''
      if (col == 1) {
        ylab <- switch(row, 'Intensity', '1st Derivative', '2nd Derivative')
      }
      xlab <- if (row == n_row) 'Chemical Shift' else ''
      main <- ' '
      if (row == 1) {
        main <- paste0(hours[col], ' Hours')
      }
      add_row_label <- function() {
        mtext(row_labels[row], side = 2, line = 5, adj = 1, las = 1, font = 2, cex = 1.5)
      }
      if (row == 1) {
        plot(x, dat$Inten, pch = 5, cex = 0.2, col = 'grey60', main = main, xlab = xlab, ylab = ylab, axes = FALSE, type = 'n')
        points(x, dat$Inten, pch = 5, cex = 0.2, col = 'grey60')
        lines(x, dat$fitted, col = 'black', lwd = 3)
        if (col == 1) add_row_label()
      } else if (row == 2) {
        idx <- seq_along(x)
        if (boundary_n > 0 && length(idx) > 2*boundary_n) {
          keep <- idx > boundary_n & idx <= (length(idx) - boundary_n)
        } else {
          keep <- rep(TRUE, length(idx))
        }
        # Resub estimates
        xk <- x[keep]
        yk <- dat$resub_deriv[keep]
        ci_low <- dat$resub_deriv_CI.lower[keep]
        ci_up  <- dat$resub_deriv_CI.upper[keep]
        # Naive estimates
        yk_naive <- dat$naive_deriv[keep]
        ci_low_naive <- dat$naive_deriv_CI.lower[keep]
        ci_up_naive  <- dat$naive_deriv_CI.upper[keep]
        # Remove NA/Inf and sort by x for both
        valid <- is.finite(xk) & is.finite(yk) & is.finite(ci_low) & is.finite(ci_up) & is.finite(yk_naive) & is.finite(ci_low_naive) & is.finite(ci_up_naive)
        xk <- xk[valid]
        yk <- yk[valid]
        ci_low <- ci_low[valid]
        ci_up <- ci_up[valid]
        yk_naive <- yk_naive[valid]
        ci_low_naive <- ci_low_naive[valid]
        ci_up_naive <- ci_up_naive[valid]
        ordx <- order(xk)
        xk <- xk[ordx]
        yk <- yk[ordx]
        ci_low <- ci_low[ordx]
        ci_up <- ci_up[ordx]
        yk_naive <- yk_naive[ordx]
        ci_low_naive <- ci_low_naive[ordx]
        ci_up_naive <- ci_up_naive[ordx]
        # Compute y limits using all data (resub, naive, and both CIs)
        y_all <- c(yk, yk_naive, ci_low, ci_up, ci_low_naive, ci_up_naive)
        ylim <- range(y_all, na.rm = TRUE)
        plot(xk, yk, type = 'n', main = main, xlab = xlab, ylab = ylab, axes = FALSE, ylim = ylim)
        # Draw CIs first
        if (length(xk) > 1) {
          polygon(c(xk, rev(xk)), c(ci_low_naive, rev(ci_up_naive)), col = rgb(0.85,0.85,0.85,0.5), border = NA)
        }
        if (length(xk) > 1) {
          polygon(c(xk, rev(xk)), c(ci_low, rev(ci_up)), col = rgb(0.7,0.7,0.7,0.5), border = NA)
        }
        # Draw resub first (thicker, dashed)
        lines(xk, yk, lty = 5, col = 'black', lwd = 2)
        # Draw naive on top (thick, solid)
        lines(xk, yk_naive, lty = 1, col = 'black', lwd = 1)
        abline(h = 0, col = 'grey70', lwd = 2)
        if (col == 1) add_row_label()
      } else if (row == 3) {
        idx <- seq_along(x)
        if (boundary_n > 0 && length(idx) > 2*boundary_n) {
          keep <- idx > boundary_n & idx <= (length(idx) - boundary_n)
        } else {
          keep <- rep(TRUE, length(idx))
        }
        # Resub 2nd derivative
        xk <- x[keep]
        yk <- dat$resub_2nd_deriv[keep]
        ci_low <- dat$resub_2nd_deriv_CI.lower[keep]
        ci_up  <- dat$resub_2nd_deriv_CI.upper[keep]
        # Naive 2nd derivative
        yk_naive <- dat$naive_2nd_deriv[keep]
        ci_low_naive <- dat$naive_2nd_deriv_CI.lower[keep]
        ci_up_naive  <- dat$naive_2nd_deriv_CI.upper[keep]
        valid <- is.finite(xk) & is.finite(yk) & is.finite(ci_low) & is.finite(ci_up) & is.finite(yk_naive) & is.finite(ci_low_naive) & is.finite(ci_up_naive)
        xk <- xk[valid]
        yk <- yk[valid]
        ci_low <- ci_low[valid]
        ci_up <- ci_up[valid]
        yk_naive <- yk_naive[valid]
        ci_low_naive <- ci_low_naive[valid]
        ci_up_naive <- ci_up_naive[valid]
        ordx <- order(xk)
        xk <- xk[ordx]
        yk <- yk[ordx]
        ci_low <- ci_low[ordx]
        ci_up <- ci_up[ordx]
        yk_naive <- yk_naive[ordx]
        ci_low_naive <- ci_low_naive[ordx]
        ci_up_naive <- ci_up_naive[ordx]
        # Compute y limits using all data (resub, naive, and both CIs)
        y_all <- c(yk, yk_naive, ci_low, ci_up, ci_low_naive, ci_up_naive)
        ylim <- range(y_all, na.rm = TRUE)
        plot(xk, yk, type = 'n', main = main, xlab = xlab, ylab = ylab, axes = FALSE, ylim = ylim)
        # Draw CIs first
        if (length(xk) > 1) {
          polygon(c(xk, rev(xk)), c(ci_low_naive, rev(ci_up_naive)), col = rgb(0.85,0.85,0.85,0.5), border = NA)
        }
        if (length(xk) > 1) {
          polygon(c(xk, rev(xk)), c(ci_low, rev(ci_up)), col = rgb(0.7,0.7,0.7,0.5), border = NA)
        }
        # Draw resub first (thicker, dotted)
        lines(xk, yk, lty = 5, col = 'black', lwd = 2)
        # Draw naive on top (thick, solid)
        lines(xk, yk_naive, lty = 1, col = 'black', lwd = 1)
        if (col == 1) add_row_label()
      }
      if (row == n_row) {
        axis(1)
      } else {
        axis(1, labels = FALSE)
      }
      axis(2)
      box(bty = 'l')
    }
  }
  par(xpd = NA)
  plt <- par('plt')
  legend_x <- grconvertX(0.5, from = 'nfc', to = 'user')
  legend_y <- grconvertY(0.98, from = 'nfc', to = 'user')
  legend(
    x = legend_x, y = legend_y,
    legend = c('Resub', 'Naive', '95% CI', 'Maxima', 'True Max'),
    lty = c(5, 1, NA, NA, NA), lwd = c(2, 1, 6, NA, NA),
    col = c('black', 'black', rgb(0.7,0.7,0.7,0.5), 1, 2),
    pch = c(NA, NA, NA, 8, 1),
    horiz = TRUE, bty = 'n',
    xjust = 0.5, yjust = 1,
    seg.len = 2, cex = 1.1, text.col = c('black', 'black', 'grey40', 'black', 'black'),
    title = 'Derivative Method'
  )
  par(xpd = FALSE)
  if (save_pdf) dev.off()
  par(old_par)
  message('Base R panel plot (resub + naive, CI) saved to ', ifelse(save_pdf, pdf_file, 'screen'))
}


plot_fnmr_panels_maxima_detection <- function(processed_results, xvar = 'chem_orig', boundary_n = 0, save_pdf = TRUE, true_max_shift = -38) {
  get_hour <- function(nm) as.numeric(sub('.*\\.', '', nm))
  time_points <- names(processed_results)
  hours <- sapply(time_points, get_hour)
  ord <- order(hours)
  time_points <- time_points[ord]
  hours <- hours[ord]
  n_col <- length(time_points)
  n_row <- 4
  
  if (save_pdf) {
    pdf_dir <- 'output/plots/FNMR_plots'
    if (!dir.exists(pdf_dir)) dir.create(pdf_dir, recursive = TRUE)
    pdf_file <- file.path(pdf_dir, 'fnmr_panel_maxima_detection.pdf')
    pdf(pdf_file, width = 2.2*n_col, height = 9)
  }
  
  old_par <- par(no.readonly = TRUE)
  
  # Modified margins and outer margins
  par(mfrow = c(n_row, n_col), 
      mar = c(4, 4.5, 2, 1),       # Increased bottom margin from 3 to 4
      oma = c(14, 3.5, 4, 0),      # Increased bottom outer margin from 12 to 14
      mgp = c(2, 0.5, 0))          # Adjusted axis label positioning
  
  row_labels <- c('(a)', '(b)', '(c)', '(d)')

  for (row in 1:n_row) {
    for (col in 1:n_col) {
      tp <- time_points[col]
      dat <- processed_results[[tp]]
      x <- dat[[xvar]]
      ylab <- ''
      if (col == 1) {
        ylab <- switch(row, 'Intensity', '1st Derivative', '2nd Derivative', 'Intensity')
      }
      xlab <- if (row == 4) 'Chemical Shifts (ppm)' else ''
      main <- if (row == 1) paste0(hours[col], ' Hours') else ' '

      add_row_label <- function() {
        mtext(row_labels[row], side = 2, line = 5, adj = 1, las = 1, font = 2, cex = 1.5)
      }

      if (row == 1) {
        plot(x, dat$Inten, pch = 5, cex = 0.2, col = 'grey60', main = main, 
             xlab = '', ylab = ylab, axes = FALSE, type = 'n')
        points(x, dat$Inten, pch = 5, cex = 0.2, col = 'grey60')
        lines(x, dat$fitted, col = 'black', lwd = 3)
        axis(1, labels = FALSE); axis(2); box(bty = 'l')
        if (col == 1) add_row_label()
      } else if (row %in% c(2, 3)) {
        deriv_type <- if (row == 2) 'resub_deriv' else 'resub_2nd_deriv'
        naive_type <- if (row == 2) 'naive_deriv' else 'naive_2nd_deriv'
        ci_low <- dat[[paste0(deriv_type, '_CI.lower')]]
        ci_up  <- dat[[paste0(deriv_type, '_CI.upper')]]
        ci_low_naive <- dat[[paste0(naive_type, '_CI.lower')]]
        ci_up_naive  <- dat[[paste0(naive_type, '_CI.upper')]]

        idx <- seq_along(x)
        keep <- if (boundary_n > 0 && length(idx) > 2 * boundary_n) idx > boundary_n & idx <= (length(idx) - boundary_n) else rep(TRUE, length(idx))
        xk <- x[keep]; yk <- dat[[deriv_type]][keep]; yk_naive <- dat[[naive_type]][keep]
        ci_low <- ci_low[keep]; ci_up <- ci_up[keep]
        ci_low_naive <- ci_low_naive[keep]; ci_up_naive <- ci_up_naive[keep]
        valid <- is.finite(xk) & is.finite(yk) & is.finite(ci_low) & is.finite(ci_up) & is.finite(yk_naive) & is.finite(ci_low_naive) & is.finite(ci_up_naive)
        xk <- xk[valid]; yk <- yk[valid]; yk_naive <- yk_naive[valid]
        ci_low <- ci_low[valid]; ci_up <- ci_up[valid]
        ci_low_naive <- ci_low_naive[valid]; ci_up_naive <- ci_up_naive[valid]
        ordx <- order(xk); xk <- xk[ordx]; yk <- yk[ordx]; yk_naive <- yk_naive[ordx]
        ci_low <- ci_low[ordx]; ci_up <- ci_up[ordx]
        ci_low_naive <- ci_low_naive[ordx]; ci_up_naive <- ci_up_naive[ordx]

        y_all <- c(yk, yk_naive, ci_low, ci_up, ci_low_naive, ci_up_naive)
        ylim <- range(y_all, na.rm = TRUE)

        plot(xk, yk, type = 'n', main = main, xlab = '', ylab = ylab, axes = FALSE, ylim = ylim)
        polygon(c(xk, rev(xk)), c(ci_low_naive, rev(ci_up_naive)), col = rgb(0.85,0.85,0.85,0.5), border = NA)
        polygon(c(xk, rev(xk)), c(ci_low, rev(ci_up)), col = rgb(0.7,0.7,0.7,0.5), border = NA)
        lines(xk, yk, lty = 5, col = 'black', lwd = 2)
        lines(xk, yk_naive, lty = 1, col = 'black', lwd = 1)
        abline(h = 0, col = 'grey70', lwd = 2)
        axis(1, labels = FALSE); axis(2); box(bty = 'l')
        if (col == 1) add_row_label()
      } else if (row == 4) {
        idx <- seq_along(x)
        keep <- if (boundary_n > 0 && length(idx) > 2 * boundary_n) idx > boundary_n & idx <= (length(idx) - boundary_n) else rep(TRUE, length(idx))
        xk <- x[keep]; yk <- dat$fitted[keep]; fd <- dat$resub_deriv[keep]; sd <- dat$resub_2nd_deriv[keep]
        valid <- is.finite(xk) & is.finite(yk) & is.finite(fd) & is.finite(sd)
        xk <- xk[valid]; yk <- yk[valid]; fd <- fd[valid]; sd <- sd[valid]
        ordx <- order(xk); xk <- xk[ordx]; yk <- yk[ordx]; fd <- fd[ordx]; sd <- sd[ordx]

        zero_crossings <- which(diff(sign(fd)) != 0)
        maxima_idx <- zero_crossings[sd[zero_crossings] < 0]
        maxima_shifts <- xk[maxima_idx]; maxima_vals <- yk[maxima_idx]

        plot(xk, yk, type = 'l', main = main, xlab = xlab, ylab = ylab, axes = FALSE)
        points(maxima_shifts, maxima_vals, pch = 8, col = 1)
        tallest_idx <- which.max(maxima_vals)
        for (i in seq_along(maxima_shifts)) {
          text(maxima_shifts[i], maxima_vals[i], labels = round(maxima_shifts[i], 2), 
               pos = if (i == tallest_idx) 2 else 3, offset = if (i == tallest_idx) 1 else 0, 
               cex = if (i == tallest_idx) 0.9 else 0.7, font = if (i == tallest_idx) 2 else 1)
        }
        points(true_max_shift, yk[which.min(abs(xk - true_max_shift))], pch = 1, cex = 2, lwd = 2)
        axis(1); axis(2); box(bty = 'l')
        if (col == 1) add_row_label()
      }
    }
  }

  # Modified legend positioning
  par(xpd = NA)  # Allow drawing outside plot area
  par(fig = c(0, 1, 0, 0.05), new = TRUE, mar = c(0, 0, 0, 0))  # Reduced height and zero margins
  plot.new()
  legend(
    "center",
    legend = c('Resub', 'Naive', '95% CI', 'Maxima', 'True Max'),
    lty = c(5, 1, NA, NA, NA), 
    lwd = c(2, 1, 6, NA, NA),
    col = c('black', 'black', rgb(0.7, 0.7, 0.7, 0.5), 1, 2),
    pch = c(NA, NA, 15, 8, 1),  # Using pch=15 for CI for better visibility
    horiz = TRUE, 
    bty = 'n',
    seg.len = 2, 
    cex = 1.1,
    text.col = c('black', 'black', 'grey40', 'black', 'black'),
    title = 'Derivative Method'
  )
  par(xpd = FALSE)

  if (save_pdf) dev.off()
  par(old_par)
  message('Base R panel plot (maxima detection) saved to ', ifelse(save_pdf, pdf_file, 'screen'))
}




### Geeting exact zero crossings for maxima detection with interpolation
#' Plot FNMR Maxima Detection with Interpolated Zero Crossings
# plot_fnmr_panels_maxima_detection_interp <- function(processed_results, xvar = 'chem_orig', boundary_n = 0, save_pdf = TRUE, true_max_shift = -38) {
#   get_hour <- function(nm) as.numeric(sub('.*\\.', '', nm))
#   time_points <- names(processed_results)
#   hours <- sapply(time_points, get_hour)
#   ord <- order(hours)
#   time_points <- time_points[ord]
#   hours <- hours[ord]
#   n_col <- length(time_points)
#   n_row <- 4
  
#   if (save_pdf) {
#     pdf_dir <- 'output/plots/FNMR_plots'
#     if (!dir.exists(pdf_dir)) dir.create(pdf_dir, recursive = TRUE)
#     pdf_file <- file.path(pdf_dir, 'fnmr_panel_maxima_detection_interp.pdf')
#     pdf(pdf_file, width = 2.2*n_col, height = 9)
#   }
  
#   old_par <- par(no.readonly = TRUE)
#   par(mfrow = c(n_row, n_col), 
#       mar = c(4, 4.5, 2, 1),
#       oma = c(14, 3.5, 4, 0),
#       mgp = c(2, 0.5, 0))
  
#   row_labels <- c('(a)', '(b)', '(c)', '(d)')

#   for (row in 1:n_row) {
#     for (col in 1:n_col) {
#       tp <- time_points[col]
#       dat <- processed_results[[tp]]
#       x <- dat[[xvar]]
#       ylab <- ''
#       if (col == 1) {
#         ylab <- switch(row, 'Intensity', '1st Derivative', '2nd Derivative', 'Intensity')
#       }
#       xlab <- if (row == 4) 'Chemical Shifts (ppm)' else ''
#       main <- if (row == 1) paste0(hours[col], ' Hours') else ' '

#       add_row_label <- function() {
#         mtext(row_labels[row], side = 2, line = 5, adj = 1, las = 1, font = 2, cex = 1.5)
#       }

#       if (row == 1) {
#         plot(x, dat$Inten, pch = 5, cex = 0.2, col = 'grey60', main = main, 
#              xlab = '', ylab = ylab, axes = FALSE, type = 'n')
#         points(x, dat$Inten, pch = 5, cex = 0.2, col = 'grey60')
#         lines(x, dat$fitted, col = 'black', lwd = 3)
#         axis(1, labels = FALSE); axis(2); box(bty = 'l')
#         if (col == 1) add_row_label()
#       } else if (row %in% c(2, 3)) {
#         deriv_type <- if (row == 2) 'resub_deriv' else 'resub_2nd_deriv'
#         naive_type <- if (row == 2) 'naive_deriv' else 'naive_2nd_deriv'
#         ci_low <- dat[[paste0(deriv_type, '_CI.lower')]]
#         ci_up  <- dat[[paste0(deriv_type, '_CI.upper')]]
#         ci_low_naive <- dat[[paste0(naive_type, '_CI.lower')]]
#         ci_up_naive  <- dat[[paste0(naive_type, '_CI.upper')]]

#         idx <- seq_along(x)
#         keep <- if (boundary_n > 0 && length(idx) > 2 * boundary_n) idx > boundary_n & idx <= (length(idx) - boundary_n) else rep(TRUE, length(idx))
#         xk <- x[keep]; yk <- dat[[deriv_type]][keep]; yk_naive <- dat[[naive_type]][keep]
#         ci_low <- ci_low[keep]; ci_up <- ci_up[keep]
#         ci_low_naive <- ci_low_naive[keep]; ci_up_naive <- ci_up_naive[keep]
#         valid <- is.finite(xk) & is.finite(yk) & is.finite(ci_low) & is.finite(ci_up) & is.finite(yk_naive) & is.finite(ci_low_naive) & is.finite(ci_up_naive)
#         xk <- xk[valid]; yk <- yk[valid]; yk_naive <- yk_naive[valid]
#         ci_low <- ci_low[valid]; ci_up <- ci_up[valid]
#         ci_low_naive <- ci_low_naive[valid]; ci_up_naive <- ci_up_naive[valid]
#         ordx <- order(xk); xk <- xk[ordx]; yk <- yk[ordx]; yk_naive <- yk_naive[ordx]
#         ci_low <- ci_low[ordx]; ci_up <- ci_up[ordx]
#         ci_low_naive <- ci_low_naive[ordx]; ci_up_naive <- ci_up_naive[ordx]

#         y_all <- c(yk, yk_naive, ci_low, ci_up, ci_low_naive, ci_up_naive)
#         ylim <- range(y_all, na.rm = TRUE)

#         plot(xk, yk, type = 'n', main = main, xlab = '', ylab = ylab, axes = FALSE, ylim = ylim)
#         polygon(c(xk, rev(xk)), c(ci_low_naive, rev(ci_up_naive)), col = rgb(0.85,0.85,0.85,0.5), border = NA)
#         polygon(c(xk, rev(xk)), c(ci_low, rev(ci_up)), col = rgb(0.7,0.7,0.7,0.5), border = NA)
#         lines(xk, yk, lty = 5, col = 'black', lwd = 2)
#         lines(xk, yk_naive, lty = 1, col = 'black', lwd = 1)
#         abline(h = 0, col = 'grey70', lwd = 2)
#         axis(1, labels = FALSE); axis(2); box(bty = 'l')
#         if (col == 1) add_row_label()
#       } else if (row == 4) {
#         idx <- seq_along(x)
#         keep <- if (boundary_n > 0 && length(idx) > 2 * boundary_n) idx > boundary_n & idx <= (length(idx) - boundary_n) else rep(TRUE, length(idx))
#         xk <- x[keep]; yk <- dat$fitted[keep]; fd <- dat$resub_deriv[keep]; sd <- dat$resub_2nd_deriv[keep]
#         valid <- is.finite(xk) & is.finite(yk) & is.finite(fd) & is.finite(sd)
#         xk <- xk[valid]; yk <- yk[valid]; fd <- fd[valid]; sd <- sd[valid]
#         ordx <- order(xk); xk <- xk[ordx]; yk <- yk[ordx]; fd <- fd[ordx]; sd <- sd[ordx]

#         # Enhanced maxima detection with interpolation
#         zero_crossings <- which(diff(sign(fd)) != 0)
#         maxima_shifts <- sapply(zero_crossings, function(i) {
#           if (i >= length(xk)) return(NA)  # Skip edge cases
#           x1 <- xk[i]
#           x2 <- xk[i + 1]
#           y1 <- fd[i]
#           y2 <- fd[i + 1]
#           # Linear interpolation for sub-resolution accuracy
#           x_interp <- x1 - y1 * (x2 - x1) / (y2 - y1)
#           return(x_interp)
#         })
#         maxima_idx <- which(!is.na(maxima_shifts) & sd[zero_crossings] < 0)
#         maxima_shifts <- maxima_shifts[maxima_idx]
#         maxima_vals <- approx(xk, yk, xout = maxima_shifts)$y  # Interpolated intensities

#         plot(xk, yk, type = 'l', main = main, xlab = xlab, ylab = ylab, axes = FALSE)
#         points(maxima_shifts, maxima_vals, pch = 8, col = 1)
#         tallest_idx <- which.max(maxima_vals)
#         for (i in seq_along(maxima_shifts)) {
#           text(maxima_shifts[i], maxima_vals[i], labels = round(maxima_shifts[i], 2), 
#                pos = if (i == tallest_idx) 2 else 3, offset = if (i == tallest_idx) 1 else 0, 
#                cex = if (i == tallest_idx) 0.9 else 0.7, font = if (i == tallest_idx) 2 else 1)
#         }
#         points(true_max_shift, yk[which.min(abs(xk - true_max_shift))], pch = 1, cex = 2, lwd = 2)
#         axis(1); axis(2); box(bty = 'l')
#         if (col == 1) add_row_label()
#       }
#     }
#   }

#   # Legend
#   par(xpd = NA)
#   par(fig = c(0, 1, 0, 0.05), new = TRUE, mar = c(0, 0, 0, 0))
#   plot.new()
#   legend(
#     "center",
#     legend = c('Resub', 'Naive', '95% CI', 'Maxima', 'True Max'),
#     lty = c(5, 1, NA, NA, NA), 
#     lwd = c(2, 1, 6, NA, NA),
#     col = c('black', 'black', rgb(0.7, 0.7, 0.7, 0.5), 1, 2),
#     pch = c(NA, NA, 15, 8, 1),
#     horiz = TRUE, 
#     bty = 'n',
#     seg.len = 2, 
#     cex = 1.1,
#     text.col = c('black', 'black', 'grey40', 'black', 'black'),
#     title = 'Derivative Method'
#   )
#   par(xpd = FALSE)

#   if (save_pdf) dev.off()
#   par(old_par)
#   message('Interpolated maxima detection plot saved to ', ifelse(save_pdf, pdf_file, 'screen'))
# }

plot_fnmr_panels_maxima_detection_interp <- function(processed_results, xvar = 'chem_orig', boundary_n = 0, save_pdf = TRUE, true_max_shift = -38, use_interpolation = TRUE) {
  get_hour <- function(nm) as.numeric(sub('.*\\.', '', nm))
  time_points <- names(processed_results)
  hours <- sapply(time_points, get_hour)
  ord <- order(hours)
  time_points <- time_points[ord]
  hours <- hours[ord]
  n_col <- length(time_points)
  n_row <- 4
  
  if (save_pdf) {
    pdf_dir <- 'output/plots/FNMR_plots'
    if (!dir.exists(pdf_dir)) dir.create(pdf_dir, recursive = TRUE)
    pdf_file <- file.path(pdf_dir, 'fnmr_panel_maxima_detection_interp.pdf')
    pdf(pdf_file, width = 2.2*n_col, height = 10)
  }
  
  old_par <- par(no.readonly = TRUE)
  par(mfrow = c(n_row, n_col), 
      mar = c(5, 4.5, 2, 1),
      oma = c(1, 3.5, 4, 0),
      mgp = c(2, 0.5, 0))
  
  row_labels <- c('(a)', '(b)', '(c)', '(d)')

  for (row in 1:n_row) {
    for (col in 1:n_col) {
      tp <- time_points[col]
      dat <- processed_results[[tp]]
      x <- dat[[xvar]]
      ylab <- ''
      if (col == 1) {
        ylab <- switch(row, 'Intensity', '1st Derivative', '2nd Derivative', 'Intensity')
      }
      xlab <- if (row == 4) 'Chemical Shifts (ppm)' else ''
      main <- if (row == 1) paste0(hours[col], ' Hours') else ' '

      add_row_label <- function() {
        mtext(row_labels[row], side = 2, line = 5, adj = 1, las = 1, font = 2, cex = 1.5)
      }

      if (row == 1) {
        plot(x, dat$Inten, pch = 5, cex = 0.2, col = 'grey60', main = main, 
             xlab = '', ylab = ylab, axes = FALSE, type = 'n')
        points(x, dat$Inten, pch = 5, cex = 0.2, col = 'grey60')
        lines(x, dat$fitted, col = 'black', lwd = 3)
        axis(1, labels = FALSE); axis(2); box(bty = 'l')
        if (col == 1) add_row_label()
      } else if (row %in% c(2, 3)) {
        deriv_type <- if (row == 2) 'resub_deriv' else 'resub_2nd_deriv'
        naive_type <- if (row == 2) 'naive_deriv' else 'naive_2nd_deriv'
        ci_low <- dat[[paste0(deriv_type, '_CI.lower')]]
        ci_up  <- dat[[paste0(deriv_type, '_CI.upper')]]
        ci_low_naive <- dat[[paste0(naive_type, '_CI.lower')]]
        ci_up_naive  <- dat[[paste0(naive_type, '_CI.upper')]]

        idx <- seq_along(x)
        keep <- if (boundary_n > 0 && length(idx) > 2 * boundary_n) idx > boundary_n & idx <= (length(idx) - boundary_n) else rep(TRUE, length(idx))
        xk <- x[keep]; yk <- dat[[deriv_type]][keep]; yk_naive <- dat[[naive_type]][keep]
        ci_low <- ci_low[keep]; ci_up <- ci_up[keep]
        ci_low_naive <- ci_low_naive[keep]; ci_up_naive <- ci_up_naive[keep]
        valid <- is.finite(xk) & is.finite(yk) & is.finite(ci_low) & is.finite(ci_up) & is.finite(yk_naive) & is.finite(ci_low_naive) & is.finite(ci_up_naive)
        xk <- xk[valid]; yk <- yk[valid]; yk_naive <- yk_naive[valid]
        ci_low <- ci_low[valid]; ci_up <- ci_up[valid]
        ci_low_naive <- ci_low_naive[valid]; ci_up_naive <- ci_up_naive[valid]
        ordx <- order(xk); xk <- xk[ordx]; yk <- yk[ordx]; yk_naive <- yk_naive[ordx]
        ci_low <- ci_low[ordx]; ci_up <- ci_up[ordx]
        ci_low_naive <- ci_low_naive[ordx]; ci_up_naive <- ci_up_naive[ordx]

        y_all <- c(yk, yk_naive, ci_low, ci_up, ci_low_naive, ci_up_naive)
        ylim <- range(y_all, na.rm = TRUE)

        plot(xk, yk, type = 'n', main = main, xlab = '', ylab = ylab, axes = FALSE, ylim = ylim)
        polygon(c(xk, rev(xk)), c(ci_low_naive, rev(ci_up_naive)), col = rgb(0.85,0.85,0.85,0.5), border = NA)
        polygon(c(xk, rev(xk)), c(ci_low, rev(ci_up)), col = rgb(0.7,0.7,0.7,0.5), border = NA)
        lines(xk, yk, lty = 5, col = 'black', lwd = 2)
        lines(xk, yk_naive, lty = 1, col = 'black', lwd = 1)
        abline(h = 0, col = 'grey70', lwd = 2)
        axis(1, labels = FALSE); axis(2); box(bty = 'l')
        if (col == 1) add_row_label()
      } else if (row == 4) {
        idx <- seq_along(x)
        keep <- if (boundary_n > 0 && length(idx) > 2 * boundary_n) idx > boundary_n & idx <= (length(idx) - boundary_n) else rep(TRUE, length(idx))
        xk <- x[keep]; yk <- dat$fitted[keep]
        fd_resub <- dat$resub_deriv[keep]; sd_resub <- dat$resub_2nd_deriv[keep]
        fd_naive <- dat$naive_deriv[keep]; sd_naive <- dat$naive_2nd_deriv[keep]
        
        valid <- is.finite(xk) & is.finite(yk) & 
                 is.finite(fd_resub) & is.finite(sd_resub) &
                 is.finite(fd_naive) & is.finite(sd_naive)
        xk <- xk[valid]; yk <- yk[valid]
        fd_resub <- fd_resub[valid]; sd_resub <- sd_resub[valid]
        fd_naive <- fd_naive[valid]; sd_naive <- sd_naive[valid]
        ordx <- order(xk); xk <- xk[ordx]; yk <- yk[ordx]
        fd_resub <- fd_resub[ordx]; sd_resub <- sd_resub[ordx]
        fd_naive <- fd_naive[ordx]; sd_naive <- sd_naive[ordx]

        # Find maxima with option for interpolation or approximation
        find_maxima <- function(fd, sd, xk, yk) {
          zero_crossings <- which(diff(sign(fd)) != 0)
          
          if (use_interpolation) {
            maxima_shifts <- sapply(zero_crossings, function(i) {
              if (i >= length(xk)) return(NA)
              x1 <- xk[i]
              x2 <- xk[i + 1]
              y1 <- fd[i]
              y2 <- fd[i + 1]
              x1 - y1 * (x2 - x1) / (y2 - y1)  # Linear interpolation
            })
          } else {
            maxima_shifts <- xk[zero_crossings]  # Simple approximation
          }
          
          maxima_idx <- which(!is.na(maxima_shifts) & sd[zero_crossings] < 0)
          maxima_shifts <- maxima_shifts[maxima_idx]
          maxima_vals <- approx(xk, yk, xout = maxima_shifts)$y
          list(shifts = maxima_shifts, vals = maxima_vals)
        }
        
        resub_maxima <- find_maxima(fd_resub, sd_resub, xk, yk)
        naive_maxima <- find_maxima(fd_naive, sd_naive, xk, yk)

        # Adjust ylim to accommodate labels
        y_range <- range(yk)
        y_expansion <- 0.15 * diff(y_range)
        ylim <- c(y_range[1], y_range[2] + y_expansion)

        plot(xk, yk, type = 'l', main = main, xlab = xlab, ylab = ylab, 
             axes = FALSE, ylim = ylim)
        
        # Plot all resub maxima (black plus signs)
        points(resub_maxima$shifts, resub_maxima$vals, pch = 3, cex=1.5, lwd=1)
        
        # Plot all naive maxima (grey40 Asterisks)
        points(naive_maxima$shifts, naive_maxima$vals, pch = 2, cex=1, lwd=1)

        # Find maxima near true max (±1.5 ppm)
        true_max_val <- yk[which.min(abs(xk - true_max_shift))]
        resub_near <- which(abs(resub_maxima$shifts - true_max_shift) < 2.5)
        naive_near <- which(abs(naive_maxima$shifts - true_max_shift) < 2.5)
        
        # Label resub maxima (right side, black)
        if (length(resub_near) > 0) {
          text(resub_maxima$shifts[resub_near], resub_maxima$vals[resub_near], 
               labels = round(resub_maxima$shifts[resub_near], 2), 
               pos = 4, cex = 1, col = 'black')
        }
        
        # Label naive maxima (left side, grey)
        if (length(naive_near) > 0) {
          text(naive_maxima$shifts[naive_near], naive_maxima$vals[naive_near], 
               labels = round(naive_maxima$shifts[naive_near], 2), 
               pos = 2, cex = 0.8, col = 'grey40')
        }
        
        # Plot true max (black open circle)
        points(true_max_shift, true_max_val, 
               pch = 1, cex = 3, lwd = 1, col = 'black')
        
        axis(1); axis(2); box(bty = 'l')
        if (col == 1) add_row_label()
      }
    }
  }

  # Legend positioned below
  par(xpd = NA)
  par(fig = c(0, 1, 0, 0.01), new = TRUE, mar = c(0, 0, 0, 0))
  plot.new()
  legend(
    "center",
    legend = c('Resub', 'Naive', '95% CI', 'Resub Max', 'Naive Max', 'True Max'),
    lty = c(5, 1, NA, NA, NA, NA), 
    lwd = c(2, 1, 6, 1, 1, 1),
    col = c('black', 'black', rgb(0.7, 0.7, 0.7, 0.5), "black", "grey40", "black"),
    pch = c(NA, NA, 15, 3, 2, 1),
    horiz = TRUE, 
    bty = 'b',
    seg.len = 2.5, 
    cex = 1.0,
    text.col = 'black'
  )
  par(xpd = FALSE)

  if (save_pdf) dev.off()
  par(old_par)
  message('Enhanced maxima detection plot saved to ', ifelse(save_pdf, pdf_file, 'screen'))
}