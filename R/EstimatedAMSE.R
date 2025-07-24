# ============================================================================
# File: EstimatedAMSE.R
# Purpose: Compare estimated and true MISE (Mean Integrated Squared Error)
#          for one-step and iterated P-spline estimators, and generate plots.
# Author: [Christopher Odoom]
# Date: October 29, 2024 (last updated: May 1, 2025)
#
# Description:
#   - Simulates data and applies P-spline smoothing using one-step and iterated methods.
#   - Computes estimated and true MISE for a range of smoothing parameters (lambda).
#   - Aggregates results over multiple simulations.
#   - Visualizes the results, including confidence intervals and minimum MISE locations.
#
# Usage:
#   - Source this script in R or run interactively.
#   - Requires: ggplot2, dplyr, gridExtra, and custom functions in RegSplineFunctions.R.
# ============================================================================

# --- Load required functions and libraries ---
library(here)
library(ggplot2)
library(dplyr)
library(cowplot)
library(gridExtra)
source("R/MainFunctions.R") # Custom spline functions

# --- Simulation and Estimation Functions ---
## Iterated
sig.fun <- function(sig, f, x, x.grid, r, nseg, pord, bdeg, lambdas, n, tol = 1e-10, ITs = 10) {
  
  y <- f + sig * rnorm(n)
  resub.fit <- resub.est(x = x, y = y, r = r, x.grid = x.grid, nseg = nseg, pord = pord,
                     bdeg = bdeg, tol = tol, ITs = ITs)
  naive.fit <- naive.est.opt(x = x, y = y, r = r, x.grid = x.grid, 
                         nseg = nseg, pord = pord, bdeg = bdeg)
  plugin.fit <- plugin.est(x=x, y=y, r=r,nseg=nseg , pord = pord,
                        bdeg=bdeg, x.grid=x,fr="Y") # nolint
  est.sig <- naive.fit$sig.hat


  mise_os <- sapply(lambdas, function(l) {
    mise.lambda.optim(lambda = l, x = x, y = y, r = r, sig = est.sig, nseg = nseg, pord = pord, 
                bdeg = bdeg, f = naive.fit$f.hat, fr = plugin.fit$fr.hat)
  })

  mise_it <- sapply(lambdas, function(l) {
    mise.lambda.optim(lambda = l, x = x.grid, y = y, r = r, sig = est.sig, nseg = nseg, pord = pord, 
                bdeg = bdeg, f = naive.fit$f.hat, fr = resub.fit$fr.hat)
  })

  return(list(mise_os = mise_os, mise_it = mise_it))
}


## This function uses sig.fun and sig.fun2
nsim_mise_all <- function(nsims, sigma, f, fr, x, x.grid, r, nseg, pord, bdeg, lambdas, n, 
                          tol = 1e-10, ITs = 10) {
  results1 <- list()
  results2 <- list()
  for (nsim in 1:nsims) {
    ## one step
    mise <- sig.fun(sig = sigma, f = f, x = x, x.grid = x.grid, r = r, nseg = nseg, pord = pord, 
                    bdeg = bdeg, lambdas = lambdas, n = n, tol = tol, ITs = ITs)

    results1[[paste0("nsim_", nsim)]] <- data.frame(lambda = lambdas, mise = mise$mise_os, nsim = nsim, r = r)
    results2[[paste0("nsim_", nsim)]] <- data.frame(lambda = lambdas, mise = mise$mise_it, nsim = nsim, r = r)
  }
  mise_true <- sapply(lambdas, function(l) {
    mise.lambda.optim(lambda = l, x = x, y = y, r = r, sig = sigma, nseg = nseg, pord = pord, 
                bdeg = bdeg, f = f, fr = fr)
  })

  results1[[paste0("nsim_", nsim + 1)]] <- data.frame(lambda = lambdas, mise = mise_true, nsim = nsim + 1, r = r)
  results2[[paste0("nsim_", nsim + 1)]] <- data.frame(lambda = lambdas, mise = mise_true, nsim = nsim + 1, r = r)

  results_df1 <- do.call(rbind, results1)
  results_df2 <- do.call(rbind, results2)

  return(list(resdf1 = results_df1, resdf2 = results_df2))
}

# --- Data Generation and Parameter Setup ---
### Define data
set.seed(1500)
n <-200
x <- seq(0, 1, length.out=n)
f<- 32 * exp(-8 * (1 - 2 * x)^2) * (1 - 2 * x)
ra <- range(f)[2]-range(f)[1]
ra
f.prime<- (4096 * x^2 - 4096 * x + 960) * exp(-8 * (1 - 2 * x)^2)
f.pprime  <- -(262144*x^3-393216*x^2+184320*x-26624)*exp(-8*(1-2*x)^2)
x.grid <- x

# Parameters for iterated
sigma <- 0.1
r<- 1
logl <- seq(-1,2.5, length.out=500)
lambdas <- 10^(logl)  # example lambda values
nseg <- 35
pord <-2
bdeg <- 4
nsims<-5
res_df<-nsim_mise_all(nsims=nsims,sigma=sigma,f=f,fr=f.prime, x=x, x.grid=x.grid, r=r, nseg=nseg,
                    pord=pord, bdeg=bdeg, lambdas=lambdas, n=n
                  ,tol=1e-10, ITs=100)

head(res_df$resdf1)
head(res_df$resdf2)
# Save the results to a file
# save(res_df, file = "rdd.RData")
#
# Load the results from the file
# load("rdd.RData")
# res_df <- get(load("rdd.RData"))

step1.data<- res_df$resdf1
# y<- step1.data %>% filter(nsim==21)
# plot(lambdas,y$mise)
step2.data<- res_df$resdf2

# --- Plotting: Compare Estimated and True MISE ---
########################## Manseging the plot ########################################

# Calculate minimum lambda for each simulation in each dataset
step1_min_lambda <- step1.data %>%filter(nsim!=(nsims+1)) %>% 
  group_by(nsim) %>%
  slice_min(log10(mise), with_ties = FALSE) %>%
  summarize(min_lambda = lambda)

step2_min_lambda <- step2.data %>%filter(nsim!=(nsims+1)) %>%
  group_by(nsim) %>%
  slice_min(log10(mise), with_ties = FALSE) %>%
  summarize(min_lambda = lambda)

# Calculate mean, standard error, and 95% confidence interval for the minimum lambdas
step1_mean_min_lambda <- mean(step1_min_lambda$min_lambda)
step1_se_lambda <- sd(step1_min_lambda$min_lambda) 
step1_ci_lambda <- c(step1_mean_min_lambda - 1.96 * step1_se_lambda, step1_mean_min_lambda + 1.96 * step1_se_lambda)

step2_mean_min_lambda <- mean(step2_min_lambda$min_lambda)
step2_se_lambda <- sd(step2_min_lambda$min_lambda) 
step2_ci_lambda <- c(step2_mean_min_lambda - 1.96 * step2_se_lambda, step2_mean_min_lambda + 1.96 * step2_se_lambda)

# Convert confidence intervals to log10 scale for plotting
step1_ci_lambda_log <- log10(step1_ci_lambda)
step2_ci_lambda_log <- log10(step2_ci_lambda)

# Set common x-axis limits to include confidence interval lines
common_x_limits <- range(log10(step1.data$lambda), log10(step2.data$lambda), step1_ci_lambda_log, step2_ci_lambda_log)

# Create the first plot with dashed vertical lines for confidence interval
png(here::here("output", "plots","estAMSE_plots", "est1.png"), width = 800, height = 400)
# Create the first plot with dashed vertical lines for confidence interval
p1 <- step1.data %>%
  ggplot(aes(x = log10(lambda), y = log10(mise), group = factor(nsim))) +
  geom_line(aes(color = factor(nsim == nsims+1))) +
  scale_color_manual(
    values = c("TRUE" = "black", "FALSE" = "grey"),
    labels = c("FALSE" = "estimate", "TRUE" = "truth"),
    name = ""
  ) +
  geom_vline(xintercept = step1_ci_lambda_log, linetype = "dashed", color = "blue") + # Add CI lines
  theme_minimal() + 
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.border = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    panel.grid = element_blank()
  ) +
  labs(x = expression(log[10](lambda)), y = expression(log[10](hat(MISE)))) +
  ggtitle("One Step") +
  xlim(common_x_limits)

#Create the second plot with dashed vertical lines for confidence interval
p2 <- step2.data %>%
  ggplot(aes(x = log10(lambda), y = log10(mise), group = factor(nsim))) +
  geom_line(aes(color = factor(nsim == nsims+1))) +
  scale_color_manual(
    values = c("TRUE" = "black", "FALSE" = "grey"),
    labels = c("FALSE" = "estimate", "TRUE" = "truth"),
    name = ""
  ) +
  geom_vline(xintercept = step2_ci_lambda_log, linetype = "dashed", color = "blue") + # Add CI lines
  theme_minimal() + 
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.border = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    panel.grid = element_blank(),
    axis.title.y = element_blank() # Remove y label for the right plot
  ) +
  labs(x = expression(log[10](lambda))) +
  ggtitle("Iterated") +
  xlim(common_x_limits)

# Arrange the plots side by side with a single y-axis label on the left and more space
grid.arrange(p1, p2, ncol = 2, widths = c(1.2, 1), left = expression(log[10](hat(MISE))))
dev.off()


################################### Plotting with CI lines  with different colors for nsim #####################################
png(here::here("output", "plots","estAMSE_plots", "est2.png"), width = 800, height = 400)
p1 <- step1.data %>%
  ggplot(aes(x = log10(lambda), y = log10(mise), group = factor(nsim))) +
  geom_line(aes(color = factor(nsim))) +
  geom_vline(xintercept = step1_ci_lambda_log, linetype = "dashed", color = "blue") + # Add CI lines
  theme_minimal() + 
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.border = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    panel.grid = element_blank()
  ) +
  labs(x = expression(log[10](lambda)), y = expression(log[10](hat(MISE)))) +
  ggtitle("One Step") +
  xlim(common_x_limits)
#Create the second plot with dashed vertical lines for confidence interval
p2 <- step2.data %>%
  ggplot(aes(x = log10(lambda), y = log10(mise), group = factor(nsim))) +
  geom_line(aes(color = factor(nsim))) +
  geom_vline(xintercept = step2_ci_lambda_log, linetype = "dashed", color = "blue") + # Add CI lines
  theme_minimal() + 
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.border = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    panel.grid = element_blank(),
    axis.title.y = element_blank() # Remove y label for the right plot
  ) +
  labs(x = expression(log[10](lambda))) +
  ggtitle("Iterated") +
  xlim(common_x_limits)

# Arrange the plots side by side with a single y-axis label on the left
grid.arrange(p1, p2, ncol = 2, left = expression(log[10](hat(MISE))))
dev.off()



# --- Additional Examples and Notes ---
# The code below is for adding a rug to the plot

######### Adding a rug ###########

# Calculate minimum lambda and corresponding mise for each simulation in each dataset
step1_min_lambda <- step1.data %>%
  group_by(nsim) %>%
  slice_min(log10(mise), with_ties = FALSE) %>%
  mutate(is_truth = (nsim == nsims + 1))  # TRUE if it's the "truth" line

step2_min_lambda <- step2.data %>%
  group_by(nsim) %>%
  slice_min(log10(mise), with_ties = FALSE) %>%
  mutate(is_truth = (nsim == nsims + 1))  #

# Convert minimum lambdas to log10 scale for plotting
step1_min_lambda_log <- log10(step1_min_lambda$lambda)
step2_min_lambda_log <- log10(step2_min_lambda$lambda)

# Set common x-axis limits
common_x_limits <- range(log10(step1.data$lambda), log10(step2.data$lambda), step1_min_lambda_log, step2_min_lambda_log)

# Create the first plot with rug for minimum lambda
png(here::here("output", "plots","estAMSE_plots", "emse.hat.png"), width = 800, height = 400)
p1 <- step1.data %>%
  ggplot(aes(x = log10(lambda), y = log10(mise), group = factor(nsim))) +
  geom_line(aes(color = factor(nsim == nsims + 1))) +
  scale_color_manual(
    values = c("TRUE" = "black", "FALSE" = "red"),
    labels = c("FALSE" = "estimate", "TRUE" = "truth"),
    name = ""
  ) +
  geom_rug(data = step1_min_lambda, aes(x = log10(lambda),color= factor(is_truth)),
           sides = "b") +
  theme_minimal() +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.border = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    panel.grid = element_blank()
  ) +
  labs(x = expression(log[10](lambda)), y = expression(log[10](hat(MISE)))) +
  ggtitle("One Step") +
  xlim(common_x_limits)

# Create the second plot with rug for minimum lambda
p2 <- step2.data %>%
  ggplot(aes(x = log10(lambda), y = log10(mise), group = factor(nsim))) +
  geom_line(aes(color = factor(nsim == nsims + 1))) +
  scale_color_manual(
    values = c("TRUE" = "black", "FALSE" = "red"),
    labels = c("FALSE" = "estimate", "TRUE" = "truth"),
    name = ""
  ) +
  geom_rug(data = step2_min_lambda, aes(x = log10(lambda),color= factor(is_truth)),
           sides = "b")+
  theme_minimal() +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.border = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    panel.grid = element_blank(),
    axis.title.y = element_blank() # Remove y label for the right plot
  ) +
  labs(x = expression(log[10](lambda))) +
  ggtitle("Iterated") +
  xlim(common_x_limits)

# Arrange the plots side by side with a single y-axis label on the left
grid.arrange(p1, p2, ncol = 2, left = expression(log[10](hat(MISE))))
dev.off()

# End of CompareMISE_PSplines.R



##----- Compare MISE  Iterated P-Spline Estimators and GCV for n=100, 1000 for two different noise levels------

# This code compares the MISE of iterated P-spline estimators and GCV for different sample sizes (n=100, n=1000) and noise levels.
## Generate data


# Function to generate results for different sample sizes and noise levels
compare_sample_sizes <- function(sample_sizes = c(100, 1000), noise_levels = c(0.1, 0.5)) {
  results_list <- list()
  
  for (n_idx in seq_along(sample_sizes)) {
    for (sig_idx in seq_along(noise_levels)) {
      n <- sample_sizes[n_idx]
      sigma <- noise_levels[sig_idx]
      
      # Generate data
      x <- seq(0, 1, length.out = n)
      f <- 32 * exp(-8 * (1 - 2 * x)^2) * (1 - 2 * x)
      f.prime <- (4096 * x^2 - 4096 * x + 960) * exp(-8 * (1 - 2 * x)^2)
      x.grid <- x
      
      # Set parameters
      r <- 1
      logl <- seq(-1, 2.5, length.out = 200)
      lambdas <- 10^(logl)
      nseg <- min(35, floor(n/5))  # Adjust number of segments based on sample size
      pord <- 2
      bdeg <- 4
      nsims <- 5
      
      # Run simulations
      cat(sprintf("Running simulations for n=%d, sigma=%.2f\n", n, sigma))
      res <- nsim_mise_all(nsims = nsims, sigma = sigma, f = f, fr = f.prime, 
                          x = x, x.grid = x.grid, r = r, nseg = nseg,
                          pord = pord, bdeg = bdeg, lambdas = lambdas, n = n,
                          tol = 1e-10, ITs = 50)
      
      # Store results
      results_list[[paste0("n", n, "_sig", sigma)]] <- res
      
      # Run GCV for comparison
      gcv_results <- lapply(1:nsims, function(sim) {
        y <- f + sigma * rnorm(n)
        naive_fit <- naive.est.opt(x = x, y = y, r = r, x.grid = x.grid, 
                               nseg = nseg, pord = pord, bdeg = bdeg)
        data.frame(
          n = n,
          sigma = sigma,
          sim = sim,
          lambda_gcv = naive_fit$lambda,
          mise_gcv = mean((naive_fit$fr.hat - f.prime)^2)
        )
      })
      
      results_list[[paste0("gcv_n", n, "_sig", sigma)]] <- do.call(rbind, gcv_results)
    }
  }
  
  return(results_list)
}

# Run the comparison
comparison_results <- compare_sample_sizes()
tail(comparison_results)





# Create visualizations
plot_comparison <- function(results = comparison_results, 
                           sample_sizes = c(1000), 
                           noise_levels = c(0.1, 0.5),
                           output_dir = "output/plots/estAMSE_plots") {
                            library(here)
library(ggplot2)
library(dplyr)
library(cowplot)
library(gridExtra)
  # Setup
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!requireNamespace("ggplot2", quietly = TRUE) || 
      !requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("cowplot", quietly = TRUE)) {
    stop("Please install required packages: ggplot2, dplyr, cowplot")
  }
  
  # Create LaTeX footnote text file
  footnote_text <- paste(
    "The figure compares MISE curves across different regularization parameters $\\lambda$.",
    "\\textbf{Solid gray lines} show estimated MISE curves from iterative computations.",
    "\\textbf{Solid black line} represents the true MISE curve.",
    "\\textbf{Blue dashed line}: true optimal $\\lambda$. \\textbf{Red dashed}: GCV-selected $\\lambda$. \\textbf{Green dotted}: Iterated P-spline $\\lambda$. \\textbf{Red point}: GCV result ($\\lambda$, MISE).",
    "All quantities are plotted on logarithmic scales (base 10)."
  )
  writeLines(footnote_text, file.path(output_dir, "figure_footnote.tex"))
  
  # Setup panel tags
  panel_tags <- c("(a)", "(b)")
  plot_list <- list()
  plot_count <- 1
  
  # Main loop
  for (n in sample_sizes) {
    for (sigma in noise_levels) {
      # Extract data
      result_key <- paste0("n", n, "_sig", sigma)
      gcv_key <- paste0("gcv_n", n, "_sig", sigma)
      
      # Skip if data not available
      if (!result_key %in% names(results) || !gcv_key %in% names(results)) next
      
      iterated_data <- results[[result_key]]$resdf2
      gcv_data <- results[[gcv_key]]
      
      # Skip if data is empty
      if (is.null(iterated_data) || nrow(iterated_data) == 0 || 
          is.null(gcv_data) || nrow(gcv_data) == 0) next
      
      # Get true MISE data
      max_nsim <- max(iterated_data$nsim, na.rm = TRUE)
      if (is.na(max_nsim)) next
      
      true_mise_data <- iterated_data[iterated_data$nsim == max_nsim, ]
      if (nrow(true_mise_data) == 0) next
      
      min_mise_idx <- which.min(true_mise_data$mise)
      if (length(min_mise_idx) == 0) next
      
      true_mise <- true_mise_data[min_mise_idx, ]
      
      # Calculate mean values
      mean_gcv_lambda <- mean(gcv_data$lambda_gcv, na.rm = TRUE)
      mean_gcv_mise <- mean(gcv_data$mise_gcv, na.rm = TRUE)
      if (is.na(mean_gcv_lambda) || is.na(mean_gcv_mise)) next
      
      # Get iterated lambda
      iterated_min_lambda <- iterated_data %>%
        dplyr::filter(nsim != max_nsim) %>%
        dplyr::group_by(nsim) %>%
        dplyr::slice_min(mise, with_ties = FALSE) %>%
        dplyr::summarize(min_lambda = mean(lambda, na.rm = TRUE))
      
      mean_iterated_lambda <- mean(iterated_min_lambda$min_lambda, na.rm = TRUE)
      
      # Create the plot
      p <- ggplot2::ggplot() +
        # Estimated MISE lines (iterations)
        ggplot2::geom_line(
          data = iterated_data[iterated_data$nsim != max_nsim,],
          ggplot2::aes(x = log10(lambda), y = log10(mise), group = factor(nsim), color = "Estimated MISE"),
          alpha = 0.7, linewidth = 0.5
        ) +
        # True MISE line
        ggplot2::geom_line(
          data = true_mise_data,
          ggplot2::aes(x = log10(lambda), y = log10(mise), color = "True MISE"),
          linewidth = 1.2
        ) +
        # Vertical lines for different lambdas
        ggplot2::geom_vline(xintercept = log10(true_mise$lambda), 
                           linetype = "dashed", color = "blue", linewidth = 1, show.legend = TRUE) +
        ggplot2::geom_vline(xintercept = log10(mean_gcv_lambda), 
                           linetype = "dashed", color = "red", linewidth = 1, show.legend = TRUE) +
        ggplot2::geom_vline(xintercept = log10(mean_iterated_lambda), 
                           linetype = "dotted", color = "purple", linewidth = 1, show.legend = TRUE) +
        # GCV result point
        ggplot2::geom_point(
          ggplot2::aes(x = log10(mean_gcv_lambda), y = log10(mean_gcv_mise)),
          color = "#D55E00", fill = "white", shape = 21, size = 3, stroke = 1, show.legend = TRUE
        ) +
        # Color scale
        ggplot2::scale_color_manual(
          name = NULL,
          values = c("Estimated MISE" = "grey60", "True MISE" = "black"),
          labels = c("Estimated MISE" = "Estimated MISE (iterations)", "True MISE" = "True MISE")
        ) +
        # Theme settings
        ggplot2::theme_bw(base_size = 13) +
        ggplot2::theme(
          panel.grid = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(hjust = 0.5, size = 13),
          axis.title = ggplot2::element_text(size = 12),
          axis.text = ggplot2::element_text(size = 10),
          panel.border = ggplot2::element_blank(),
          axis.line.x.bottom = ggplot2::element_line(color = "black", linewidth = 0.7),
          axis.line.y.left = ggplot2::element_line(color = "black", linewidth = 0.7),
          axis.line.y.right = ggplot2::element_blank(),
          axis.line.x.top = ggplot2::element_blank(),
          plot.margin = ggplot2::unit(c(2, 2, 2, 2), "mm"),
          legend.position = if (plot_count == 2) "none" else "none"
        ) +
        # Labels
        ggplot2::labs(
          x = if (plot_count > 2) expression(log[10](lambda)) else NULL,
          y = if (plot_count %% 2 == 1) expression(log[10](MISE)) else NULL,
          title = bquote(n == .(n) ~ ", " ~ sigma == .(sigma))
        )
      
      # Add panel tag
      p <- cowplot::ggdraw(p) + 
           cowplot::draw_plot_label(panel_tags[plot_count], 
                                  x = 0.04, y = 0.96, 
                                  hjust = 0, vjust = 1, 
                                  fontface = "bold", size = 14)
      
      # Add to list
      plot_list[[plot_count]] <- p
      plot_count <- plot_count + 1
    }
  }
  
  # Arrange and save
  combined <- cowplot::plot_grid(plotlist = plot_list, ncol = 2, align = "hv")
  pdf_file <- file.path(output_dir, "mise_comparison_all_panels.pdf")
  cowplot::save_plot(pdf_file, combined, base_width = 6, base_height = 4)
  
  message(sprintf("Combined panel plot saved to %s", pdf_file))
  message("LaTeX footnote saved as figure_footnote.tex")
  
  return(combined) # Return the plot object
}

# Generate and save the comparison plots
comparison_plot <- plot_comparison(comparison_results)
# Summarize results

summarize_comparison <- function(results, sample_sizes = c(100, 1000), noise_levels = c(0.05, 0.2)) {
  summary_df <- data.frame(
    sample_size = integer(),
    noise_level = numeric(),
    true_optimal_lambda = numeric(),
    mean_gcv_lambda = numeric(),
    mean_iterated_lambda = numeric(),
    true_optimal_mise = numeric(),
    mean_gcv_mise = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (n in sample_sizes) {
    for (sigma in noise_levels) {
      # Extract data
      result_key <- paste0("n", n, "_sig", sigma)
      gcv_key <- paste0("gcv_n", n, "_sig", sigma)
      
      # Skip if data not available
      if (!result_key %in% names(results) || !gcv_key %in% names(results)) next
      
      iterated_data <- results[[result_key]]$resdf2
      gcv_data <- results[[gcv_key]]
      
      # Skip if data is empty
      if (is.null(iterated_data) || nrow(iterated_data) == 0 || 
          is.null(gcv_data) || nrow(gcv_data) == 0) next
      
      # Get minimum MISE values from true function
      true_mise <- iterated_data %>% 
        filter(nsim == max(nsim)) %>%
        slice_min(mise, with_ties = FALSE)
      
      # Get mean optimal lambda from iterated method
      iterated_min_lambda <- iterated_data %>%
        filter(nsim != max(nsim)) %>%
        group_by(nsim) %>%
        slice_min(mise, with_ties = FALSE) %>%
        summarize(min_lambda = mean(lambda))
      
      # Calculate mean values
      mean_gcv_lambda <- mean(gcv_data$lambda_gcv)
      mean_gcv_mise <- mean(gcv_data$mise_gcv)
      
      # Add to summary dataframe
      summary_df <- rbind(summary_df, data.frame(
        sample_size = n,
        noise_level = sigma,
        true_optimal_lambda = true_mise$lambda,
        mean_gcv_lambda = mean_gcv_lambda,
        mean_iterated_lambda = iterated_min_lambda$min_lambda,
        true_optimal_mise = true_mise$mise,
        mean_gcv_mise = mean_gcv_mise
      ))
    }
  }
  
  print(summary_df)
  return(summary_df)
}

# Generate and print summary statistics
comparison_summary <- summarize_comparison(comparison_results)
print("Summary of comparison between iterated P-splines and GCV:")
print(comparison_summary)

# Calculate percentage improvements
comparison_summary$lambda_error_gcv <- abs(comparison_summary$mean_gcv_lambda - comparison_summary$true_optimal_lambda) / comparison_summary$true_optimal_lambda * 100
comparison_summary$lambda_error_iterated <- abs(comparison_summary$mean_iterated_lambda - comparison_summary$true_optimal_lambda) / comparison_summary$true_optimal_lambda * 100
print("Error percentages in lambda selection:")
print(comparison_summary[, c("sample_size", "noise_level", "lambda_error_gcv", "lambda_error_iterated")])

# Save the results for future reference
save(comparison_results, comparison_summary, file = here::here("output", "data", "estAMSE_data", "mise_comparison_results.RData"))








# Create grayscale visualizations for Biometrika paper
plot_comparison_grayscale <- function(results = comparison_results, 
                                     sample_sizes = c(1000), 
                                     noise_levels = c(0.1, 0.5),
                                     output_dir = "output/plots/estAMSE_plots_grayscale") {
  # Setup
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!requireNamespace("ggplot2", quietly = TRUE) || 
      !requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("cowplot", quietly = TRUE)) {
    stop("Please install required packages: ggplot2, dplyr, cowplot")
  }
  
  # Create LaTeX footnote text file for grayscale version
  footnote_text <- paste(
    "The figure compares MISE curves across different regularization parameters $\\lambda$.",
    "\\textbf{Thin gray lines} show estimated MISE curves from iterative computations.",
    "\\textbf{Thick solid line} represents the true MISE curve.",
    "\\textbf{Dashed line}: true optimal $\\lambda$. \\textbf{Dash-dot line}: GCV-selected $\\lambda$. \\textbf{Dotted line}: Iterated P-spline $\\lambda$. \\textbf{Circle}: GCV result ($\\lambda$, MISE).",
    "All quantities are plotted on logarithmic scales (base 10)."
  )
  writeLines(footnote_text, file.path(output_dir, "figure_footnote_grayscale.tex"))
  
  # Setup panel tags
  panel_tags <- c("(a)", "(b)")
  plot_list <- list()
  plot_count <- 1
  
  # First pass: collect y-axis ranges for consistent scaling
  all_y_values <- c()
  
  for (n in sample_sizes) {
    for (sigma in noise_levels) {
      result_key <- paste0("n", n, "_sig", sigma)
      gcv_key <- paste0("gcv_n", n, "_sig", sigma)
      
      if (!result_key %in% names(results) || !gcv_key %in% names(results)) next
      
      iterated_data <- results[[result_key]]$resdf2
      gcv_data <- results[[gcv_key]]
      
      if (is.null(iterated_data) || nrow(iterated_data) == 0 || 
          is.null(gcv_data) || nrow(gcv_data) == 0) next
      
      # Collect y values for scaling
      all_y_values <- c(all_y_values, log10(iterated_data$mise))
      all_y_values <- c(all_y_values, log10(gcv_data$mise_gcv))
    }
  }
  
  # Calculate common y-axis limits
  if (length(all_y_values) > 0) {
    y_range <- range(all_y_values, na.rm = TRUE)
    y_margin <- diff(y_range) * 0.05  # 5% margin
    common_y_limits <- c(y_range[1] - y_margin, y_range[2] + y_margin)
  } else {
    common_y_limits <- NULL
  }
  
  # Main loop
  for (n in sample_sizes) {
    for (sigma in noise_levels) {
      # Extract data
      result_key <- paste0("n", n, "_sig", sigma)
      gcv_key <- paste0("gcv_n", n, "_sig", sigma)
      
      # Skip if data not available
      if (!result_key %in% names(results) || !gcv_key %in% names(results)) next
      
      iterated_data <- results[[result_key]]$resdf2
      gcv_data <- results[[gcv_key]]
      
      # Skip if data is empty
      if (is.null(iterated_data) || nrow(iterated_data) == 0 || 
          is.null(gcv_data) || nrow(gcv_data) == 0) next
      
      # Get true MISE data
      max_nsim <- max(iterated_data$nsim, na.rm = TRUE)
      if (is.na(max_nsim)) next
      
      true_mise_data <- iterated_data[iterated_data$nsim == max_nsim, ]
      if (nrow(true_mise_data) == 0) next
      
      min_mise_idx <- which.min(true_mise_data$mise)
      if (length(min_mise_idx) == 0) next
      
      true_mise <- true_mise_data[min_mise_idx, ]
      
      # Calculate mean values
      mean_gcv_lambda <- mean(gcv_data$lambda_gcv, na.rm = TRUE)
      mean_gcv_mise <- mean(gcv_data$mise_gcv, na.rm = TRUE)
      if (is.na(mean_gcv_lambda) || is.na(mean_gcv_mise)) next
      
      # Get iterated lambda
      iterated_min_lambda <- iterated_data %>%
        filter(nsim != max(nsim)) %>%
        group_by(nsim) %>%
        slice_min(mise, with_ties = FALSE) %>%
        summarize(min_lambda = lambda)
      
      mean_iterated_lambda <- mean(iterated_min_lambda$min_lambda, na.rm = TRUE)
      
      # Create the grayscale plot
      p <- ggplot2::ggplot() +
        # Plot simulation MISE curves
        ggplot2::geom_line(
          data = iterated_data %>% filter(nsim != max(nsim)),
          ggplot2::aes(x = log10(lambda), y = log10(mise), group = factor(nsim)),
          color = "gray70", linewidth = 0.5, alpha = 0.7
        ) +
        # Plot true MISE curve
        ggplot2::geom_line(
          data = true_mise_data,
          ggplot2::aes(x = log10(lambda), y = log10(mise)),
          color = "black", linewidth = 1, alpha = 0.9
        ) +
        # Vertical lines for different lambdas with specified line types and widths
        ggplot2::geom_vline(xintercept = log10(true_mise$lambda), 
                           linetype = 1, linewidth = 1, color = "black") +  # lty=1, lwd=1: thin solid
        ggplot2::geom_vline(xintercept = log10(mean_gcv_lambda), 
                           linetype = 1, linewidth = 1, color = "gray30") +  # lty=1, lwd=3: thick solid
        ggplot2::geom_vline(xintercept = log10(mean_iterated_lambda), 
                           linetype = 4, linewidth = 1, color = "black") +   # lty=4, lwd=3: thick dashed
        # GCV result point
        ggplot2::geom_point(
          data = data.frame(x = log10(mean_gcv_lambda), y = log10(mean_gcv_mise)),
          ggplot2::aes(x = x, y = y),
          color = "black", fill = "white", shape = 21, size = 3, stroke = 1
        ) +
        # Theme settings
        ggplot2::theme_bw(base_size = 13) +
        ggplot2::theme(
          legend.position = "none",
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white"),
          plot.title = element_text(hjust = 0.5, size = 14),
          axis.title = element_text(size = 13),
          axis.text = element_text(size = 12)
        ) +
        # Labels
        ggplot2::labs(
          title = sprintf("n = %d, σ = %.2f", n, sigma),
          x = expression(log[10](lambda)),
          y = if (plot_count == 1) expression(log[10](MISE)) else ""
        )
      
      # Apply consistent y-axis limits if available
      if (!is.null(common_y_limits)) {
        p <- p + ggplot2::ylim(common_y_limits)
      }
      
      # Add panel tag with updated size of 2.5
      p <- cowplot::ggdraw(p) + 
           cowplot::draw_plot_label(panel_tags[plot_count], 
                                   x = 0.1, y = 0.9, 
                                   size = 2.5,  # Updated panel label size
                                   fontface = "bold")
      
      # Add to list
      plot_list[[plot_count]] <- p
      plot_count <- plot_count + 1
    }
  }
  
  # Arrange and save
  combined <- cowplot::plot_grid(plotlist = plot_list, ncol = 2, align = "hv")
  pdf_file <- file.path(output_dir, "mise_comparison_all_panels_grayscale.pdf")
  cowplot::save_plot(pdf_file, combined, base_width = 6, base_height = 2.5)
  
  message(sprintf("Grayscale combined panel plot saved to %s", pdf_file))
  message("Grayscale LaTeX footnote saved as figure_footnote_grayscale.tex")
  
  return(combined) # Return the plot object
}

# Generate and save the grayscale comparison plots
comparison_plot_grayscale <- plot_comparison_grayscale(comparison_results, output_dir = here::here("output", "plots", "estAMSE_plots"))

# Create a version of the plot using base R graphics to use lty and lwd directly
plot_comparison_base_r <- function(results = comparison_results,
                                  sample_sizes = c(1000),
                                  noise_levels = c(0.1, 0.5),
                                  output_dir = here::here("output", "plots", "estAMSE_plots")) {
  # Setup output directory
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Setup PDF device
  pdf_file <- file.path(output_dir, "mise_comparison_base_r.pdf")
  pdf(pdf_file, width = 7, height = 3)
  
  # Setup layout for two plots side by side
  layout(matrix(1:2, ncol = 2, byrow = TRUE))
  # Adjust margins to make space for panel labels on the left
  # mar = c(bottom, left, top, right), oma = c(bottom, left, top, right)
  par(mar = c(4, 5.5, 2.5, 1)) # Increase left margin for y-axis label
  
  # Create LaTeX footnote
  footnote_text <- paste(
    "The figure compares MISE curves across different regularization parameters $\\lambda$.",
    "\\textbf{Thin gray lines} show estimated MISE curves.",
    "\\textbf{Thick solid line} represents the true MISE curve.",
    "\\textbf{Thin solid line}: true optimal $\\lambda$. \\textbf{Thick solid line}: GCV-selected $\\lambda$. \\textbf{Thick dashed line}: Iterated P-spline $\\lambda$. \\textbf{Circle}: GCV result ($\\lambda$, MISE).",
    "All quantities are plotted on logarithmic scales (base 10)."
  )
  writeLines(footnote_text, file.path(output_dir, "figure_footnote_base_r.tex"))
  
  # Panel labels and counter
  panel_labels <- c("(a)", "(b)")
  plot_count <- 1
  
  # First pass to collect y-axis limits for consistent scaling
  all_y_values <- c()
  
  for (n in sample_sizes) {
    for (sigma in noise_levels) {
      result_key <- paste0("n", n, "_sig", sigma)
      gcv_key <- paste0("gcv_n", n, "_sig", sigma)
      
      if (!result_key %in% names(results) || !gcv_key %in% names(results)) next
      
      iterated_data <- results[[result_key]]$resdf2
      gcv_data <- results[[gcv_key]]
      
      if (is.null(iterated_data) || nrow(iterated_data) == 0 || 
          is.null(gcv_data) || nrow(gcv_data) == 0) next
      
      # Collect y values for scaling
      all_y_values <- c(all_y_values, log10(iterated_data$mise))
      all_y_values <- c(all_y_values, log10(gcv_data$mise_gcv))
    }
  }
  
  # Calculate common y-axis limits
  if (length(all_y_values) > 0) {
    y_range <- range(all_y_values, na.rm = TRUE)
    y_margin <- diff(y_range) * 0.05  # 5% margin
    y_limits <- c(y_range[1] - y_margin, y_range[2] + y_margin)
  } else {
    y_limits <- NULL
  }
  
  # Main plotting loop
  for (n in sample_sizes) {
    for (sigma in noise_levels) {
      # Extract data
      result_key <- paste0("n", n, "_sig", sigma)
      gcv_key <- paste0("gcv_n", n, "_sig", sigma)
      
      if (!result_key %in% names(results) || !gcv_key %in% names(results)) next
      
      iterated_data <- results[[result_key]]$resdf2
      gcv_data <- results[[gcv_key]]
      
      if (is.null(iterated_data) || nrow(iterated_data) == 0 || 
          is.null(gcv_data) || nrow(gcv_data) == 0) next
      
      # Get true MISE data
      max_nsim <- max(iterated_data$nsim, na.rm = TRUE)
      if (is.na(max_nsim)) next
      
      true_mise_data <- iterated_data[iterated_data$nsim == max_nsim, ]
      if (nrow(true_mise_data) == 0) next
      
      min_mise_idx <- which.min(true_mise_data$mise)
      if (length(min_mise_idx) == 0) next
      
      true_mise <- true_mise_data[min_mise_idx, ]
      
      # Calculate mean values
      mean_gcv_lambda <- mean(gcv_data$lambda_gcv, na.rm = TRUE)
      mean_gcv_mise <- mean(gcv_data$mise_gcv, na.rm = TRUE)
      if (is.na(mean_gcv_lambda) || is.na(mean_gcv_mise)) next
      
      # Get iterated lambda
      iterated_min_lambda <- iterated_data %>%
        filter(nsim != max(nsim)) %>%
        group_by(nsim) %>%
        slice_min(mise, with_ties = FALSE) %>%
        summarize(min_lambda = lambda)
      
      mean_iterated_lambda <- mean(iterated_min_lambda$min_lambda, na.rm = TRUE)
      
      # Extract x and y values for plotting
      x_values <- log10(true_mise_data$lambda)
      y_values <- log10(true_mise_data$mise)
      
      # Create the plot with truly default settings and no default box
      plot(x_values, y_values, 
           type = "n", 
           xlim = range(x_values), 
           ylim = y_limits,
           xlab = expression(log[10](lambda)),
           ylab = ifelse(plot_count == 1, expression(log[10](EMSE)),""),
           cex.lab = 1.2, # Increase label size
           main = substitute(paste("n = ", n_val, ", ", sigma, " = ", sig_val), 
                            list(n_val = n, sig_val = sigma)),
           bty = "l")  # No default box
      
      # Add simulation curves (thin gray lines)
      nsim_count <- 0
      for (sim in unique(iterated_data$nsim)) {
        if (sim == max_nsim) next
        nsim_count <- nsim_count + 1
        sim_data <- iterated_data[iterated_data$nsim == sim, ]
        lines(log10(sim_data$lambda), log10(sim_data$mise), 
              col = "gray70", lwd = 2)
      }
      
      # Add true MISE curve (thick black line)
      lines(log10(true_mise_data$lambda), log10(true_mise_data$mise), 
            col = "black", lwd = 3)
      
      # Get axis ranges for later use
      axis_range_x <- par("usr")[1:2]
      axis_range_y <- par("usr")[3:4]
      # # Draw L-shaped box (left and bottom axes only)
      # lines(c(axis_range_x[1], axis_range_x[2]), c(axis_range_y[1], axis_range_y[1]), lwd = 1)  # Bottom
      # lines(c(axis_range_x[1], axis_range_x[1]), c(axis_range_y[1], axis_range_y[2]), lwd = 1)  # Left
      
      # Add vertical lines with specified lty and lwd, but clipped at the max y-value
      max_y <- max(log10(iterated_data$mise), na.rm = TRUE)
      segments(x0 = log10(true_mise$lambda), y0 = axis_range_y[1], 
               x1 = log10(true_mise$lambda), y1 = max_y, 
               lty = 1, lwd = 1, col = "black")      # thick solid
      segments(x0 = log10(mean_gcv_lambda), y0 = axis_range_y[1], 
               x1 = log10(mean_gcv_lambda), y1 = max_y, 
               lty = 1, lwd = 3, col = "gray30")      # thick dashed
      segments(x0 = log10(mean_iterated_lambda), y0 = axis_range_y[1], 
               x1 = log10(mean_iterated_lambda), y1 = max_y, 
               lty = 4, lwd = 3, col = "black")      # thin solid
      
      # # Add GCV result point
      # points(log10(mean_gcv_lambda), log10(mean_gcv_mise), 
      #        pch = 21, col = "black", bg = "white", cex = 1.2, lwd = 1)
      
      # # # Add panel label in the outer left margin, positioned at top of each plot
      # plot_center_y <- if(plot_count == 1) 0.75 else 0.25  # Top of each plot
      # mtext(panel_labels[plot_count], 
      #       side = 2,      # 2 = left side
      #       outer = TRUE,  # Use outer margin
      #       line = 0.5,    # Closer to plot edge in outer margin
      #       at = plot_center_y,  # Position at top of each plot
      #       adj = 0.5,     # Center horizontally
      #       cex = 1.5,  
      #       font = 2)      # bold

      # Increment plot counter
      plot_count <- plot_count + 1
    }
  }
  
  # Close PDF device
  dev.off()
  
  # Return the file path
  message(sprintf("Base R plot saved to %s", pdf_file))
  message("LaTeX footnote saved as figure_footnote_base_r.tex")
  
  # Caption for the plot (commented out for future reference)
  # The caption below could be used in publications or presentations
  # -----------------------------------------------------------------
  # Figure X: Comparison of MISE curves for different regularization parameters.
  # The plot shows the logarithm (base 10) of Mean Integrated Squared Error (MISE)
  # against the logarithm (base 10) of the regularization parameter λ for two settings:
  # Panel (a): n=1000, σ=0.1 and Panel (b): n=1000, σ=0.5.
  # Thin gray lines represent estimated MISE curves from different simulation iterations.
  # The thick black line shows the true MISE curve. Vertical lines indicate:
  # - Black thin solid line: true optimal λ value
  # - Gray thick solid line: GCV-selected λ value
  # - Black thick dashed line: Iterated P-spline selected λ value
  # -----------------------------------------------------------------
  
  return(pdf_file)
}

# Generate the base R plot with direct lty and lwd control
base_r_plot <- plot_comparison_base_r(comparison_results)
