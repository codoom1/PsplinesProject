# File: Application.R
# Project: Application of the improved derivative estimator
# Created Date: 2024-09-01
## This script is design to show application of the improved derivative
## estimator. The data used in this example is the stock price of IONQ
## (Quantum computing company) from Yahoo Finance. The data is used to
## estimate the derivative of the stock price using the proposed method
## and compare it with the naive method. 

## The script is divided into the following sections:
## 1. Load the required libraries
## 2. Pull IONQ stock data from Yahoo Finance
## 3. Prepare the data
## 4. Data Preprocessing
## 5. Estimate the derivative using the proposed method
## 6. Plot the results
## 7. Save the plot for your paper

## 1. Load the required libraries
## Note: You may need to install the required packages if you haven't already
## Note: The libraries used in this script are:
library(here)
library(JOPS)
library(dplyr)
library(ggplot2)
library(MASS)
library(lattice)
library(gridExtra)
#install.packages(c("quantmod", "patchwork"))
library(quantmod)
library(patchwork)
library(scales)  # for alpha transparency

### use functions from the main functions file
# Note: The MainFunctions.R file contains the functions used in this script.
source("R/MainFunctions.R")



######### Application of the proposed method to  Quantum computing stock price #######################################
# 1. Try to fetch new IONQ stock data from Yahoo Finance, fallback to existing data
# First try to fetch from Yahoo Finance, if it fails use existing data
if (!file.exists(here::here("data", "processed", "IONQ_stock_data.txt"))) {
  cat("Fetching new IONQ stock data from Yahoo Finance...\n")
  # Get the stock data for IONQ from Yahoo Finance
  getSymbols("IONQ", src = "yahoo", from = "2024-09-01", to = "2025-07-01")

  # 2. Prepare the data
  # Convert the IONQ data to a data frame and select relevant columns
  stock_data <- IONQ %>%
    as.data.frame() %>%
    mutate(
      time = as.Date(rownames(.)),
      price = IONQ.Adjusted
    ) %>%
    dplyr::select(time, price) %>%
    arrange(time)

    ## save the stock data to a text file
write.table(stock_data, file = here::here("data", "processed", "IONQ_stock_data.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  cat("IONQ stock data fetched and saved to 'data/processed/IONQ_stock_data.txt'.\n")
  
  # Check if the data was fetched successfully
  if (nrow(stock_data) == 0) {
    stop("No data fetched from Yahoo Finance. Please check the stock symbol or date range.")
  }
} else {
  cat("Loading existing IONQ stock data...\n")
  stock_data <- read.table(here::here("data", "processed", "IONQ_stock_data.txt"), 
                          sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  stock_data$time <- as.Date(stock_data$time)
  # Select only the basic columns for reprocessing
  stock_data <- stock_data %>%
    dplyr::select(time, price) %>%
    arrange(time)
}


# Check the structure of the data
str(stock_data)
# Check the first few rows of the data
head(stock_data)
# Check the last few rows of the data
tail(stock_data)
# Check for missing values
missing_values <- sum(is.na(stock_data))
if (missing_values > 0) {
  cat("There are", missing_values, "missing values in the data.\n")
} else {
  cat("No missing values in the data.\n")
}
# Check the range of dates
date_range <- range(stock_data$time)
cat("Date range:", date_range[1], "to", date_range[2], "\n")
# Check the range of prices
price_range <- range(stock_data$price)
cat("Price range:", price_range[1], "to", price_range[2], "\n")
# Check the number of rows
num_rows <- nrow(stock_data)
cat("Number of rows:", num_rows, "\n")
# Check the number of columns
num_cols <- ncol(stock_data)
cat("Number of columns:", num_cols, "\n")
# Check the column names
col_names <- colnames(stock_data)
cat("Column names:", paste(col_names, collapse = ", "), "\n")
# Check the summary of the data
summary(stock_data)
# Check the distribution of prices
png(here::here("output", "plots", "Appl_plots","IONQ_stock_prices_distribution.png"), width = 800, height = 600)
hist(stock_data$price, main = "Distribution of IONQ Stock Prices", xlab = "Price ($)", col = "lightblue", border = "black")

dev.off()

## 3. Data Preprocessing
# Normalize time to [0, 1] range
# Note: This is important for the derivative estimation.
time_numeric <- as.numeric(stock_data$time)
time_range <- max(time_numeric) - min(time_numeric)
time_normalized <- (time_numeric - min(time_numeric)) / time_range
# check common difference between normalized time
diff_time <- diff(time_normalized)
stock_data$time_normalized <- time_normalized
## find the normalizing constant to rescale the derivative
normalizing_constant <- time_range
## Note: The normalizing constant is the range of the time variable



##### Estimating derivative with our method resub and comparing with naive method  ###
x<- stock_data$time_normalized
y <- stock_data$price
nseg <- 25 ## number of segments, this is the number of knots-1
bdeg <- 3 ## degree of the spline
pord <- 1 ## order of the penalty
r <- 1 ## the order of the derivative we want to estimate

# Estimate the derivative using the naive derivative estimator
naive.est <- naive.est.opt(x=x, y=y,r=r
                       , nseg = nseg, bdeg = bdeg,pord = pord,
                       x.grid = x)

# Estimate the derivative using the proposed method
# Note: The resub.est function is a wrapper for the resub function
resub.fr<- resub.est(x=x, y=y, r=r, x.grid = x,nseg=nseg, pord = pord,
                     bdeg=bdeg,tol=1e-20, ITs=100)

                     
 ## Check the smoothing parameter
# Check the smoothing parameter for the naive estimator
naive.est$lambda
# Check the smoothing parameter for the resub(iterative) estimator
resub.fr$lambda

## Add the estimated derivatives to the data frame
# Note: The naive estimator is already in the correct format
# Note: The resub estimator is a list, so we need to extract the fitted values
stock_data <- stock_data %>%
  mutate(
    price_fit = naive.est$f.hat, ##fitted mean function
    dt_days = as.numeric(difftime(time, lag(time), units = "days")),  # Time diff in days
    derivative_emp = c(NA, diff(price) / dt_days[-1]),  # The empirical derivative estimate using finite difference Already in $/day
    derivative_resub = (resub.fr$fr.hat/normalizing_constant), ### The resub derivative estimator
    derivative_naive = (naive.est$fr.hat/normalizing_constant) ### The naive derivative estimator
  )

plot(stock_data$derivative_naive, stock_data$derivative_emp, xlab = "Naive Derivative", ylab = "Empirical Derivative", main = "Naive vs Empirical Derivative")

# Save the IONQ stock data to a text file
write.table(stock_data, file = here::here("data", "processed", "IONQ_stock_data.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

# --- Print a summary table comparing the derivatives ---
deriv_summary <- data.frame(
  Estimator = c("Empirical", "Resubstitution", "Naive"),
  Mean = c(mean(stock_data$derivative_emp, na.rm = TRUE),
           mean(stock_data$derivative_resub, na.rm = TRUE),
           mean(stock_data$derivative_naive, na.rm = TRUE)),
  SD = c(sd(stock_data$derivative_emp, na.rm = TRUE),
         sd(stock_data$derivative_resub, na.rm = TRUE),
         sd(stock_data$derivative_naive, na.rm = TRUE)),
  Min = c(min(stock_data$derivative_emp, na.rm = TRUE),
          min(stock_data$derivative_resub, na.rm = TRUE),
          min(stock_data$derivative_naive, na.rm = TRUE)),
  Max = c(max(stock_data$derivative_emp, na.rm = TRUE),
          max(stock_data$derivative_resub, na.rm = TRUE),
          max(stock_data$derivative_naive, na.rm = TRUE))
)

print(deriv_summary)

# 4. Plot the results
# Set date range limits
date_limits <- as.Date(c("2024-09-01", Sys.Date()))

# --- Plot 1: Price + Smooth Trend ---
p1 <- ggplot(stock_data, aes(x = time)) +
  geom_point(aes(y = price), color = "gray50", alpha = 0.6, size = 1.5) +
  geom_line(aes(y = price_fit), color = "#3366FF", linewidth = 0.8) +
  scale_x_date(
    limits = date_limits,
    date_labels = "%b %Y",  # "Sep 2024" format
    date_breaks = "1 months"  # Tick every 1 months
  ) +
  labs(title = "IONQ Stock Price with P-spline fit", y = "Price ($)", x = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(size = 12, face = "bold"))

# --- Plot 2: Empirical Derivative ---
p2 <- ggplot(stock_data[-1, ], aes(x = time)) +
  geom_line(aes(y = derivative_emp), color = "#FF4E3F", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  scale_x_date(
    limits = date_limits,
    date_labels = "%b %Y",
    date_breaks = "1 months"
  ) +
  labs(title = "Empirical (Noisy) Derivative", y = "Derivative ($/day)", x = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(size = 12, face = "bold"))

# --- Plot 3: Smoothed Derivatives ---
p3 <- ggplot(stock_data[-1, ], aes(x = time)) +
  geom_line(aes(y = derivative_resub, color = "Resub"), linewidth = 0.6) +
  geom_line(aes(y = derivative_naive, color = "Naive"), linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  scale_x_date(
    limits = date_limits,
    date_labels = "%b %Y",
    date_breaks = "1 months"
  ) +
  scale_color_manual(
    name = NULL,
    values = c("Resub" = "#48B023", "Naive" = "deeppink"),
    guide = guide_legend(
      direction = "vertical",
      keyheight = unit(0.5, "cm"),
      override.aes = list(linewidth = 2)
    )
  ) +
  labs(
    title = "Smoothed Derivative Comparison",
    y = "Derivative ($/day)",
    x = "Date"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position.inside = c(1, 0.35),
    legend.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 12, face = "bold")
  )

# --- Combine Plots ---
combined_plots <- (p1 / p2 / p3) + 
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(
    plot.margin = margin(10, 15, 10, 15),
    plot.tag = element_text(face = "bold")
  )

print(combined_plots)
# Save the plot for your paper
ggsave(here::here("output", "plots","Appl_plots", "ionq_stock_analysis.png"), combined_plots, width = 8, height = 5, dpi = 300)
ggsave(here::here("output", "plots","Appl_plots", "ionq_stock_analysis.pdf"), combined_plots, width = 8, height = 5, device = cairo_pdf)

# --- Compare Distributions of Derivative Estimates ---

# Prepare data for histogram plot
hist_data <- stock_data[-1, ] %>%
  dplyr::select(derivative_emp, derivative_resub, derivative_naive) %>%
  tidyr::pivot_longer(cols = everything(), names_to = "method", values_to = "value")

# Plot overlaid histograms
hist_plot <- ggplot(hist_data, aes(x = value, fill = method)) +
  geom_histogram(alpha = 0.4, position = "identity", bins = 40) +
  scale_fill_manual(
    values = c(
      "derivative_emp" = "#FF4E3F",
      "derivative_resub" = "#48B023",
      "derivative_naive" = "deeppink"
    ),
    labels = c(
      "derivative_emp" = "Empirical",
      "derivative_resub" = "Resub",
      "derivative_naive" = "Naive"
    ),
    name = "Estimator"
  ) +
  labs(
    title = "Distribution of Derivative Estimates",
    x = "Derivative ($/day)",
    y = "Count"
  ) +
  theme_minimal(base_size = 12)

print(hist_plot)
ggsave(here::here("output", "plots","Appl_plots", "IONQ_derivative_histograms.png"), hist_plot, width = 7, height = 4, dpi = 300)

# --- Compare Distributions of Derivative Estimates with base R hist() ---

png(here::here("output", "plots","Appl_plots", "IONQ_derivative_histograms_baseR.png"), width = 800, height = 500)

# Remove NAs for fair comparison
d_emp <- na.omit(stock_data$derivative_emp)
d_resub <- na.omit(stock_data$derivative_resub)
d_naive <- na.omit(stock_data$derivative_naive)

# Find common range for all histograms
all_vals <- c(d_emp, d_resub, d_naive)
xlim <- range(all_vals, na.rm = TRUE)

# Plot the empirical derivative histogram first
hist(d_emp, breaks = 40, col = rgb(1, 0.3, 0.25, 0.4), border = "red", xlim = xlim,
     main = "Distribution of Derivative Estimates", xlab = "Derivative ($/day)", ylab = "Count")
# Overlay the resubstitution estimator
hist(d_resub, breaks = 40, col = rgb(0.2, 0.7, 0.2, 0.4), border = "darkgreen", add = TRUE)
# Overlay the naive estimator
hist(d_naive, breaks = 40, col = rgb(0.2, 0.3, 1, 0.4), border = "blue", add = TRUE)

legend("topright", legend = c("Empirical", "Resubstitution", "Naive"),
       fill = c(rgb(1, 0.3, 0.25, 0.4), rgb(0.2, 0.7, 0.2, 0.4), rgb(0.2, 0.3, 1, 0.4)),
       border = c("red", "darkgreen", "blue"),
       bty = "n", cex = 1.1)

dev.off()

# =============================================================================
# NEW GRAYSCALE PLOTS FOR BIOMETRIKA SUBMISSION
# =============================================================================
# --- Biometrika version: Grayscale with line types, 1 row 3 columns ---

# Set date range limits
date_limits <- as.Date(c("2024-09-01", "2025-07-01"))

# --- Plot 1: Price + Smooth Trend (Grayscale) ---
p1_bw <- ggplot(stock_data, aes(x = time)) +
  geom_point(aes(y = price), color = "gray50", alpha = 0.6, size = 1.5) +
  geom_line(aes(y = price_fit), color = "black", linewidth = 0.8) +
  scale_x_date(
    limits = date_limits,
    date_labels = "%b %Y",
    date_breaks = "2 months",
    expand = c(0.02, 0)
  ) +
  labs(title = "IONQ Stock Price with P-spline fit", y = "Price ($)", x = "Date") +
  theme_classic(base_size = 11) +  # Use theme_classic for L-shaped box
  theme(
    plot.title = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    panel.grid = element_blank(),  # Remove all grid lines
    axis.line.y.right = element_blank(),  # Remove right axis line
    axis.line.x.top = element_blank(),    # Remove top axis line
    axis.ticks.length = unit(0.2, "cm"),  # Adjust tick length
    axis.line = element_line(linewidth = 0.5)  # Make axis lines slightly thicker
  )

# --- Plot 2: Empirical Derivative (Grayscale) ---
p2_bw <- ggplot(stock_data[-1, ], aes(x = time)) +
  geom_line(aes(y = derivative_emp), color = "black", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_date(
    limits = date_limits,
    date_labels = "%b %Y",
    date_breaks = "2 months",
    expand = c(0.02, 0)
  ) +
  labs(title = "Empirical (Noisy) Derivative", y = "Derivative ($/day)", x = "Date") +
  theme_classic(base_size = 11) +  # Use theme_classic for L-shaped box
  theme(
    plot.title = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    panel.grid = element_blank(),  # Remove all grid lines
    axis.line.y.right = element_blank(),  # Remove right axis line
    axis.line.x.top = element_blank(),    # Remove top axis line
    axis.ticks.length = unit(0.2, "cm"),  # Adjust tick length
    axis.line = element_line(linewidth = 0.5)  # Make axis lines slightly thicker
  )

# --- Plot 3: Smoothed Derivatives (Grayscale with line types) ---
p3_bw <- ggplot(stock_data[-1, ], aes(x = time)) +
  geom_line(aes(y = derivative_resub, linetype = "Resub"), color = "black", linewidth = 0.7) +
  geom_line(aes(y = derivative_naive, linetype = "Naive"), color = "black", linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_date(
    limits = date_limits,
    date_labels = "%b %Y",
    date_breaks = "2 months",
    expand = c(0.02, 0)
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c("Resub" = "solid", "Naive" = "dotted"),
    guide = guide_legend(
      direction = "vertical",
      keyheight = unit(0.5, "cm"),
      override.aes = list(linewidth = 1.2)
    )
  ) +
  labs(
    title = "Smoothed Derivative Comparison",
    y = "Derivative ($/day)",
    x = "Date"
  ) +
  theme_classic(base_size = 11) +  # Use theme_classic for L-shaped box
  theme(
    legend.position.inside = c(0.85, 0.8),
    legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.3),
    legend.margin = margin(4, 6, 4, 6),
    plot.title = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    panel.grid = element_blank(),  # Remove all grid lines
    axis.line.y.right = element_blank(),  # Remove right axis line
    axis.line.x.top = element_blank(),    # Remove top axis line
    axis.ticks.length = unit(0.2, "cm"),  # Adjust tick length
    axis.line = element_line(linewidth = 0.5)  # Make axis lines slightly thicker
  )

# --- Combine Plots in 1 row, 3 columns ---
combined_plots_bw <- (p1_bw | p2_bw | p3_bw) + 
  plot_layout(widths = c(1, 1, 1)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(
    plot.margin = margin(5, 10, 5, 10),
    plot.tag = element_text(face = "bold", size = 12)
  )

print(combined_plots_bw)

# Save the grayscale plot for Biometrika submission
ggsave(here::here("output", "plots","Appl_plots", "ionq_stock_analysis_biometrika.png"), 
       combined_plots_bw, width = 12, height = 4, dpi = 300)
ggsave(here::here("output", "plots","Appl_plots", "ionq_stock_analysis_biometrika.pdf"), 
       combined_plots_bw, width = 12, height = 4, device = cairo_pdf)

# Caption for the Biometrika plot:
# Figure 1. Analysis of IONQ stock price data using P-spline derivative estimation. 
# (a) Daily adjusted closing prices (gray points) with P-spline smooth fit (solid black line).
# (b) Empirical derivative estimates computed using finite differences (solid black line).
# (c) Comparison of smoothed derivative estimates: resubstitution method (solid line) 
# versus naive method (dotted line). The resubstitution method produces less variable 
# derivative estimates while preserving the underlying trend structure.

# =============================================================================
# END OF BIOMETRIKA PLOTS
# =============================================================================
