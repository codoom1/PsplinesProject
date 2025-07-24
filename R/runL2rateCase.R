


# =============================================
# Script: runL2rateCase.R
# Purpose: Run L2 convergence simulations for selected cases
# Author: Christopher Odoom
# Date: July 14, 2025
# =============================================
source(here::here("R", "L2rates.R"))
cat("\n============================================\n")
cat("Running L2 Convergence Simulation: slowr2\n")
cat("Sample size range: 1000 to 5000\n")
cat("Noise level (sig): 0.1\n")
cat("Number of simulations per sample size: 10\n")
cat("============================================\n\n")


samp <- samp(range_start=1000, range_end=5000, points_per_range=4)
# Sort sample sizes from largest to smallest
samp <- sort(samp, decreasing = TRUE)
cat("Sample sizes used (largest to smallest): ", paste(samp, collapse=", "), "\n")

run_l2_simulations(sample_size = samp, sig = 0.1, nsim = 100, cases = "fastr2")
cat("\nSimulation for 'fastr2' completed. Results saved in output/data/L2rates_data/fastr2.RData\n")

