# # Install necessary pacnsegages if not already installed
# if (!require(ggplot2)) install.pacnsegages("ggplot2")
# if (!require(gridExtra)) install.pacnsegages("gridExtra")

library(here)

library(ggplot2)
library(gridExtra)
library(dplyr)
source("R/mainFunctions.R") # Load the functions from the provided file
######## Investigating the L2 rates of convergence ###########

# Define the function
choose_grid <- function(x.grid1, xtype, percentage = 0.1) {
  n <- length(x.grid1)
  # Calculate the number of points based on the given percentage
  points <- floor(n * percentage)
  
  if (xtype == "int") {
    start_index <- ceiling((n - points) / 2)
    end_index <- start_index + points - 1
  } else {
    start_index <- n - points
    end_index <- n - 1
  }
  
  x.grid <- x.grid1[start_index:end_index]
  return(x.grid)
}

# This script contains functions and simulations to investigate the L2 rates of convergence
# and compare different derivative estimation methods using penalized splines.
# The main function, `investigate_L2_convergence`, performs simulations to estimate the Mean Integrated Squared Error (MISE)
# for various methods (naive, plug-in, re-substitution, and oracle) and plots the results.

# Rename the main function to reflect the purpose of the file
investigate_L2_convergence <- function(sample_sizes = seq(100, 1000, 100), sig = 0.001, r = 1, 
                                       pord = 2, bdeg = 3, nsim = 1,
                                       type = "slow") {
  results <- data.frame()
  slopes <- data.frame(method = character(), slope = numeric(), stringsAsFactors = FALSE)

  for (n in sample_sizes) {
    x <- seq(from = 0, to = 1, length = n)
    f <- 32 * exp(-8 * (1 - 2 * x)^2) * (1 - 2 * x)
    x.grid <- x
    if (r == 0) {
      fr.grid <- 32 * exp(-8 * (1 - 2 * x.grid)^2) * (1 - 2 * x.grid)
    } else if (r == 1) {
      fr.grid <- (4096 * x.grid^2 - 4096 * x.grid + 960) * exp(-8 * (1 - 2 * x.grid)^2)
    } else {
      fr.grid <- -(262144 * x.grid^3 - 393216 * x.grid^2 + 184320 * x.grid - 26624) * exp(-8 * (1 - 2 * x.grid)^2)
    }

    m <- pord
    v <- (4 * m + 1)

    if (type == "slow") {
      nseg<- ceiling(25 * n^(1 / v))
    } else {
      nseg<- ceiling(25 * n^(1 / v + 0.1))
    }
    keep.emise <- data.frame(
      naive = rep(NA, nsim),
      plug.in = rep(NA, nsim),
      resub = rep(NA, nsim),
      oracle = rep(NA, nsim)
    )

    for (i in 1:nsim) {
      # Generate the data
      start_time <- Sys.time()  # Start time for the simulation
      y <- f + sig * rnorm(n)
      naive.fit <- naive.est.opt(x = x, y = y, r = r, x.grid = x.grid, nseg = nseg, pord = pord, bdeg = bdeg) # nolint

      plugin.fit <- plugin.est(x=x, y=y, r=r,nseg=nseg , pord = pord ,
                        bdeg=bdeg, x.grid=x.grid,fr="Y") # nolint

      resub.fit <- resub.est(x = x, y = y, r = r, x.grid = x.grid, nseg= nseg, pord = pord,bdeg = bdeg, tol = 1e-5, ITs = 100) # nolint

      oracle.fit <- oracle.est(initial.lambda = resub.fit$lambda, x = x, y = y, r = r, fr.grid = fr.grid, nseg= nseg, pord = pord, bdeg = bdeg, x.grid = x.grid) # nolint

      keep.emise$naive[i] <- sqrt(mean((fr.grid - naive.fit$fr.hat)^2))
      keep.emise$plug.in[i] <- sqrt(mean((fr.grid - plugin.fit$fr.hat)^2))
      keep.emise$resub[i] <- sqrt(mean((fr.grid - resub.fit$fr.hat)^2))
      keep.emise$oracle[i] <- sqrt(mean((fr.grid - oracle.fit$fr.hat)^2))
      print("*******")
      print(keep.emise[i,])
      print("*******")
      print(paste("Simulation", i, "completed for sample size", n))
      print("*******")

        end_time <- Sys.time()  # End time for the simulation
        elapsed_time <- end_time - start_time  # Calculate elapsed time
        print(paste("Elapsed time for sample size", n, ":", elapsed_time))
        print("*******")
    }

    # Ensure that `keep.emise` is not empty before applying `na.omit`
    if (!is.null(keep.emise) && nrow(keep.emise) > 0) {
        keep.emise <- na.omit(keep.emise)  # Remove rows with NA values
        emise <- apply(keep.emise, 2, mean)  # Calculate the mean for each column
        log10emise <- log10(emise)  # Compute the log10 of the means
        print("*******")
        print(emise)
        print(n)
    } else {
        warning("keep.emise is empty or NULL. Skipping calculations.")
    }

    temp_results <- data.frame(
      method = names(emise),
      log10n = log10(n),
      log10emise = log10emise
    )

    results <- rbind(results, temp_results)

    for (m in unique(temp_results$method)) {
      method_results <- results[results$method == m, ]
      if (nrow(method_results) > 1) {
        lm_fit <- lm(log10emise ~ log10n, data = method_results)
        slope <- coef(lm_fit)[2]
        print(paste("Slope for method", m, ":", slope))

        slopes <- rbind(slopes, data.frame(method = m, slope = slope, stringsAsFactors = FALSE))
      }
    }

    # # Plotting
    # p <- results %>% ggplot(aes(x = log10n, y = log10emise, color = method)) +
    #   geom_point() +
    #   geom_smooth(method = "lm", se = FALSE) +
    #   labs(title = "Log10(n) vs Log10(emise)", x = "log10(n)", y = "log10(emise)") +
    #   theme_minimal() +
    #   theme(
    #     panel.background = element_rect(fill = "gray95"),
    #     plot.background = element_rect(fill = "gray90"),
    #     panel.grid.major = element_line(color = "white"),
    #     panel.grid.minor = element_line(color = "white")
    #   )

    # # Ensure the plot is displayed by explicitly calling print() and wrapping the ggplot object in a print statement
    # print(p)  # Display the plot

    # # Optionally, save the plot to a file for later review
    # ggsave(filename = here::here("output", "plots","L2rates_plots", paste0("plot_L2_convergence_n_r", n,r, ".png")), plot = p, width = 8, height = 6)

    # Sys.sleep(0.0005)
    
  }

  return(results)
}

# the first 
## 0.67
## 0.44

samp <- function(range_start=100,range_end=1000,points_per_range=2 ){
  sample_size <- c()
  # Loop over each 1000-unit range and generate points
  for (i in seq(range_start, range_end - range_start, by = range_start)) {
    log_start <- log10(i)
    log_end <- log10(i + range_start)
    logn <- seq(log_start, log_end, length.out = points_per_range + 1)
    # Convert the log scale sequence to the normal scale and append to sample_size
    sample_size <- c(sample_size, round(10^(logn), 0))
  }
  # Ensure unique and sorted values
  sample_size <- sort(unique(sample_size))
  return(sample_size)
}


# ---- Helper function to run and save L2 convergence simulations for different settings ----
run_l2_simulations <- function(sample_size, sig, nsim = 10, cases = c("slowr1", "fastr1", "slowr2", "fastr2")) {
  # Define r, pord, and bdeg for each type
  params <- list(
    slowr1 = list(r = 1, pord = 2, bdeg = 4, type = "slow"),
    fastr1 = list(r = 1, pord = 2, bdeg = 4, type = "fast"),
    slowr2 = list(r = 2, pord = 2, bdeg = 4, type = "slow"),
    fastr2 = list(r = 2, pord = 2, bdeg = 4, type = "fast")
  )
  for (case in cases) {
    par <- params[[case]]
    res <- investigate_L2_convergence(
      sample_sizes = sample_size,
      sig = sig,
      r = par$r,
      pord = par$pord,
      bdeg = par$bdeg,
      nsim = nsim,
      type = par$type
    )
    save(res, file = here::here("output", "data", "L2rates_data", paste0(case, ".RData")))
    print(paste("Completed:", case))
  }
  cat("All selected L2 convergence simulations completed and saved.\n")
}

# Example call to run all L2 simulations in one line
#samp1 <- samp(range_start=500,range_end=3000,points_per_range=1)
#run_l2_simulations(sample_size = samp1, sig = 0.1, nsim = 10, cases = "slowr1")


