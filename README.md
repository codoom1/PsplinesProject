# Clean_projectCodes: FNMR Data Analysis and P-spline Derivative Estimation


## Overview
This project is a comprehensive research codebase for the development, evaluation, and application of advanced derivative estimation methods using P-splines. The central contribution is the **resub** (resubstitution) method for estimating derivatives of smooth functions from noisy data. The project is structured to:

- **Develop the resub derivative estimation method:**
  - Introduces and implements the resub method for estimating first and second derivatives using penalized B-splines (P-splines).
  - Provides robust R functions for fitting, tuning, and extracting derivatives from data.

- **Compare with existing methods:**
  - Benchmarks the resub method against other state-of-the-art derivative estimation techniques, including naive and alternative smoothing approaches.
  - Includes scripts for systematic comparison, visualization, and tabulation of results.

- **Theoretical analysis:**
  - Empirically investigates the L2 rates of convergence for the resub and competing methods.
  - Provides simulation studies and plots to illustrate theoretical properties and practical performance.

- **Applications:**
  - Applies the developed methods to real-world datasets, including:
    - **FNMR (Fast Nuclear Magnetic Resonance) spectral data:** For chemical peak detection and derivative analysis across time points.
    - **Stock market data (e.g., IONQ):** For financial time series analysis and feature extraction.

The codebase is organized for reproducible research, supporting both methodological development and practical data analysis. It includes utilities for data preprocessing, method evaluation, and publication-quality figure generation. All scripts are modular and can be run independently or as part of a full workflow.

## Project Structure
```
Clean_projectCodes.Rproj           # RStudio project file
renv.lock                         # renv environment lock file for reproducibility
submit_job.sh                     # Shell script for job submission (if applicable)
data/
  raw/                            # Raw input data (e.g., FNMRData.csv)
  processed/                      # Processed data (e.g., IONQ_stock_data.txt)
logs/                             # Log files
output/
  data/                           # Output data from analyses
  plots/                          # Generated plots (PDF, PNG)
R/
  Application_spectra.R           # Main script for FNMR spectra application
  Applications_stock.R            # Script for stock data application
  compareMethods_plot.R           # Plotting for method comparison
  compareMethods.R                # Method comparison logic
  competingMethods.R              # Competing methods implementation
  EstimatedAMSE.R                 # AMSE estimation
  FNMR.R                          # Core FNMR data processing and derivative estimation
  L2rates.R                       # L2 convergence rate analysis
  L2ratesplots.R                  # L2 rate plotting
  MainFunctions.R                 # Main utility functions
  oldprojectScript.R              # Legacy script
  projectPrelims.R                # Preliminary analysis
  runL2rateCase.R                 # Run L2 rate case studies
renv/                             # renv environment for package management
```

## Getting Started

### 1. Prerequisites
- **R (>= 4.0)**
- **RStudio** (recommended)
- The following R packages (managed via `renv`):
  - `here`, `ggplot2`, `dplyr`, `tidyr`, and others as specified in `renv.lock`

To restore the project environment:
```r
install.packages("renv")
renv::restore()
```

### 2. Data
- Place your raw FNMR data (e.g., `FNMRData.csv`) in `data/raw/`.
- Processed and output data will be saved in `data/processed/` and `output/data/` respectively.

## Main Scripts and Their Usage

### 1. `R/Application_spectra.R`
**Purpose:**
- Loads FNMR data, preprocesses chemical shifts, estimates derivatives for all time points, and generates a variety of plots for analysis.

**How to Run:**
Open in RStudio and run all lines, or in the terminal:
```r
source("R/Application_spectra.R")
```

**Key Steps Performed:**
- Loads and normalizes FNMR data.
- Estimates mean, first, and second derivatives using P-splines.
- Visualizes spectra, derivatives, and peak detection for all time points.
- Compares derivative values at the true peak.
- Generates publication-ready plots (saved in `output/plots/`).

**Notable Functions Used:**
- `load_data()`, `preprocess_shift()`, `process_all_fnmr()`, `plot_fnmr_peaks()`, `plot_fnmr_panels_base()`, `plot_fnmr_deriv_zoom()`, `get_deriv_at_true_peak()`, `plot_fnmr_panels_maxima_detection_interp()`

### 2. `R/Applications_stock.R`
**Purpose:**
- Applies similar P-spline analysis to stock data (e.g., IONQ stock prices).

**How to Run:**
```r
source("R/Applications_stock.R")
```

### 3. `R/compareMethods.R` and `R/compareMethods_plot.R`
**Purpose:**
- Implements and visualizes comparisons between different derivative estimation methods.

**How to Run:**
```r
source("R/compareMethods.R")
source("R/compareMethods_plot.R")
```

### 4. `R/EstimatedAMSE.R`
**Purpose:**
- Estimates the Asymptotic Mean Squared Error (AMSE) for derivative estimators.

**How to Run:**
```r
source("R/EstimatedAMSE.R")
```

### 5. `R/L2rates.R` and `R/L2ratesplots.R`
**Purpose:**
- Analyzes and visualizes L2 convergence rates for the estimators.

**How to Run:**
```r
source("R/L2rates.R")
source("R/L2ratesplots.R")
```

### 6. `R/FNMR.R`
**Purpose:**
- Contains core functions for FNMR data processing, P-spline fitting, and derivative estimation.
- Not intended to be run directly; sourced by other scripts.

### 7. `R/MainFunctions.R`, `R/competingMethods.R`, `R/projectPrelims.R`, `R/oldprojectScript.R`, `R/runL2rateCase.R`
**Purpose:**
- Utility functions, alternative methods, preliminary analyses, and legacy code.
- Run or source as needed for specific analyses.

## Output
- **Plots:** Saved in `output/plots/` (PDF, PNG)
- **Processed Data:** Saved in `output/data/`
- **Logs:** Saved in `logs/`

## Reproducibility
- The project uses `renv` for package management. Always run `renv::restore()` before analysis to ensure package versions match those used in development.

## Troubleshooting
- If you encounter missing package errors, ensure you have run `renv::restore()`.
- For file path issues, ensure you are running scripts from the project root or using RStudio's project functionality.

## Contact
For questions or contributions, please contact the project maintainer or open an issue in your version control system.
