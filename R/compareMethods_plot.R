# ============================================================================
# File: compareMethods_plot.R
# Purpose: Load simulation results and generate plots and tables for analysis
# Author: [Christopher Odoom]
# Date: July 15, 2025
#
# Description:
#   - Loads saved simulation results from compareMethods_data.
#   - Prepares data for plotting and analysis.
#   - Generates MISE boxplots and LaTeX tables.
# ============================================================================

library(here)
library(ggplot2)
library(dplyr)
library(tidyr)



# ---- Load saved simulation results before analysis or plotting ----
load(here::here("output", "data", "compareMethods_data", "emsea1.RData")) # loads sima1
load(here::here("output", "data", "compareMethods_data", "emseb1.RData")) # loads simb1
load(here::here("output", "data", "compareMethods_data", "emsea2.RData")) # loads sima2
load(here::here("output", "data", "compareMethods_data", "emseb2.RData")) # loads simb2

# ---- Prepare data for plotting ----
get_mise_df <- function(res_list, sigs, fun, deriv) {
  do.call(rbind, lapply(seq_along(sigs), function(i) {
    n_sim <- nrow(res_list[[i]])
    do.call(rbind, lapply(seq_len(n_sim), function(j) {
      data.frame(
        sigma = sigs[i],
        sim = j,
        estimator = estimators,
        MISE = sapply(estimators, function(est) res_list[[i]][j, est]),
        FunType = fun,
        derivative = deriv
      )
    }))
  }))
}

sigs <- c(0.1, 0.5, 2)
estimators <- c("naive", "simple", "plug.in", "resub", "oracle", "GCPDAvg", "Dailc")
mise_df <- rbind(
  get_mise_df(sima1, sigs, "a", 1),
  get_mise_df(simb1, sigs, "b", 1),
  get_mise_df(sima2, sigs, "a", 2),
  get_mise_df(simb2, sigs, "b", 2)
)

# ---- Set up parameters and color palette ----

col_vec <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d")

# ---- Plot: All MISE boxplots on the same graph, faceted by function, derivative order, and noise level ----
library(ggplot2)
p <- ggplot(mise_df, aes(x = estimator, y = MISE, fill = estimator)) +
  geom_boxplot(outlier.size = 0.7, lwd = 0.3) +
  facet_wrap(~ derivative + FunType + sigma, labeller = label_both, scales = "free_y") +
  scale_fill_manual(values = col_vec) +
  labs(title = "MISE for All Estimators, Functions, Derivatives, and Noise Levels",
       y = "MISE", x = "Estimator") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")

ggsave(here::here("output", "plots", "compareMethods_plots", "MISE_all_boxplots.png"), p, width = 14, height = 6, dpi = 300)
print(p)




# ---- Helper function to format mean (se) for LaTeX ----
format_mise <- function(x) {
  m <- mean(x)
  s <- sd(x) / sqrt(length(x))
  sprintf("%.3f (%.3f)", m, s)
}

# ---- Function to extract and format results for a given noise level, function, and derivative order ----
extract_results <- function(res_df, est_names) {
  sapply(est_names, function(est) format_mise(res_df[[est]]))
}

# ---- Function to generate LaTeX table for MISE results ----
generate_mise_latex_table <- function(sigs, sima1, simb1, sima2, simb2, estimators) {
  cat(
    "\\begin{table}[htbp]\n",
    "\\centering\n",
    "\\begin{threeparttable}\n",
    "\\caption{Finite-sample MISE (mean and standard error) for different estimators under varying noise levels ($n = 500$).}\n",
    "\\label{tab:finiteP300}\n",
    "\\begin{tabular}{llcccc}\n",
    "\\toprule\n",
    " &  & \\multicolumn{2}{c}{First Derivative} & \\multicolumn{2}{c}{Second Derivative} \\\n",
    "\\$\\sigma$ & Method & $f_a$ & $f_b$ & $f_a$ & $f_b$ \\\n",
    "\\midrule\n",
    paste(
      unlist(
        lapply(seq_along(sigs), function(i) {
          sigma <- sigs[i]
          fa1 <- extract_results(sima1[[i]], estimators)
          fb1 <- extract_results(simb1[[i]], estimators)
          fa2 <- extract_results(sima2[[i]], estimators)
          fb2 <- extract_results(simb2[[i]], estimators)
          c(
            sapply(seq_along(estimators), function(j) {
              method <- estimators[j]
              sprintf("%s & %s & %s & %s & %s & %s \\\\",
                ifelse(j == 1, sprintf("\\\\multirow{7}{*}{%.1f}", sigma), ""),
                switch(method,
                  naive = "Naive",
                  simple = "Simple",
                  plug.in = "Plug-in",
                  resub = "Resub",
                  oracle = "Oracle",
                  GCPDAvg = "GCPDAvg",
                  Dailc = "DaiM"
                ),
                fa1[j], fb1[j], fa2[j], fb2[j]
              )
            }),
            "\\midrule"
          )
        })
      ),
      collapse = "\n"
    ),
    "\\bottomrule\n",
    "\\end{tabular}\n",
    "\\begin{tablenotes}\n",
    "\\small\n",
    "\\item Note: Values are mean integrated squared errors (MISE) with Monte Carlo standard errors in parentheses, based on 500 simulations. Bolded values indicate the best performance (excluding oracle).\n",
    "\\end{tablenotes}\n",
    "\\end{threeparttable}\n",
    "\\end{table}\n"
  )
}

# ---- Call the LaTeX table generation function ----
generate_mise_latex_table(sigs, sima1, simb1, sima2, simb2, estimators)


