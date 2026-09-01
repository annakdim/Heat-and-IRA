################################################################################
################################################################################
####                            PROJECT: Heat-and-IRA                       ####
####            DEMO: End-to-end run on simulated data                     ####
####                                                                        ####
####  This script demonstrates the core modelling approach used in         ####
####  "3_Analysis case-crossover/3a_Analysis_ERA5_percentile_ths.R" on the ####
####  small synthetic dataset created by demo/simulate_demo_data.R.        ####
####  It is a minimal, fast-running illustration of the pipeline, not a    ####
####  reproduction of the manuscript's results (which require the         ####
####  restricted clinical + ERA5 data described in the paper).             ####
####                                                                        ####
####  Expected run time on a normal desktop computer: < 1 minute           ####
################################################################################
################################################################################

rm(list = ls())

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(survival)   # clogit()

# ---------------------------------------------------------------------------
# 0. Simulate the demo data if it doesn't exist yet
# ---------------------------------------------------------------------------
if (!file.exists("demo/demo_data.RData")) {
  source("demo/simulate_demo_data.R")
}
load("demo/demo_data.RData")   # loads object: apgar_temp
data <- apgar_temp

# ---------------------------------------------------------------------------
# 1. Build binary heat-exposure ("heatwave") indicators
#    (same logic as step 3a, restricted to a few thresholds/lags for speed)
# ---------------------------------------------------------------------------
ths  <- c(80, 90)
lags <- 0:2
metric <- "tmean"

for (L in lags) {
  in_col <- paste0(metric, "_pct_lag", L)
  for (th in ths) {
    out_col <- paste0(metric, "_pc", th, "_lag", L)
    data[[out_col]] <- as.integer(data[[in_col]] > th)
  }
}

data <- data %>%
  mutate(HW1 = ifelse(tmean_pc80_lag0 == 1, 1, 0)) %>%
  mutate(HW2 = ifelse(tmean_pc80_lag0 == 1 & tmean_pc80_lag1 == 1, 1, 0)) %>%
  mutate(HW3 = ifelse(tmean_pc90_lag0 == 1, 1, 0)) %>%
  mutate(HW4 = ifelse(tmean_pc90_lag0 == 1 & tmean_pc90_lag1 == 1, 1, 0))

hw_vars <- c("HW1", "HW2", "HW3", "HW4")

# ---------------------------------------------------------------------------
# 2. Fit one conditional logistic regression per heat-exposure definition
#    model: case ~ HWk + strata(strata_id)
# ---------------------------------------------------------------------------
fit_one_local <- function(df, v) {
  dat <- df %>%
    dplyr::select(case, strata_id, dplyr::all_of(v)) %>%
    dplyr::filter(stats::complete.cases(.))

  fml <- stats::as.formula(paste0("case ~ ", v, " + strata(strata_id)"))

  tryCatch({
    m  <- survival::clogit(fml, data = dat)
    sm <- summary(m)
    beta <- sm$coefficients[1, "coef"]
    se   <- sm$coefficients[1, "se(coef)"]
    pval <- sm$coefficients[1, "Pr(>|z|)"]

    tibble(
      heatwave = v, n = stats::nobs(m),
      OR  = exp(beta),
      LCL = exp(beta - 1.96 * se),
      UCL = exp(beta + 1.96 * se),
      p   = pval
    )
  }, error = function(e) {
    tibble(heatwave = v, n = nrow(dat), OR = NA, LCL = NA, UCL = NA, p = NA)
  })
}

results <- purrr::map_dfr(hw_vars, ~ fit_one_local(data, .x))

# ---------------------------------------------------------------------------
# 3. Expected output
# ---------------------------------------------------------------------------
cat("\n=== Demo results (simulated data - illustrative only) ===\n")
print(results)

dir.create("demo/output", showWarnings = FALSE)
write.csv(results, "demo/output/demo_results.csv", row.names = FALSE)
cat("\nSaved to demo/output/demo_results.csv\n")
cat("Demo completed successfully.\n")
