################################################################################
################################################################################
####                            PROJECT: Heat-and-IRA                       ####
####            DEMO: Simulate a small synthetic dataset                    ####
####                                                                        ####
####  This script does NOT use any real participant data. It generates a   ####
####  small, fully synthetic dataset with the same variable names and      ####
####  structure expected by the analysis scripts in                        ####
####  "3_Analysis case-crossover/", so that the code can be installed,     ####
####  run, and checked end-to-end without access to the restricted         ####
####  clinical data described in the manuscript's Data Availability        ####
####  statement.                                                           ####
####                                                                        ####
####  Output: demo/demo_data.RData (object: apgar_temp)                    ####
################################################################################
################################################################################

rm(list = ls())
set.seed(123)

library(dplyr)

# ---------------------------------------------------------------------------
# 1. Simulation parameters
# ---------------------------------------------------------------------------
n_cases          <- 40     # number of simulated IRA "case" days
referents_per_case <- 3    # referent (control) days matched to each case
samples <- c("Burkina Faso - PregnAnZI-2 trial",
             "The Gambia - PregnAnZI-2 trial",
             "The Gambia - PRECISE",
             "The Gambia - PregnAnZI-1")

metrics <- c("tmax", "tmean", "tmin", "dtr")
lags    <- 0:4

# ---------------------------------------------------------------------------
# 2. Build one row per case/referent observation
# ---------------------------------------------------------------------------
rows <- list()
row_i <- 1

for (i in seq_len(n_cases)) {

  strata_id <- i
  sample_i  <- sample(samples, 1)
  country_i <- ifelse(grepl("Burkina Faso", sample_i), "Burkina Faso", "Gambia")

  # simulate one case day + N referent days sharing the same strata_id,
  # mirroring the time-stratified case-crossover design used in the study
  n_obs <- 1 + referents_per_case
  case_flag <- c(1, rep(0, referents_per_case))

  for (j in seq_len(n_obs)) {

    row <- list(
      sample      = sample_i,
      country     = country_i,
      strata_id   = strata_id,
      case        = case_flag[j],
      stillbirth  = sample(c("No", "Yes"), 1, prob = c(0.85, 0.15)),
      apgar_score = ifelse(case_flag[j] == 1,
                            sample(0:6, 1),      # cases: Apgar < 7 by definition
                            sample(7:10, 1)),     # referents: healthy Apgar
      birth_weight = round(rnorm(1, mean = 3.0, sd = 0.5), 2)
    )

    # simulated temperature percentile exposures (0-100 scale) for each
    # metric and lag, analogous to the *_pct_lag* columns produced by
    # "2_Generate temperature exposure-ERA5/"
    for (m in metrics) {
      for (L in lags) {
        col <- paste0(m, "_pct_lag", L)
        # cases get a modest upward shift to create a plausible demo signal
        shift <- ifelse(case_flag[j] == 1 & m %in% c("tmax", "tmean"), 12, 0)
        val <- min(100, max(0, rnorm(1, mean = 55 + shift, sd = 20)))
        row[[col]] <- val
      }
    }

    rows[[row_i]] <- row
    row_i <- row_i + 1
  }
}

apgar_temp <- dplyr::bind_rows(rows)

# ---------------------------------------------------------------------------
# 3. Save in the same format expected by the analysis scripts
# ---------------------------------------------------------------------------
dir.create("demo", showWarnings = FALSE)
save(apgar_temp, file = "demo/demo_data.RData")

cat("Simulated demo dataset created: demo/demo_data.RData\n")
cat("Rows:", nrow(apgar_temp), " | Cases:", sum(apgar_temp$case == 1),
    " | Referents:", sum(apgar_temp$case == 0), "\n")
