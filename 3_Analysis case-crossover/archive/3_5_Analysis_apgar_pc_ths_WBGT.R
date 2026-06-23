################################################################################
################################################################################
####                            PROJECT:                                    ####
####                       PregnAnZI-2 trial                                ####  
####                       STEP 3. Analysis                                 ####
################################################################################
################################################################################

rm(list =ls())

library(here)
library(dplyr)
library(rgdal)
library(tidyverse)
library(lubridate)
library(ggplot2)
library(dlnm)
library(survival)
library(purrr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)

here::here()    #shows where it thinks the project root is.
here::dr_here() #tells you why 
getwd()         #current working directory when R started.


load("data/PregnAnZI-2 trial/apgar_data_for_analysis_WBGT.RData")

score = 7

data <- apgar_wbgt_temp %>% 
  mutate(country = as.character(country)) %>% 
  #generate low apgar score outcome
  mutate(del_apgar_n = as.numeric(del_apgar_n)) %>% 
  mutate(case = ifelse(del_apgar_n <= score, 1, 0)) %>% 
  #filter to cases with low apgar score
  filter(case == 1) %>% 
  mutate(case = if_else(case_date==date, 1, 0)) %>% 
  #restrict to low-risk births
  subset(del_plural_n == 1) %>%                  #singleton births
  #subset(gesage >= 37) %>%                      #term births (between 37 and 41 weeks) - a lot of missing
  subset(del_wgt_n >= 2.5 & del_wgt_n <= 4) %>%  #optimal birth weight (>=2.5 kg),
  subset(congen == "No") %>%                     #absence of congenital anomalies
  #only keep live birth and fresh stillbirth (removing MSB)
  subset(del_outc_n == "LB" | del_outc_n == "FSB") %>% 
  #generate unique identifier variable for each case and control observation
  mutate(strata_id = as.numeric(factor(paste0(country, randno, idneonate)))) 


# Percentile: exposures parameters
ths  <- c(75, 90, 95, 98)
lags <- 0:2
metric <- c("wbgt_tmax", "wbgt_tmean")

for (v in metric) {
  for (L in lags) {
    in_col <- paste0(v, "_pct_lag", L)
    for (th in ths) {
      out_col <- paste0(v, "_pc", th, "_lag", L)
      data[[out_col]] <- as.integer(data[[in_col]] > th)
    }
  }
}



############ Case-cross over design + conditional logistic regression ##########

# grid of exposure indicators: tmax/tmean x pc x lag
grid <- tidyr::expand_grid(
  pc   = ths,
  lag  = lags,
  temp = metric
) %>%
  mutate(var = sprintf("%s_pc%d_lag%d", temp, pc, lag))


# ---- helpers ----------------------------------------------------------------
# Fit one conditional logit per variable name, using a provided data frame
fit_one_local <- function(df, v) {
  if (!v %in% names(df)) {
    return(tibble::tibble(var = v, coef = NA_real_, se = NA_real_, OR = NA_real_,
                          LCL = NA_real_, UCL = NA_real_, p = NA_real_,
                          n = NA_integer_, converged = FALSE, error = "var not found"))
  }
  
  dat <- df %>%
    dplyr::select(case, strata_id, dplyr::all_of(v)) %>%
    dplyr::filter(stats::complete.cases(.))
  
  if (nrow(dat) == 0L) {
    return(tibble::tibble(var = v, coef = NA_real_, se = NA_real_, OR = NA_real_,
                          LCL = NA_real_, UCL = NA_real_, p = NA_real_,
                          n = 0L, converged = FALSE, error = "no complete cases"))
  }
  
  fml <- stats::as.formula(paste0("case ~ ", v, " + strata(strata_id)"))
  
  tryCatch({
    m  <- survival::clogit(fml, data = dat)
    sm <- summary(m)
    beta <- sm$coefficients[1, "coef"]
    se   <- sm$coefficients[1, "se(coef)"]
    pval <- sm$coefficients[1, "Pr(>|z|)"]
    
    tibble::tibble(
      var = v,
      coef = beta,
      se   = se,
      OR   = exp(beta),
      LCL  = exp(beta - 1.96 * se),
      UCL  = exp(beta + 1.96 * se),
      p    = pval,
      n    = stats::nobs(m),
      converged = isTRUE(m$converged),
      error = NA_character_
    )
  }, error = function(e) {
    tibble::tibble(
      var = v,
      coef = NA_real_, se = NA_real_, OR = NA_real_, LCL = NA_real_, UCL = NA_real_,
      p = NA_real_, n = nrow(dat), converged = FALSE, error = conditionMessage(e)
    )
  })
}

# ---- run_models_for: applies fit_one_local over grid for a given sample ----
run_models_for <- function(df, sample_label) {
  grid %>%
    dplyr::mutate(res = purrr::map(var, ~ fit_one_local(df, .x))) %>%
    dplyr::select(-var) %>%                     # avoid duplicate 'var' on unnest
    tidyr::unnest(cols = c(res)) %>%
    dplyr::rename(exposure = temp, percentile = pc) %>%
    dplyr::mutate(sample = sample_label) %>%
    dplyr::select(sample, exposure, percentile, lag, var, n, coef, se, OR, LCL, UCL, p, converged, error)
}

# ---- define samples ---------------------------------------------------------
# Assumes your main analysis data is in `data` and has a country variable called `country`.
# The filters below are robust to common encodings (name or ISO codes).
data_full <- data

data_bfa <- data %>%
  filter(str_detect(country, regex("^burkina|^bfa$", ignore_case = TRUE)))

data_gmb <- data %>%
  filter(str_detect(country, regex("^gambia|^gmb$", ignore_case = TRUE)))

# ---- run analyses -----------------------------------------------------------
results_full <- run_models_for(data_full, "Full sample")
results_bfa  <- run_models_for(data_bfa,  "Burkina Faso")
results_gmb  <- run_models_for(data_gmb,  "Gambia")

# Combine for convenience
results_all <- bind_rows(results_full, results_bfa, results_gmb)

write.csv(results_all, "results/Results_APGAR_WBGT.csv")

# ---- save results (CSV) -----------------------------------------------------
#write_csv(results_full, "results/results_full_sample_pct.csv")
#write_csv(results_bfa,  "results/results_burkina_faso_pct.csv")
#write_csv(results_gmb,  "results/results_gambia_pct.csv")
#write_csv(results_all,  "results/results_all_samples_combined_pct.csv")

# ---- plotting dataset -------------------------------------------------------
plot_dat <- results_all %>%
  filter(!is.na(OR), !is.na(LCL), !is.na(UCL)) %>%
  # apply requested exclusions
  filter(!(percentile %in% c(80, 85)),
         !(lag %in% 3:5)) %>%
  mutate(
    exposure   = factor(exposure, levels = c("wbgt_tmax", "wbgt_tmean")),
    lag        = factor(lag, levels = sort(unique(lag))),
    percentile = as.numeric(percentile),
    sample     = factor(sample, levels = c("Full sample", "Burkina Faso", "Gambia"))
  )


# ---- plot: separate panels per sample (rows) and exposure (columns) ---------
library(ggplot2)
library(dplyr)

# Set shared limits across panels
y_limits <- c(0.5, 5.0)


# Prepare truncated CI columns but no arrows
plot_dat_limited <- plot_dat %>%
  mutate(
    LCL_plot = pmax(LCL, y_limits[1]),
    UCL_plot = pmin(UCL, y_limits[2])
  )

pd <- position_dodge(width = 0.5)

gg_alt <- ggplot(plot_dat_limited,
                 aes(x = factor(percentile), y = OR, ymin = LCL_plot, ymax = UCL_plot, color = lag)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_errorbar(position = pd, width = 0) +
  geom_point(position = pd, size = 2) +
  facet_grid(sample ~ exposure, switch = "y", scales = "fixed") +
  scale_y_log10(
    position = "right",
    limits = y_limits,
    breaks = c(0.5, 1,2, 5)
  ) +
  labs(
    title = paste("APGAR score \u2264 ", score, sep=""),
    x = "Percentile threshold",
    y = "Odds ratio (log scale)",
    color = "Lag (days)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    axis.title.y.right = element_text(angle = 90)
  )

gg_alt

ggsave(
  filename = paste("figures/Results_apgarU", score, "_pct_WBGT.png", sep=""),        # File name for the image
  plot = gg_alt,                 # The plot object to save
  width = 10,                    # Width of the image in inches
  height = 6,                    # Height of the image in inches
  dpi = 500                      # Resolution in dots per inch
)

