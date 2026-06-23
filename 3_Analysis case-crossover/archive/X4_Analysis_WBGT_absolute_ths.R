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
data_pregnanzi <- apgar_wbgt_temp
rm(list = setdiff(ls(), "data_pregnanzi"))

load("data/PRECISE/apgar_data_for_analysis_WBGT.RData")
data_precise <- apgar_temp
rm(list = setdiff(ls(), c("data_pregnanzi", "data_precise")))

data <- bind_rows(data_pregnanzi, data_precise)
unique(data$sample)

data <- data %>% 
  #remove stillbirths 
  subset(stillbirth == "No") %>% 
  #generate unique identifier variable for each case and control observation
  mutate(strata_id = as.numeric(factor(paste0(sample, strata_id)))) 


################################################################################
### WBGT TMEAN absolute thresholds
################################################################################

# identify the absolute temperatures correspinding to specific percentile thresholds
load("data/PregnAnZI-2 trial/WBGT_data.RData")

ths_tmean_75 <- wbgt_temp %>% filter(wbgt_tmean_pct == 75) %>% 
  group_by(country) %>% summarise(wbgt_tmean = min(wbgt_tmean))

ths_tmean_90 <- wbgt_temp %>% filter(wbgt_tmean_pct == 90) %>% 
  group_by(country) %>% summarise(wbgt_tmean = min(wbgt_tmean))

ths_tmean_95 <- wbgt_temp %>% filter(wbgt_tmean_pct == 95) %>% 
  group_by(country) %>% summarise(wbgt_tmean = min(wbgt_tmean))

# exposures parameters
ths    <- c(25, 26, 27, 28)
lags   <- 0:4
metric <- c("wbgt_tmean")

for (v in metric) {
  for (L in lags) {
    in_col <- paste0(v, "_lag", L)
    for (th in ths) {
      out_col <- paste0(v, "_above", th, "_lag", L)
      data[[out_col]] <- as.integer(data[[in_col]] >= th)
    }
  }
}



### Generate heatwave events
data <- data %>% 
  #above 25 degrees
  mutate(HW1 = ifelse(wbgt_tmean_above25_lag0==1, 1, 0)) %>% 
  mutate(HW2 = ifelse(wbgt_tmean_above25_lag0==1 & wbgt_tmean_above25_lag1==1, 1, 0)) %>% 
  mutate(HW3 = ifelse(wbgt_tmean_above25_lag0==1 & wbgt_tmean_above25_lag1==1 & wbgt_tmean_above25_lag2==1, 1, 0)) %>% 
  mutate(HW4 = ifelse(wbgt_tmean_above25_lag0==1 & wbgt_tmean_above25_lag1==1 & wbgt_tmean_above25_lag2==1 & wbgt_tmean_above25_lag3==1, 1, 0)) %>% 
  mutate(HW5 = ifelse(wbgt_tmean_above25_lag0==1 & wbgt_tmean_above25_lag1==1 & wbgt_tmean_above25_lag2==1 & wbgt_tmean_above25_lag3==1 & wbgt_tmean_above25_lag4==1, 1, 0)) %>%  
  #above 26 degrees
  mutate(HW6 = ifelse(wbgt_tmean_above26_lag0==1, 1, 0)) %>% 
  mutate(HW7 = ifelse(wbgt_tmean_above26_lag0==1 & wbgt_tmean_above26_lag1==1, 1, 0)) %>% 
  mutate(HW8 = ifelse(wbgt_tmean_above26_lag0==1 & wbgt_tmean_above26_lag1==1 & wbgt_tmean_above26_lag2==1, 1, 0)) %>% 
  mutate(HW9 = ifelse(wbgt_tmean_above26_lag0==1 & wbgt_tmean_above26_lag1==1 & wbgt_tmean_above26_lag2==1 & wbgt_tmean_above26_lag3==1, 1, 0)) %>% 
  mutate(HW10 = ifelse(wbgt_tmean_above26_lag0==1 & wbgt_tmean_above26_lag1==1 & wbgt_tmean_above26_lag2==1 & wbgt_tmean_above26_lag3==1 & wbgt_tmean_above26_lag4==1, 1, 0)) %>% 
  #above 27 degrees
  mutate(HW11 = ifelse(wbgt_tmean_above27_lag0==1, 1, 0)) %>% 
  mutate(HW12 = ifelse(wbgt_tmean_above27_lag0==1 & wbgt_tmean_above27_lag1==1, 1, 0)) %>% 
  mutate(HW13 = ifelse(wbgt_tmean_above27_lag0==1 & wbgt_tmean_above27_lag1==1 & wbgt_tmean_above27_lag2==1, 1, 0)) %>% 
  mutate(HW14 = ifelse(wbgt_tmean_above27_lag0==1 & wbgt_tmean_above27_lag1==1 & wbgt_tmean_above27_lag2==1 & wbgt_tmean_above27_lag3==1, 1, 0)) %>% 
  mutate(HW15 = ifelse(wbgt_tmean_above27_lag0==1 & wbgt_tmean_above27_lag1==1 & wbgt_tmean_above27_lag2==1 & wbgt_tmean_above27_lag3==1 & wbgt_tmean_above27_lag4==1, 1, 0)) %>% 
  #above 28 degrees
  mutate(HW16 = ifelse(wbgt_tmean_above28_lag0==1, 1, 0)) %>% 
  mutate(HW17 = ifelse(wbgt_tmean_above28_lag0==1 & wbgt_tmean_above28_lag1==1, 1, 0)) %>% 
  mutate(HW18 = ifelse(wbgt_tmean_above28_lag0==1 & wbgt_tmean_above28_lag1==1 & wbgt_tmean_above28_lag2==1, 1, 0)) %>% 
  mutate(HW19 = ifelse(wbgt_tmean_above28_lag0==1 & wbgt_tmean_above28_lag1==1 & wbgt_tmean_above28_lag2==1 & wbgt_tmean_above28_lag3==1, 1, 0)) %>% 
  mutate(HW20 = ifelse(wbgt_tmean_above28_lag0==1 & wbgt_tmean_above28_lag1==1 & wbgt_tmean_above28_lag2==1 & wbgt_tmean_above28_lag3==1 & wbgt_tmean_above28_lag4==1, 1, 0))

model <- clogit(case ~ HW4 + strata(strata_id), data = data)
summary(model)

######################################################################## 
### Case-crossover design + conditional logistic regression 
### Heatwave indicators: HW1–HW30
### TMEAN
########################################################################

library(dplyr)
library(tidyr)
library(purrr)
library(survival)
library(tibble)

# Vector of heatwave exposure variables
hw_vars <- paste0("HW", 1:20)   # adjust if your names differ

# ---- helpers ----------------------------------------------------------------
# Fit one conditional logit per heatwave variable, using a provided data frame
# Model: case ~ HWk + pre_lag0 + strata(strata_id)
fit_one_local <- function(df, v) {
  
  # make sure both the heatwave variable and pre_lag0 exist
  needed_vars <- c(v, "pre_lag0")
  if (!all(needed_vars %in% names(df))) {
    return(tibble(
      heatwave  = v,
      coef      = NA_real_,
      se        = NA_real_,
      OR        = NA_real_,
      LCL       = NA_real_,
      UCL       = NA_real_,
      p         = NA_real_,
      n         = NA_integer_,
      converged = FALSE,
      error     = paste("var not found:",
                        paste(needed_vars[!needed_vars %in% names(df)], collapse = ", "))
    ))
  }
  
  dat <- df %>%
    dplyr::select(case, strata_id, dplyr::all_of(needed_vars)) %>%
    dplyr::filter(stats::complete.cases(.))
  
  if (nrow(dat) == 0L) {
    return(tibble(
      heatwave  = v,
      coef      = NA_real_,
      se        = NA_real_,
      OR        = NA_real_,
      LCL       = NA_real_,
      UCL       = NA_real_,
      p         = NA_real_,
      n         = 0L,
      converged = FALSE,
      error     = "no complete cases"
    ))
  }
  
  # model: case ~ HWk + pre_lag0 + strata(strata_id)
  fml <- stats::as.formula(paste0("case ~ ", v, " + pre_lag0 + strata(strata_id)"))
  
  tryCatch({
    m  <- survival::clogit(fml, data = dat)
    sm <- summary(m)
    
    # first row is the heatwave variable (since it's first in the formula)
    beta <- sm$coefficients[1, "coef"]
    se   <- sm$coefficients[1, "se(coef)"]
    pval <- sm$coefficients[1, "Pr(>|z|)"]
    
    tibble(
      heatwave  = v,
      coef      = beta,
      se        = se,
      OR        = exp(beta),
      LCL       = exp(beta - 1.96 * se),
      UCL       = exp(beta + 1.96 * se),
      p         = pval,
      n         = stats::nobs(m),
      converged = isTRUE(m$converged),
      error     = NA_character_
    )
  }, error = function(e) {
    tibble(
      heatwave  = v,
      coef      = NA_real_,
      se        = NA_real_,
      OR        = NA_real_,
      LCL       = NA_real_,
      UCL       = NA_real_,
      p         = NA_real_,
      n         = nrow(dat),
      converged = FALSE,
      error     = conditionMessage(e)
    )
  })
}

# ---- run_models_for: applies fit_one_local over all HW variables -----------
run_models_for <- function(df, sample_label) {
  tibble(heatwave = hw_vars) %>%
    dplyr::mutate(res = purrr::map(heatwave, ~ fit_one_local(df, .x))) %>%
    dplyr::select(-heatwave) %>%
    tidyr::unnest(cols = c(res)) %>%
    dplyr::mutate(sample = sample_label) %>%
    dplyr::select(sample, heatwave, n, coef, se, OR, LCL, UCL, p, converged, error)
}

# ---- define samples ---------------------------------------------------------
data_full <- data

data_bfa <- data %>%
  dplyr::filter(sample == "Burkina Faso - PregnAnZI-2 trial")

data_gmb <- data %>%
  dplyr::filter(sample %in% c("The Gambia - PregnAnZI-2 trial",
                              "The Gambia - PRECISE",
                              "The Gambia - PregnAnZI-1"))

# ---- run analyses -----------------------------------------------------------
results_full <- run_models_for(data_full, "Full sample")
results_bfa  <- run_models_for(data_bfa,  "Burkina Faso")
results_gmb  <- run_models_for(data_gmb,  "Gambia")

# Combine for convenience
results_all <- dplyr::bind_rows(results_full, results_bfa, results_gmb)
results_all <- results_all %>%
  mutate(
    duration = case_when(
      heatwave %in% c("HW1", "HW6", "HW11", "HW16") ~ "1 day",
      heatwave %in% c("HW2", "HW7", "HW12", "HW17") ~ "2 days",
      heatwave %in% c("HW3", "HW8", "HW13", "HW18") ~ "3 days",
      heatwave %in% c("HW4", "HW9", "HW14", "HW19") ~ "4 days",
      heatwave %in% c("HW5", "HW10", "HW15", "HW20") ~ "5 days",
      TRUE ~ NA_character_
    )
  ) %>% 
  mutate(
    percentile = case_when(
      heatwave %in% c("HW1", "HW2", "HW3", "HW4", "HW5") ~ "\u226525°C",
      heatwave %in% c("HW6", "HW7", "HW8", "HW9", "HW10") ~ "\u226526°C",
      heatwave %in% c("HW11", "HW12", "HW13", "HW14", "HW15") ~ "\u226527°C",
      heatwave %in% c("HW16", "HW17", "HW18", "HW19", "HW20") ~ "\u226528°C",
      
      TRUE ~ NA_character_
    )
  )

# ---- plotting dataset -------------------------------------------------------
plot_dat <- results_all %>%
  filter(!is.na(OR), !is.na(LCL), !is.na(UCL)) %>%
  mutate(
    heatwave   = factor(heatwave, levels = c("HW1", "HW2", "HW3", "HW4", "HW5", "HW6", "HW7", "HW8", "HW9", "HW10", "HW11", "HW12", "HW13", "HW14", "HW15", "HW16", "HW17", "HW18", "HW19", "HW20")),
    duration   = factor(duration, levels = c("1 day", "2 days", "3 days", "4 days", "5 days")),
    percentile = factor(percentile,levels = c("\u226525°C", "\u226526°C", "\u226527°C", "\u226528°C")),
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
                 aes(x = factor(heatwave), y = OR, ymin = LCL_plot, ymax = UCL_plot, color = duration)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_errorbar(position = pd, width = 0) +
  geom_point(position = pd, size = 2) +
  facet_grid(sample ~ percentile, switch = "y", scales = "free") +  
  scale_y_log10(
    position = "right",
    limits = y_limits,
    breaks = c(0.5, 1, 2, 5)
  ) +
  labs(
    title = paste("APGAR score \u2264 ", 7, sep = ""),
    x = "Heatwave definition",
    y = "Odds ratio (log scale)",
    color = "Heat duration"
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
  filename = paste("figures/Results_apgarU7_abs_WBGTtmean_heatwave.png", sep=""),       
  plot = gg_alt,                 # The plot object to save
  width = 12,                    # Width of the image in inches
  height = 6,                    # Height of the image in inches
  dpi = 500                      # Resolution in dots per inch
)

write.csv(results_all, "results/Results_agbarU7_abs_WBGTtmean_heatwave.csv")

################################################################################
### TMAX absolute thresholds
################################################################################

# identify the absolute temperatures correspinding to specific percentile thresholds
load("data/PregnAnZI-2 trial/WBGT_data.RData")

ths_tmax_75 <- wbgt_temp %>% filter(wbgt_tmax_pct == 75) %>% 
  group_by(country) %>% summarise(wbgt_tmax = min(wbgt_tmax))

ths_tmax_90 <- wbgt_temp %>% filter(wbgt_tmax_pct == 90) %>% 
  group_by(country) %>% summarise(wbgt_tmax = min(wbgt_tmax))

ths_tmax_95 <- wbgt_temp %>% filter(wbgt_tmax_pct == 95) %>% 
  group_by(country) %>% summarise(wbgt_tmax = min(wbgt_tmax))

# exposures parameters
ths    <- c(28, 29, 30, 31)
lags   <- 0:4
metric <- c("wbgt_tmax")

for (v in metric) {
  for (L in lags) {
    in_col <- paste0(v, "_lag", L)
    for (th in ths) {
      out_col <- paste0(v, "_above", th, "_lag", L)
      data[[out_col]] <- as.integer(data[[in_col]] >= th)
    }
  }
}



### Generate heatwave events
data <- data %>% 
  #above 28 degrees
  mutate(HW1 = ifelse(wbgt_tmax_above28_lag0==1, 1, 0)) %>% 
  mutate(HW2 = ifelse(wbgt_tmax_above28_lag0==1 & wbgt_tmax_above28_lag1==1, 1, 0)) %>% 
  mutate(HW3 = ifelse(wbgt_tmax_above28_lag0==1 & wbgt_tmax_above28_lag1==1 & wbgt_tmax_above28_lag2==1, 1, 0)) %>% 
  mutate(HW4 = ifelse(wbgt_tmax_above28_lag0==1 & wbgt_tmax_above28_lag1==1 & wbgt_tmax_above28_lag2==1 & wbgt_tmax_above28_lag3==1, 1, 0)) %>% 
  mutate(HW5 = ifelse(wbgt_tmax_above28_lag0==1 & wbgt_tmax_above28_lag1==1 & wbgt_tmax_above28_lag2==1 & wbgt_tmax_above28_lag3==1 & wbgt_tmax_above28_lag4==1, 1, 0)) %>%  
  #above 29 degrees
  mutate(HW6 = ifelse(wbgt_tmax_above29_lag0==1, 1, 0)) %>% 
  mutate(HW7 = ifelse(wbgt_tmax_above29_lag0==1 & wbgt_tmax_above29_lag1==1, 1, 0)) %>% 
  mutate(HW8 = ifelse(wbgt_tmax_above29_lag0==1 & wbgt_tmax_above29_lag1==1 & wbgt_tmax_above29_lag2==1, 1, 0)) %>% 
  mutate(HW9 = ifelse(wbgt_tmax_above29_lag0==1 & wbgt_tmax_above29_lag1==1 & wbgt_tmax_above29_lag2==1 & wbgt_tmax_above29_lag3==1, 1, 0)) %>% 
  mutate(HW10 = ifelse(wbgt_tmax_above29_lag0==1 & wbgt_tmax_above29_lag1==1 & wbgt_tmax_above29_lag2==1 & wbgt_tmax_above29_lag3==1 & wbgt_tmax_above29_lag4==1, 1, 0)) %>% 
  #above 30 degrees
  mutate(HW11 = ifelse(wbgt_tmax_above30_lag0==1, 1, 0)) %>% 
  mutate(HW12 = ifelse(wbgt_tmax_above30_lag0==1 & wbgt_tmax_above30_lag1==1, 1, 0)) %>% 
  mutate(HW13 = ifelse(wbgt_tmax_above30_lag0==1 & wbgt_tmax_above30_lag1==1 & wbgt_tmax_above30_lag2==1, 1, 0)) %>% 
  mutate(HW14 = ifelse(wbgt_tmax_above30_lag0==1 & wbgt_tmax_above30_lag1==1 & wbgt_tmax_above30_lag2==1 & wbgt_tmax_above30_lag3==1, 1, 0)) %>% 
  mutate(HW15 = ifelse(wbgt_tmax_above30_lag0==1 & wbgt_tmax_above30_lag1==1 & wbgt_tmax_above30_lag2==1 & wbgt_tmax_above30_lag3==1 & wbgt_tmax_above30_lag4==1, 1, 0)) %>% 
  #above 31 degrees
  mutate(HW16 = ifelse(wbgt_tmax_above31_lag0==1, 1, 0)) %>% 
  mutate(HW17 = ifelse(wbgt_tmax_above31_lag0==1 & wbgt_tmax_above31_lag1==1, 1, 0)) %>% 
  mutate(HW18 = ifelse(wbgt_tmax_above31_lag0==1 & wbgt_tmax_above31_lag1==1 & wbgt_tmax_above31_lag2==1, 1, 0)) %>% 
  mutate(HW19 = ifelse(wbgt_tmax_above31_lag0==1 & wbgt_tmax_above31_lag1==1 & wbgt_tmax_above31_lag2==1 & wbgt_tmax_above31_lag3==1, 1, 0)) %>% 
  mutate(HW20 = ifelse(wbgt_tmax_above31_lag0==1 & wbgt_tmax_above31_lag1==1 & wbgt_tmax_above31_lag2==1 & wbgt_tmax_above31_lag3==1 & wbgt_tmax_above31_lag4==1, 1, 0))


model <- clogit(case ~ HW4 + pre_lag0 + strata(strata_id), data = data)
summary(model)


######################################################################## 
### Case-crossover design + conditional logistic regression 
### Heatwave indicators: HW1–HW30
### TMAX
########################################################################

library(dplyr)
library(tidyr)
library(purrr)
library(survival)
library(tibble)

# Vector of heatwave exposure variables
hw_vars <- paste0("HW", 1:20)   # adjust if your names differ

# ---- helpers ----------------------------------------------------------------
# Fit one conditional logit per heatwave variable, using a provided data frame
# Model: case ~ HWk + pre_lag0 + strata(strata_id)
fit_one_local <- function(df, v) {
  
  # make sure both the heatwave variable and pre_lag0 exist
  needed_vars <- c(v, "pre_lag0")
  if (!all(needed_vars %in% names(df))) {
    return(tibble(
      heatwave  = v,
      coef      = NA_real_,
      se        = NA_real_,
      OR        = NA_real_,
      LCL       = NA_real_,
      UCL       = NA_real_,
      p         = NA_real_,
      n         = NA_integer_,
      converged = FALSE,
      error     = paste("var not found:",
                        paste(needed_vars[!needed_vars %in% names(df)], collapse = ", "))
    ))
  }
  
  dat <- df %>%
    dplyr::select(case, strata_id, dplyr::all_of(needed_vars)) %>%
    dplyr::filter(stats::complete.cases(.))
  
  if (nrow(dat) == 0L) {
    return(tibble(
      heatwave  = v,
      coef      = NA_real_,
      se        = NA_real_,
      OR        = NA_real_,
      LCL       = NA_real_,
      UCL       = NA_real_,
      p         = NA_real_,
      n         = 0L,
      converged = FALSE,
      error     = "no complete cases"
    ))
  }
  
  # model: case ~ HWk + pre_lag0 + strata(strata_id)
  fml <- stats::as.formula(paste0("case ~ ", v, " + pre_lag0 + strata(strata_id)"))
  
  tryCatch({
    m  <- survival::clogit(fml, data = dat)
    sm <- summary(m)
    
    # first row is the heatwave variable (since it's first in the formula)
    beta <- sm$coefficients[1, "coef"]
    se   <- sm$coefficients[1, "se(coef)"]
    pval <- sm$coefficients[1, "Pr(>|z|)"]
    
    tibble(
      heatwave  = v,
      coef      = beta,
      se        = se,
      OR        = exp(beta),
      LCL       = exp(beta - 1.96 * se),
      UCL       = exp(beta + 1.96 * se),
      p         = pval,
      n         = stats::nobs(m),
      converged = isTRUE(m$converged),
      error     = NA_character_
    )
  }, error = function(e) {
    tibble(
      heatwave  = v,
      coef      = NA_real_,
      se        = NA_real_,
      OR        = NA_real_,
      LCL       = NA_real_,
      UCL       = NA_real_,
      p         = NA_real_,
      n         = nrow(dat),
      converged = FALSE,
      error     = conditionMessage(e)
    )
  })
}

# ---- run_models_for: applies fit_one_local over all HW variables -----------
run_models_for <- function(df, sample_label) {
  tibble(heatwave = hw_vars) %>%
    dplyr::mutate(res = purrr::map(heatwave, ~ fit_one_local(df, .x))) %>%
    dplyr::select(-heatwave) %>%
    tidyr::unnest(cols = c(res)) %>%
    dplyr::mutate(sample = sample_label) %>%
    dplyr::select(sample, heatwave, n, coef, se, OR, LCL, UCL, p, converged, error)
}

# ---- define samples ---------------------------------------------------------
data_full <- data

data_bfa <- data %>%
  dplyr::filter(sample == "Burkina Faso - PregnAnZI-2 trial")

data_gmb <- data %>%
  dplyr::filter(sample %in% c("The Gambia - PregnAnZI-2 trial",
                              "The Gambia - PRECISE",
                              "The Gambia - PregnAnZI-1"))

# ---- run analyses -----------------------------------------------------------
results_full <- run_models_for(data_full, "Full sample")
results_bfa  <- run_models_for(data_bfa,  "Burkina Faso")
results_gmb  <- run_models_for(data_gmb,  "Gambia")

# Combine for convenience
results_all <- dplyr::bind_rows(results_full, results_bfa, results_gmb)

results_all <- results_all %>%
  mutate(
    duration = case_when(
      heatwave %in% c("HW1", "HW6", "HW11", "HW16") ~ "1 day",
      heatwave %in% c("HW2", "HW7", "HW12", "HW17") ~ "2 days",
      heatwave %in% c("HW3", "HW8", "HW13", "HW18") ~ "3 days",
      heatwave %in% c("HW4", "HW9", "HW14", "HW19") ~ "4 days",
      heatwave %in% c("HW5", "HW10", "HW15", "HW20") ~ "5 days",
      TRUE ~ NA_character_
    )
  ) %>% 
  mutate(
    percentile = case_when(
      heatwave %in% c("HW1", "HW2", "HW3", "HW4", "HW5") ~ "\u226528°C",
      heatwave %in% c("HW6", "HW7", "HW8", "HW9", "HW10") ~ "\u226529°C",
      heatwave %in% c("HW11", "HW12", "HW13", "HW14", "HW15") ~ "\u226530°C",
      heatwave %in% c("HW16", "HW17", "HW18", "HW19", "HW20") ~ "\u226531°C",
      TRUE ~ NA_character_
    )
  )

# ---- plotting dataset -------------------------------------------------------
plot_dat <- results_all %>%
  filter(!is.na(OR), !is.na(LCL), !is.na(UCL)) %>%
  mutate(
    heatwave   = factor(heatwave, levels = c("HW1", "HW2", "HW3", "HW4", "HW5", "HW6", "HW7", "HW8", "HW9", "HW10", "HW11", "HW12", "HW13", "HW14", "HW15", "HW16", "HW17", "HW18", "HW19", "HW20")),
    duration   = factor(duration, levels = c("1 day", "2 days", "3 days", "4 days", "5 days")),
    percentile = factor(percentile,levels = c("\u226528°C", "\u226529°C", "\u226530°C", "\u226531°C")),
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
                 aes(x = factor(heatwave), y = OR, ymin = LCL_plot, ymax = UCL_plot, color = duration)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_errorbar(position = pd, width = 0) +
  geom_point(position = pd, size = 2) +
  facet_grid(sample ~ percentile, switch = "y", scales = "free") +  
  scale_y_log10(
    position = "right",
    limits = y_limits,
    breaks = c(0.5, 1, 2, 5)
  ) +
  labs(
    title = paste("APGAR score \u2264 ", 7, sep = ""),
    x = "Heatwave definition",
    y = "Odds ratio (log scale)",
    color = "Heat duration"
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
  filename = paste("figures/Results_apgarU7_abs_WBGTtmax_heatwave.png", sep=""),       
  plot = gg_alt,                 # The plot object to save
  width = 12,                    # Width of the image in inches
  height = 6,                    # Height of the image in inches
  dpi = 500                      # Resolution in dots per inch
)

write.csv(results_all, "results/Results_agbarU7_abs_WBGTtmax_heatwave.csv")
