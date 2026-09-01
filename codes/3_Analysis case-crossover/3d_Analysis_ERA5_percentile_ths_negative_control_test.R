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
library(tibble)
library(stringr)
library(readr)
library(tidyr)
library(purrr)
library(survival)


here::here()    #shows where it thinks the project root is.
here::dr_here() #tells you why 
getwd()         #current working directory when R started.

load("data/PregnAnZI-2 trial/apgar_data_for_analysis_DBT_ERA5_negative_control.RData")
data_pregnanzi2 <- apgar_temp
rm(list = setdiff(ls(), "data_pregnanzi2"))

load("data/PRECISE/apgar_data_for_analysis_DBT_ERA5_negative_control.RData")
data_precise <- apgar_temp
rm(list = setdiff(ls(), c("data_pregnanzi2", "data_precise")))

load("data/Pregnanzi-1/apgar_data_for_analysis_DBT_ERA5_negative_control.RData")
data_pregnanzi1 <- apgar_temp
rm(list = setdiff(ls(), c("data_pregnanzi2", "data_precise", "data_pregnanzi1")))

data <- bind_rows(data_pregnanzi2, data_precise, data_pregnanzi1)
unique(data$sample)

data <- data %>% 
  #remove stillbirths 
  #subset(stillbirth == "No") %>% 
  subset(apgar_score<7) %>% 
  subset(country!="Kenya") %>% 
  #subset(birth_weight >= 2.5 & birth_weight <= 4) %>%    #optimal birth weight (>=2.5 kg),
  #generate unique identifier variable for each case and control observation
  mutate(strata_id = as.numeric(factor(paste0(sample, strata_id)))) 


### Generate heatwave events
# exposures parameters
ths  <- c(75, 80, 85, 90, 95, 98)
lags <- 0:2
metric <- c("tmax", "tmean", "tmin")

for (v in metric) {
  for (L in lags) {
    in_col <- paste0(v, "_pct_lag", L)
    for (th in ths) {
      out_col <- paste0(v, "_pc", th, "_lag", L)
      data[[out_col]] <- as.integer(data[[in_col]] > th)
    }
  }
}


# DTR
ths  <- c(2, 5, 10, 15, 20, 25)
lags <- 0:2
metric <- c("dtr")

for (v in metric) {
  for (L in lags) {
    in_col <- paste0(v, "_pct_lag", L)
    for (th in ths) {
      out_col <- paste0(v, "_pc", th, "_lag", L)
      data[[out_col]] <- as.integer(data[[in_col]] < th)
    }
  }
}

########################################################################
### TMEAN percentile thresholds
########################################################################

data <- data %>% 
  #75th percentile
  mutate(HW1 = ifelse(tmean_pc75_lag0==1, 1, 0)) %>% 
  mutate(HW2 = ifelse(tmean_pc75_lag0==1 & tmean_pc75_lag1==1, 1, 0)) %>% 
  mutate(HW3 = ifelse(tmean_pc75_lag0==1 & tmean_pc75_lag1==1 & tmean_pc75_lag2==1, 1, 0)) %>% 
  #mutate(HW4 = ifelse(tmean_pc75_lag0==1 & tmean_pc75_lag1==1 & tmean_pc75_lag2==1 & tmean_pc75_lag3==1, 1, 0)) %>% 
  #mutate(HW5 = ifelse(tmean_pc75_lag0==1 & tmean_pc75_lag1==1 & tmean_pc75_lag2==1 & tmean_pc75_lag3==1 & tmean_pc75_lag4==1, 1, 0)) %>%  
  #80th percentile
  mutate(HW6 = ifelse(tmean_pc80_lag0==1, 1, 0)) %>% 
  mutate(HW7 = ifelse(tmean_pc80_lag0==1 & tmean_pc80_lag1==1, 1, 0)) %>% 
  mutate(HW8 = ifelse(tmean_pc80_lag0==1 & tmean_pc80_lag1==1 & tmean_pc80_lag2==1, 1, 0)) %>% 
  #mutate(HW9 = ifelse(tmean_pc80_lag0==1 & tmean_pc80_lag1==1 & tmean_pc80_lag2==1 & tmean_pc80_lag3==1, 1, 0)) %>% 
  #mutate(HW10 = ifelse(tmean_pc80_lag0==1 & tmean_pc80_lag1==1 & tmean_pc80_lag2==1 & tmean_pc80_lag3==1 & tmean_pc80_lag4==1, 1, 0)) %>%  
  #85th percentile
  mutate(HW11 = ifelse(tmean_pc85_lag0==1, 1, 0)) %>% 
  mutate(HW12 = ifelse(tmean_pc85_lag0==1 & tmean_pc85_lag1==1, 1, 0)) %>% 
  mutate(HW13 = ifelse(tmean_pc85_lag0==1 & tmean_pc85_lag1==1 & tmean_pc85_lag2==1, 1, 0)) %>% 
  #mutate(HW14 = ifelse(tmean_pc85_lag0==1 & tmean_pc85_lag1==1 & tmean_pc85_lag2==1 & tmean_pc85_lag3==1, 1, 0)) %>% 
  #mutate(HW15 = ifelse(tmean_pc85_lag0==1 & tmean_pc85_lag1==1 & tmean_pc85_lag2==1 & tmean_pc85_lag3==1 & tmean_pc85_lag4==1, 1, 0)) %>%  
  #90h percentile
  mutate(HW16 = ifelse(tmean_pc90_lag0==1, 1, 0)) %>% 
  mutate(HW17 = ifelse(tmean_pc90_lag0==1 & tmean_pc90_lag1==1, 1, 0)) %>% 
  mutate(HW18 = ifelse(tmean_pc90_lag0==1 & tmean_pc90_lag1==1 & tmean_pc90_lag2==1, 1, 0)) %>% 
  #mutate(HW19 = ifelse(tmean_pc90_lag0==1 & tmean_pc90_lag1==1 & tmean_pc90_lag2==1 & tmean_pc90_lag3==1, 1, 0)) %>% 
  #mutate(HW20 = ifelse(tmean_pc90_lag0==1 & tmean_pc90_lag1==1 & tmean_pc90_lag2==1 & tmean_pc90_lag3==1 & tmean_pc90_lag4==1, 1, 0)) %>% 
  #95h percentile
  mutate(HW21 = ifelse(tmean_pc95_lag0==1, 1, 0)) %>% 
  mutate(HW22 = ifelse(tmean_pc95_lag0==1 & tmean_pc95_lag1==1, 1, 0)) %>% 
  mutate(HW23 = ifelse(tmean_pc95_lag0==1 & tmean_pc95_lag1==1 & tmean_pc95_lag2==1, 1, 0)) %>% 
  #mutate(HW24 = ifelse(tmean_pc95_lag0==1 & tmean_pc95_lag1==1 & tmean_pc95_lag2==1 & tmean_pc95_lag3==1, 1, 0)) %>% 
  #mutate(HW25 = ifelse(tmean_pc95_lag0==1 & tmean_pc95_lag1==1 & tmean_pc95_lag2==1 & tmean_pc95_lag3==1 & tmean_pc95_lag4==1, 1, 0)) %>% 
  #98h percentile
  mutate(HW26 = ifelse(tmean_pc98_lag0==1, 1, 0)) %>% 
  mutate(HW27 = ifelse(tmean_pc98_lag0==1 & tmean_pc98_lag1==1, 1, 0)) %>% 
  mutate(HW28 = ifelse(tmean_pc98_lag0==1 & tmean_pc98_lag1==1 & tmean_pc98_lag2==1, 1, 0)) 
  #mutate(HW29 = ifelse(tmean_pc98_lag0==1 & tmean_pc98_lag1==1 & tmean_pc98_lag2==1 & tmean_pc98_lag3==1, 1, 0)) %>% 
  #mutate(HW30 = ifelse(tmean_pc98_lag0==1 & tmean_pc98_lag1==1 & tmean_pc98_lag2==1 & tmean_pc98_lag3==1 & tmean_pc98_lag4==1, 1, 0))


# Vector of heatwave exposure variables
hw_vars <- paste0("HW", 1:28)   # adjust if your names differ

# ---- helpers ----------------------------------------------------------------
# Fit one conditional logit per heatwave variable, using a provided data frame
fit_one_local <- function(df, v) {
  if (!v %in% names(df)) {
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
      error     = "var not found"
    ))
  }
  
  dat <- df %>%
    dplyr::select(case, strata_id, dplyr::all_of(v)) %>%
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
  
  # model: case ~ HWk + strata(strata_id)
  fml <- stats::as.formula(paste0("case ~ ", v, " + strata(strata_id)"))
  
  tryCatch({
    m  <- survival::clogit(fml, data = dat)
    sm <- summary(m)
    
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
data_ken <- data %>%
  dplyr::filter(sample == "Kenya - PRECISE")

data_urban <- data %>%
  dplyr::filter(sample %in% c("The Gambia - PregnAnZI-2 trial", "The Gambia - PregnAnZI-1") | (sample=="Kenya - PRECISE" & survey_site=="Kenya  - Mariakani"))

data_rural <- data %>% 
  dplyr::filter(sample %in% c("The Gambia - PRECISE", "Burkina Faso - PregnAnZI-2 trial") | (sample=="Kenya - PRECISE" & survey_site=="Kenya  - Rabai")) 
  
data_season_wet <- data %>% dplyr::filter(season == "Wet")
data_season_dry <- data %>% dplyr::filter(season == "Dry")

data_day <- data %>% dplyr::filter(day_night == "Day")
data_night <- data %>% dplyr::filter(day_night == "Night")

data_twin_no <- data %>% dplyr::filter(twin == "No")
data_twin_yes <- data %>% dplyr::filter(twin == "Yes")

data_cs_no <- data %>% dplyr::filter(cs == "No")
data_cs_yes <- data %>% dplyr::filter(cs == "Yes")

data_lbw_no <- data %>% dplyr::filter(birth_weight >= 2.5)
data_lbw_yes <- data %>% dplyr::filter(birth_weight < 2.5)

# ---- run analyses -----------------------------------------------------------
results_full <- run_models_for(data_full, "Full sample")

results_bfa  <- run_models_for(data_bfa,  "Burkina Faso")
results_gmb  <- run_models_for(data_gmb,  "Gambia")
results_ken  <- run_models_for(data_ken,  "Kenya")

results_urban  <- run_models_for(data_urban,  "Urban")
results_rural  <- run_models_for(data_rural,  "Rural")

results_season_wet  <- run_models_for(data_season_wet,  "Wet season")
results_season_dry  <- run_models_for(data_season_dry,  "Dry season")

results_day  <- run_models_for(data_day,  "Daytime delivery")
results_night  <- run_models_for(data_night,  "Nighttime delivery")

results_twin_no  <- run_models_for(data_twin_no,  "Singleton birth")
results_twin_yes  <- run_models_for(data_twin_yes,  "Multiple birth")

results_cs_no  <- run_models_for(data_cs_no,  "Vaginal birth")
results_cs_yes  <- run_models_for(data_cs_yes,  "Cesarean section")

results_lbw_no  <- run_models_for(data_lbw_no,  "Normal birthweight")
results_lbw_yes  <- run_models_for(data_lbw_yes,  "Low birthweight")

# Combine for convenience
results_all <- dplyr::bind_rows(results_full, results_bfa, results_gmb, results_ken, 
                                results_urban, results_rural, 
                                results_season_wet, results_season_dry, 
                                results_day, results_night,
                                results_twin_no, results_twin_yes,
                                results_cs_no, results_cs_yes,
                                results_lbw_no, results_lbw_yes)

results_all <- results_all %>%
  mutate(
    duration = case_when(
      heatwave %in% c("HW1", "HW6", "HW11", "HW16", "HW21", "HW26") ~ "1 day",
      heatwave %in% c("HW2", "HW7", "HW12", "HW17", "HW22", "HW27") ~ "2 days",
      heatwave %in% c("HW3", "HW8", "HW13", "HW18", "HW23", "HW28") ~ "3 days",
      heatwave %in% c("HW4", "HW9", "HW14", "HW19", "HW24", "HW29") ~ "4 days",
      heatwave %in% c("HW5", "HW10", "HW15", "HW20", "HW25", "HW30") ~ "5 days",
      TRUE ~ NA_character_
    )
  ) %>% 
  mutate(
    percentile = case_when(
      heatwave %in% c("HW1", "HW2", "HW3", "HW4", "HW5") ~ "75th percentile",
      heatwave %in% c("HW6", "HW7", "HW8", "HW9", "HW10") ~ "80th percentile",
      heatwave %in% c("HW11", "HW12", "HW13", "HW14", "HW15") ~ "85th percentile",
      heatwave %in% c("HW16", "HW17", "HW18", "HW19", "HW20") ~ "90th percentile",
      heatwave %in% c("HW21", "HW22", "HW23", "HW24", "HW25") ~ "95th percentile",
      heatwave %in% c("HW26", "HW27", "HW28", "HW29", "HW30") ~ "98th percentile",
      TRUE ~ NA_character_
    )
  )

# ---- plotting dataset -------------------------------------------------------
plot_dat <- results_all %>%
  filter(!is.na(OR), !is.na(LCL), !is.na(UCL)) %>%
  mutate(
    heatwave   = factor(heatwave, levels = c("HW1", "HW2", "HW3", "HW4", "HW5", "HW6", "HW7", "HW8", "HW9", "HW10", "HW11", "HW12", "HW13", "HW14", "HW15", "HW16", "HW17", "HW18", "HW19", "HW20", "HW21", "HW22", "HW23", "HW24", "HW25", "HW26", "HW27", "HW28", "HW29", "HW30")),
    duration   = factor(duration, levels = c("1 day", "2 days", "3 days", "4 days", "5 days")),
    percentile = factor(percentile,levels = c("75th percentile", "80th percentile", "85th percentile","90th percentile", "95th percentile", "98th percentile")),
    sample     = factor(sample, levels = c("Full sample", "Burkina Faso", "Gambia", "Kenya", 
                                           "Urban", "Rural", 
                                           "Wet season", "Dry season", 
                                           "Daytime delivery",  "Nighttime delivery", 
                                           "Singleton birth", "Multiple birth",
                                           "Vaginal birth", "Cesarean section",
                                           "Normal birthweight", "Low birthweight"))
  )


# ---- plot: separate panels per sample (rows) and exposure (columns) ---------

# Set shared limits across panels
y_limits <- c(0.5, 5.0)

# Prepare truncated CI columns but no arrows
plot_dat_limited <- plot_dat %>%
  mutate(
    LCL_plot = pmax(LCL, y_limits[1]),
    UCL_plot = pmin(UCL, y_limits[2])
  ) %>% 
  filter(percentile!="98th percentile")

pd <- position_dodge(width = 0.5)

gg_alt <- ggplot(plot_dat_limited,
                 aes(x = factor(heatwave),
                     y = OR,
                     ymin = LCL_plot,
                     ymax = UCL_plot,
                     color = duration)) +
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
    title = "Daily mean temperature",
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

write.csv(results_all, "results/ApgarU7/ERA5_Tmean_pct_heatwave_negative_control_test.csv")

######################################################################## 
### TMAX percentile thresholds
########################################################################

data <- data %>% 
  #75th percentile
  mutate(HW1 = ifelse(tmax_pc75_lag0==1, 1, 0)) %>% 
  mutate(HW2 = ifelse(tmax_pc75_lag0==1 & tmax_pc75_lag1==1, 1, 0)) %>% 
  mutate(HW3 = ifelse(tmax_pc75_lag0==1 & tmax_pc75_lag1==1 & tmax_pc75_lag2==1, 1, 0)) %>% 
  #mutate(HW4 = ifelse(tmax_pc75_lag0==1 & tmax_pc75_lag1==1 & tmax_pc75_lag2==1 & tmax_pc75_lag3==1, 1, 0)) %>% 
  #mutate(HW5 = ifelse(tmax_pc75_lag0==1 & tmax_pc75_lag1==1 & tmax_pc75_lag2==1 & tmax_pc75_lag3==1 & tmax_pc75_lag4==1, 1, 0)) %>%  
  #80th percentile
  mutate(HW6  = ifelse(tmax_pc80_lag0==1, 1, 0)) %>% 
  mutate(HW7  = ifelse(tmax_pc80_lag0==1 & tmax_pc80_lag1==1, 1, 0)) %>% 
  mutate(HW8  = ifelse(tmax_pc80_lag0==1 & tmax_pc80_lag1==1 & tmax_pc80_lag2==1, 1, 0)) %>% 
  #mutate(HW9  = ifelse(tmax_pc80_lag0==1 & tmax_pc80_lag1==1 & tmax_pc80_lag2==1 & tmax_pc80_lag3==1, 1, 0)) %>% 
  #mutate(HW10 = ifelse(tmax_pc80_lag0==1 & tmax_pc80_lag1==1 & tmax_pc80_lag2==1 & tmax_pc80_lag3==1 & tmax_pc80_lag4==1, 1, 0)) %>%  
  #85th percentile
  mutate(HW11 = ifelse(tmax_pc85_lag0==1, 1, 0)) %>% 
  mutate(HW12 = ifelse(tmax_pc85_lag0==1 & tmax_pc85_lag1==1, 1, 0)) %>% 
  mutate(HW13 = ifelse(tmax_pc85_lag0==1 & tmax_pc85_lag1==1 & tmax_pc85_lag2==1, 1, 0)) %>% 
  #mutate(HW14 = ifelse(tmax_pc85_lag0==1 & tmax_pc85_lag1==1 & tmax_pc85_lag2==1 & tmax_pc85_lag3==1, 1, 0)) %>% 
  #mutate(HW15 = ifelse(tmax_pc85_lag0==1 & tmax_pc85_lag1==1 & tmax_pc85_lag2==1 & tmax_pc85_lag3==1 & tmax_pc85_lag4==1, 1, 0)) %>%  
  #90th percentile
  mutate(HW16 = ifelse(tmax_pc90_lag0==1, 1, 0)) %>% 
  mutate(HW17 = ifelse(tmax_pc90_lag0==1 & tmax_pc90_lag1==1, 1, 0)) %>% 
  mutate(HW18 = ifelse(tmax_pc90_lag0==1 & tmax_pc90_lag1==1 & tmax_pc90_lag2==1, 1, 0)) %>% 
  #mutate(HW19 = ifelse(tmax_pc90_lag0==1 & tmax_pc90_lag1==1 & tmax_pc90_lag2==1 & tmax_pc90_lag3==1, 1, 0)) %>% 
  #mutate(HW20 = ifelse(tmax_pc90_lag0==1 & tmax_pc90_lag1==1 & tmax_pc90_lag2==1 & tmax_pc90_lag3==1 & tmax_pc90_lag4==1, 1, 0)) %>% 
  #95th percentile
  mutate(HW21 = ifelse(tmax_pc95_lag0==1, 1, 0)) %>% 
  mutate(HW22 = ifelse(tmax_pc95_lag0==1 & tmax_pc95_lag1==1, 1, 0)) %>% 
  mutate(HW23 = ifelse(tmax_pc95_lag0==1 & tmax_pc95_lag1==1 & tmax_pc95_lag2==1, 1, 0)) %>% 
  #mutate(HW24 = ifelse(tmax_pc95_lag0==1 & tmax_pc95_lag1==1 & tmax_pc95_lag2==1 & tmax_pc95_lag3==1, 1, 0)) %>% 
  #mutate(HW25 = ifelse(tmax_pc95_lag0==1 & tmax_pc95_lag1==1 & tmax_pc95_lag2==1 & tmax_pc95_lag3==1 & tmax_pc95_lag4==1, 1, 0)) %>% 
  #98th percentile
  mutate(HW26 = ifelse(tmax_pc98_lag0==1, 1, 0)) %>% 
  mutate(HW27 = ifelse(tmax_pc98_lag0==1 & tmax_pc98_lag1==1, 1, 0)) %>% 
  mutate(HW28 = ifelse(tmax_pc98_lag0==1 & tmax_pc98_lag1==1 & tmax_pc98_lag2==1, 1, 0)) 
  #mutate(HW29 = ifelse(tmax_pc98_lag0==1 & tmax_pc98_lag1==1 & tmax_pc98_lag2==1 & tmax_pc98_lag3==1, 1, 0)) %>% 
  #mutate(HW30 = ifelse(tmax_pc98_lag0==1 & tmax_pc98_lag1==1 & tmax_pc98_lag2==1 & tmax_pc98_lag3==1 & tmax_pc98_lag4==1, 1, 0))

# Vector of heatwave exposure variables
hw_vars <- paste0("HW", 1:30)

# ---- helpers ----------------------------------------------------------------
fit_one_local <- function(df, v) {
  if (!v %in% names(df)) {
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
      error     = "var not found"
    ))
  }
  
  dat <- df %>%
    dplyr::select(case, strata_id, dplyr::all_of(v)) %>%
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
  
  fml <- stats::as.formula(paste0("case ~ ", v, " + strata(strata_id)"))
  
  tryCatch({
    m  <- survival::clogit(fml, data = dat)
    sm <- summary(m)
    
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

run_models_for <- function(df, sample_label) {
  tibble(heatwave = hw_vars) %>%
    dplyr::mutate(res = purrr::map(heatwave, ~ fit_one_local(df, .x))) %>%
    dplyr::select(-heatwave) %>%
    tidyr::unnest(cols = c(res)) %>%
    dplyr::mutate(sample = sample_label) %>%
    dplyr::select(sample, heatwave, n, coef, se, OR, LCL, UCL, p, converged, error)
}

# ---- define samples (unchanged) ---------------------------------------------
data_full <- data
data_bfa <- data %>% dplyr::filter(sample == "Burkina Faso - PregnAnZI-2 trial")
data_gmb <- data %>%
  dplyr::filter(sample %in% c("The Gambia - PregnAnZI-2 trial",
                              "The Gambia - PRECISE",
                              "The Gambia - PregnAnZI-1"))
data_ken <- data %>% dplyr::filter(sample == "Kenya - PRECISE")
data_urban <- data %>%
  dplyr::filter(sample %in% c("The Gambia - PregnAnZI-2 trial", "The Gambia - PregnAnZI-1") |
                  (sample=="Kenya - PRECISE" & survey_site=="Kenya  - Mariakani"))
data_rural <- data %>% 
  dplyr::filter(sample %in% c("The Gambia - PRECISE", "Burkina Faso - PregnAnZI-2 trial") |
                  (sample=="Kenya - PRECISE" & survey_site=="Kenya  - Rabai")) 

data_season_wet <- data %>% dplyr::filter(season == "Wet")
data_season_dry <- data %>% dplyr::filter(season == "Dry")

data_day <- data %>% dplyr::filter(day_night == "Day")
data_night <- data %>% dplyr::filter(day_night == "Night")

data_twin_no <- data %>% dplyr::filter(twin == "No")
data_twin_yes <- data %>% dplyr::filter(twin == "Yes")

data_cs_no <- data %>% dplyr::filter(cs == "No")
data_cs_yes <- data %>% dplyr::filter(cs == "Yes")

data_lbw_no <- data %>% dplyr::filter(birth_weight >= 2.5)
data_lbw_yes <- data %>% dplyr::filter(birth_weight < 2.5)

# ---- run analyses -----------------------------------------------------------
results_full   <- run_models_for(data_full,  "Full sample")
results_bfa    <- run_models_for(data_bfa,   "Burkina Faso")
results_gmb    <- run_models_for(data_gmb,   "Gambia")
results_ken    <- run_models_for(data_ken,   "Kenya")

results_urban  <- run_models_for(data_urban, "Urban")
results_rural  <- run_models_for(data_rural, "Rural")

results_season_wet  <- run_models_for(data_season_wet,  "Wet season")
results_season_dry  <- run_models_for(data_season_dry,  "Dry season")

results_day  <- run_models_for(data_day,  "Daytime delivery")
results_night  <- run_models_for(data_night,  "Nighttime delivery")

results_twin_no  <- run_models_for(data_twin_no,  "Singleton birth")
results_twin_yes  <- run_models_for(data_twin_yes,  "Multiple birth")

results_cs_no  <- run_models_for(data_cs_no,  "Vaginal birth")
results_cs_yes  <- run_models_for(data_cs_yes,  "Cesarean section")

results_lbw_no  <- run_models_for(data_lbw_no,  "Normal birthweight")
results_lbw_yes  <- run_models_for(data_lbw_yes,  "Low birthweight")

# Combine for convenience
results_all <- dplyr::bind_rows(results_full, results_bfa, results_gmb, results_ken, 
                                results_urban, results_rural, 
                                results_season_wet, results_season_dry, 
                                results_day, results_night,
                                results_twin_no, results_twin_yes,
                                results_cs_no, results_cs_yes,
                                results_lbw_no, results_lbw_yes)

results_all <- results_all %>%
  mutate(
    duration = case_when(
      heatwave %in% c("HW1", "HW6", "HW11", "HW16", "HW21", "HW26") ~ "1 day",
      heatwave %in% c("HW2", "HW7", "HW12", "HW17", "HW22", "HW27") ~ "2 days",
      heatwave %in% c("HW3", "HW8", "HW13", "HW18", "HW23", "HW28") ~ "3 days",
      heatwave %in% c("HW4", "HW9", "HW14", "HW19", "HW24", "HW29") ~ "4 days",
      heatwave %in% c("HW5", "HW10", "HW15", "HW20", "HW25", "HW30") ~ "5 days",
      TRUE ~ NA_character_
    )
  ) %>% 
  mutate(
    percentile = case_when(
      heatwave %in% c("HW1",  "HW2",  "HW3",  "HW4",  "HW5")  ~ "75th percentile",
      heatwave %in% c("HW6",  "HW7",  "HW8",  "HW9",  "HW10") ~ "80th percentile",
      heatwave %in% c("HW11", "HW12", "HW13", "HW14", "HW15") ~ "85th percentile",
      heatwave %in% c("HW16", "HW17", "HW18", "HW19", "HW20") ~ "90th percentile",
      heatwave %in% c("HW21", "HW22", "HW23", "HW24", "HW25") ~ "95th percentile",
      heatwave %in% c("HW26", "HW27", "HW28", "HW29", "HW30") ~ "98th percentile",
      TRUE ~ NA_character_
    )
  )

# ---- plotting dataset -------------------------------------------------------
plot_dat <- results_all %>%
  filter(!is.na(OR), !is.na(LCL), !is.na(UCL)) %>%
  mutate(
    heatwave   = factor(heatwave, levels = c("HW1", "HW2", "HW3", "HW4", "HW5", "HW6", "HW7", "HW8", "HW9", "HW10", "HW11", "HW12", "HW13", "HW14", "HW15", "HW16", "HW17", "HW18", "HW19", "HW20", "HW21", "HW22", "HW23", "HW24", "HW25", "HW26", "HW27", "HW28", "HW29", "HW30")),
    duration   = factor(duration, levels = c("1 day", "2 days", "3 days", "4 days", "5 days")),
    percentile = factor(percentile,levels = c("75th percentile", "80th percentile", "85th percentile","90th percentile", "95th percentile", "98th percentile")),
    sample     = factor(sample, levels = c("Full sample", "Burkina Faso", "Gambia", "Kenya", 
                                           "Urban", "Rural", 
                                           "Wet season", "Dry season", 
                                           "Daytime delivery",  "Nighttime delivery", 
                                           "Singleton birth", "Multiple birth",
                                           "Vaginal birth", "Cesarean section",
                                           "Normal birthweight", "Low birthweight"))
  )

# ---- plot -------------------------------------------------------------------
y_limits <- c(0.5, 5.0)

plot_dat_limited <- plot_dat %>%
  mutate(
    LCL_plot = pmax(LCL, y_limits[1]),
    UCL_plot = pmin(UCL, y_limits[2])
  ) %>% 
  filter(percentile!="98th percentile")

pd <- position_dodge(width = 0.5)

gg_alt <- ggplot(plot_dat_limited,
                 aes(x = factor(heatwave),
                     y = OR,
                     ymin = LCL_plot,
                     ymax = UCL_plot,
                     color = duration)) +
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
    title = "Daily maximum temperature",
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


write.csv(results_all, "results/ApgarU7/ERA5_Tmax_pct_heatwave_negative_control_test.csv")

######################################################################## 
### TMIN percentile thresholds
########################################################################

data <- data %>% 
  # 75th percentile
  mutate(HW1 = ifelse(tmin_pc75_lag0==1, 1, 0)) %>% 
  mutate(HW2 = ifelse(tmin_pc75_lag0==1 & tmin_pc75_lag1==1, 1, 0)) %>% 
  mutate(HW3 = ifelse(tmin_pc75_lag0==1 & tmin_pc75_lag1==1 & tmin_pc75_lag2==1, 1, 0)) %>% 
  #mutate(HW4 = ifelse(tmin_pc75_lag0==1 & tmin_pc75_lag1==1 & tmin_pc75_lag2==1 & tmin_pc75_lag3==1, 1, 0)) %>% 
  #mutate(HW5 = ifelse(tmin_pc75_lag0==1 & tmin_pc75_lag1==1 & tmin_pc75_lag2==1 & tmin_pc75_lag3==1 & tmin_pc75_lag4==1, 1, 0)) %>%  
  # 80th percentile
  mutate(HW6  = ifelse(tmin_pc80_lag0==1, 1, 0)) %>% 
  mutate(HW7  = ifelse(tmin_pc80_lag0==1 & tmin_pc80_lag1==1, 1, 0)) %>% 
  mutate(HW8  = ifelse(tmin_pc80_lag0==1 & tmin_pc80_lag1==1 & tmin_pc80_lag2==1, 1, 0)) %>% 
  #mutate(HW9  = ifelse(tmin_pc80_lag0==1 & tmin_pc80_lag1==1 & tmin_pc80_lag2==1 & tmin_pc80_lag3==1, 1, 0)) %>% 
  #mutate(HW10 = ifelse(tmin_pc80_lag0==1 & tmin_pc80_lag1==1 & tmin_pc80_lag2==1 & tmin_pc80_lag3==1 & tmin_pc80_lag4==1, 1, 0)) %>%  
  # 85th percentile
  mutate(HW11 = ifelse(tmin_pc85_lag0==1, 1, 0)) %>% 
  mutate(HW12 = ifelse(tmin_pc85_lag0==1 & tmin_pc85_lag1==1, 1, 0)) %>% 
  mutate(HW13 = ifelse(tmin_pc85_lag0==1 & tmin_pc85_lag1==1 & tmin_pc85_lag2==1, 1, 0)) %>% 
  #mutate(HW14 = ifelse(tmin_pc85_lag0==1 & tmin_pc85_lag1==1 & tmin_pc85_lag2==1 & tmin_pc85_lag3==1, 1, 0)) %>% 
  #mutate(HW15 = ifelse(tmin_pc85_lag0==1 & tmin_pc85_lag1==1 & tmin_pc85_lag2==1 & tmin_pc85_lag3==1 & tmin_pc85_lag4==1, 1, 0)) %>%  
  # 90th percentile
  mutate(HW16 = ifelse(tmin_pc90_lag0==1, 1, 0)) %>% 
  mutate(HW17 = ifelse(tmin_pc90_lag0==1 & tmin_pc90_lag1==1, 1, 0)) %>% 
  mutate(HW18 = ifelse(tmin_pc90_lag0==1 & tmin_pc90_lag1==1 & tmin_pc90_lag2==1, 1, 0)) %>% 
  #mutate(HW19 = ifelse(tmin_pc90_lag0==1 & tmin_pc90_lag1==1 & tmin_pc90_lag2==1 & tmin_pc90_lag3==1, 1, 0)) %>% 
  #mutate(HW20 = ifelse(tmin_pc90_lag0==1 & tmin_pc90_lag1==1 & tmin_pc90_lag2==1 & tmin_pc90_lag3==1 & tmin_pc90_lag4==1, 1, 0)) %>% 
  # 95th percentile
  mutate(HW21 = ifelse(tmin_pc95_lag0==1, 1, 0)) %>% 
  mutate(HW22 = ifelse(tmin_pc95_lag0==1 & tmin_pc95_lag1==1, 1, 0)) %>% 
  mutate(HW23 = ifelse(tmin_pc95_lag0==1 & tmin_pc95_lag1==1 & tmin_pc95_lag2==1, 1, 0)) %>% 
  #mutate(HW24 = ifelse(tmin_pc95_lag0==1 & tmin_pc95_lag1==1 & tmin_pc95_lag2==1 & tmin_pc95_lag3==1, 1, 0)) %>% 
  #mutate(HW25 = ifelse(tmin_pc95_lag0==1 & tmin_pc95_lag1==1 & tmin_pc95_lag2==1 & tmin_pc95_lag3==1 & tmin_pc95_lag4==1, 1, 0)) %>% 
  # 98th percentile
  mutate(HW26 = ifelse(tmin_pc98_lag0==1, 1, 0)) %>% 
  mutate(HW27 = ifelse(tmin_pc98_lag0==1 & tmin_pc98_lag1==1, 1, 0)) %>% 
  mutate(HW28 = ifelse(tmin_pc98_lag0==1 & tmin_pc98_lag1==1 & tmin_pc98_lag2==1, 1, 0)) 
  #mutate(HW29 = ifelse(tmin_pc98_lag0==1 & tmin_pc98_lag1==1 & tmin_pc98_lag2==1 & tmin_pc98_lag3==1, 1, 0)) %>% 
  #mutate(HW30 = ifelse(tmin_pc98_lag0==1 & tmin_pc98_lag1==1 & tmin_pc98_lag2==1 & tmin_pc98_lag3==1 & tmin_pc98_lag4==1, 1, 0))

# Vector of heatwave exposure variables
hw_vars <- paste0("HW", 1:30)

# ---- helpers ----------------------------------------------------------------
# Fit one conditional logit per heatwave variable, using a provided data frame
fit_one_local <- function(df, v) {
  if (!v %in% names(df)) {
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
      error     = "var not found"
    ))
  }
  
  dat <- df %>%
    dplyr::select(case, strata_id, dplyr::all_of(v)) %>%
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
  
  # model: case ~ HWk + strata(strata_id)
  fml <- stats::as.formula(paste0("case ~ ", v, " + strata(strata_id)"))
  
  tryCatch({
    m  <- survival::clogit(fml, data = dat)
    sm <- summary(m)
    
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

data_ken <- data %>%
  dplyr::filter(sample == "Kenya - PRECISE")

data_urban <- data %>%
  dplyr::filter(sample %in% c("The Gambia - PregnAnZI-2 trial", "The Gambia - PregnAnZI-1") |
                  (sample=="Kenya - PRECISE" & survey_site=="Kenya  - Mariakani"))

data_rural <- data %>% 
  dplyr::filter(sample %in% c("The Gambia - PRECISE", "Burkina Faso - PregnAnZI-2 trial") |
                  (sample=="Kenya - PRECISE" & survey_site=="Kenya  - Rabai")) 


data_season_wet <- data %>% dplyr::filter(season == "Wet")
data_season_dry <- data %>% dplyr::filter(season == "Dry")

data_day <- data %>% dplyr::filter(day_night == "Day")
data_night <- data %>% dplyr::filter(day_night == "Night")

data_twin_no <- data %>% dplyr::filter(twin == "No")
data_twin_yes <- data %>% dplyr::filter(twin == "Yes")

data_cs_no <- data %>% dplyr::filter(cs == "No")
data_cs_yes <- data %>% dplyr::filter(cs == "Yes")

data_lbw_no <- data %>% dplyr::filter(birth_weight >= 2.5)
data_lbw_yes <- data %>% dplyr::filter(birth_weight < 2.5)

# ---- run analyses -----------------------------------------------------------
results_full   <- run_models_for(data_full,  "Full sample")
results_bfa    <- run_models_for(data_bfa,   "Burkina Faso")
results_gmb    <- run_models_for(data_gmb,   "Gambia")
results_ken    <- run_models_for(data_ken,   "Kenya")
results_urban  <- run_models_for(data_urban, "Urban")
results_rural  <- run_models_for(data_rural, "Rural")

results_season_wet  <- run_models_for(data_season_wet,  "Wet season")
results_season_dry  <- run_models_for(data_season_dry,  "Dry season")

results_day  <- run_models_for(data_day,  "Daytime delivery")
results_night  <- run_models_for(data_night,  "Nighttime delivery")

results_twin_no  <- run_models_for(data_twin_no,  "Singleton birth")
results_twin_yes  <- run_models_for(data_twin_yes,  "Multiple birth")

results_cs_no  <- run_models_for(data_cs_no,  "Vaginal birth")
results_cs_yes  <- run_models_for(data_cs_yes,  "Cesarean section")

results_lbw_no  <- run_models_for(data_lbw_no,  "Normal birthweight")
results_lbw_yes  <- run_models_for(data_lbw_yes,  "Low birthweight")

# Combine for convenience
results_all <- dplyr::bind_rows(results_full, results_bfa, results_gmb, results_ken, 
                                results_urban, results_rural, 
                                results_season_wet, results_season_dry, 
                                results_day, results_night,
                                results_twin_no, results_twin_yes,
                                results_cs_no, results_cs_yes,
                                results_lbw_no, results_lbw_yes)

results_all <- results_all %>%
  mutate(
    duration = case_when(
      heatwave %in% c("HW1", "HW6", "HW11", "HW16", "HW21", "HW26") ~ "1 day",
      heatwave %in% c("HW2", "HW7", "HW12", "HW17", "HW22", "HW27") ~ "2 days",
      heatwave %in% c("HW3", "HW8", "HW13", "HW18", "HW23", "HW28") ~ "3 days",
      heatwave %in% c("HW4", "HW9", "HW14", "HW19", "HW24", "HW29") ~ "4 days",
      heatwave %in% c("HW5", "HW10", "HW15", "HW20", "HW25", "HW30") ~ "5 days",
      TRUE ~ NA_character_
    )
  ) %>% 
  mutate(
    percentile = case_when(
      heatwave %in% c("HW1",  "HW2",  "HW3",  "HW4",  "HW5")  ~ "75th percentile",
      heatwave %in% c("HW6",  "HW7",  "HW8",  "HW9",  "HW10") ~ "80th percentile",
      heatwave %in% c("HW11", "HW12", "HW13", "HW14", "HW15") ~ "85th percentile",
      heatwave %in% c("HW16", "HW17", "HW18", "HW19", "HW20") ~ "90th percentile",
      heatwave %in% c("HW21", "HW22", "HW23", "HW24", "HW25") ~ "95th percentile",
      heatwave %in% c("HW26", "HW27", "HW28", "HW29", "HW30") ~ "98th percentile",
      TRUE ~ NA_character_
    )
  )

# ---- plotting dataset -------------------------------------------------------
plot_dat <- results_all %>%
  filter(!is.na(OR), !is.na(LCL), !is.na(UCL)) %>%
  mutate(
    heatwave   = factor(heatwave, levels = c("HW1", "HW2", "HW3", "HW4", "HW5", "HW6", "HW7", "HW8", "HW9", "HW10", "HW11", "HW12", "HW13", "HW14", "HW15", "HW16", "HW17", "HW18", "HW19", "HW20", "HW21", "HW22", "HW23", "HW24", "HW25", "HW26", "HW27", "HW28", "HW29", "HW30")),
    duration   = factor(duration, levels = c("1 day", "2 days", "3 days", "4 days", "5 days")),
    percentile = factor(percentile,levels = c("75th percentile", "80th percentile", "85th percentile","90th percentile", "95th percentile", "98th percentile")),
    sample     = factor(sample, levels = c("Full sample", "Burkina Faso", "Gambia", "Kenya", 
                                           "Urban", "Rural", 
                                           "Wet season", "Dry season", 
                                           "Daytime delivery",  "Nighttime delivery", 
                                           "Singleton birth", "Multiple birth",
                                           "Vaginal birth", "Cesarean section",
                                           "Normal birthweight", "Low birthweight"))
  )

# ---- plot: separate panels per sample (rows) and percentile (columns) ----

y_limits <- c(0.5, 5.0)

plot_dat_limited <- plot_dat %>%
  mutate(
    LCL_plot = pmax(LCL, y_limits[1]),
    UCL_plot = pmin(UCL, y_limits[2])
  ) %>% 
  filter(percentile!="98th percentile")

pd <- position_dodge(width = 0.5)

gg_alt <- ggplot(plot_dat_limited,
                 aes(x = factor(heatwave),
                     y = OR,
                     ymin = LCL_plot,
                     ymax = UCL_plot,
                     color = duration)) +
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
    title = "Daily minimum temperature",
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

write.csv(results_all, "results/ApgarU7/ERA5_Tmin_pct_heatwave_negative_control_test.csv")

