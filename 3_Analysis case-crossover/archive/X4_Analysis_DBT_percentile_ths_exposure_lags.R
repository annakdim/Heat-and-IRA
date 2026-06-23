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

load("data/PregnAnZI-2 trial/apgar_data_for_analysis_DBT.RData")
data_pregnanzi2 <- apgar_temp
rm(list = setdiff(ls(), "data_pregnanzi2"))

load("data/PRECISE/apgar_data_for_analysis_DBT.RData")
data_precise <- apgar_temp
rm(list = setdiff(ls(), c("data_pregnanzi2", "data_precise")))

load("data/Pregnanzi-1/apgar_data_for_analysis_DBT.RData")
data_pregnanzi1 <- apgar_temp
rm(list = setdiff(ls(), c("data_pregnanzi2", "data_precise", "data_pregnanzi1")))

data <- bind_rows(data_pregnanzi2, data_precise, data_pregnanzi1)
unique(data$country)

data <- data %>% 
  #remove stillbirths 
  #subset(stillbirth == "No") %>% 
  subset(apgar_score<7) %>% 
  #generate unique identifier variable for each case and control observation
  mutate(strata_id = as.numeric(factor(paste0(sample, strata_id))))

### Generate heatwave events
# exposures parameters
ths  <- c(75, 80, 85, 90, 95, 98)
lags <- 0:4
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


# grid of all combinations of exposure metric, percentile and lag
grid <- expand_grid(
  temp = metric,
  pc   = ths,
  lag  = lags
) %>%
  mutate(
    var = paste0(temp, "_pc", pc, "_lag", lag)  # matches your binary exposure names
  )

# ---- helpers ----------------------------------------------------------------
# Fit one conditional logit per variable name, using a provided data frame
# Now controlling for precipitation at the same lag (e.g. pre_lag0)
fit_one_local <- function(df, v, lag) {
  
  # exposure var (already in grid)
  temp_var <- v
  
  # precipitation control at same lag (adjust name to your data if needed)
  precip_var <- sprintf("pre_lag%d", lag)
  
  # check both variables exist
  needed_vars <- c(temp_var, precip_var)
  if (!all(needed_vars %in% names(df))) {
    return(tibble::tibble(
      var = v, coef = NA_real_, se = NA_real_, OR = NA_real_,
      LCL = NA_real_, UCL = NA_real_, p = NA_real_,
      n = NA_integer_, converged = FALSE,
      error = paste("missing:", paste(needed_vars[!needed_vars %in% names(df)], collapse = ", "))
    ))
  }
  
  dat <- df %>%
    dplyr::select(case, strata_id, dplyr::all_of(needed_vars)) %>%
    dplyr::filter(stats::complete.cases(.))
  
  if (nrow(dat) == 0L) {
    return(tibble::tibble(
      var = v, coef = NA_real_, se = NA_real_, OR = NA_real_,
      LCL = NA_real_, UCL = NA_real_, p = NA_real_,
      n = 0L, converged = FALSE, error = "no complete cases"
    ))
  }
  
  # model: temperature exposure + precipitation covariate
  fml <- stats::as.formula(
    sprintf("case ~ %s + %s + strata(strata_id)", temp_var, precip_var)
  )
  
  tryCatch({
    m  <- survival::clogit(fml, data = dat)
    sm <- summary(m)
    coefs <- sm$coefficients
    
    # explicitly pull the row for the temperature variable
    idx <- which(rownames(coefs) == temp_var)
    if (length(idx) != 1L) {
      stop("temperature coefficient not found or not unique in model")
    }
    
    beta <- coefs[idx, "coef"]
    se   <- coefs[idx, "se(coef)"]
    pval <- coefs[idx, "Pr(>|z|)"]
    
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
    dplyr::mutate(
      # pass both var and lag so we can get pre_lag<lag>
      res = purrr::map2(var, lag, ~ fit_one_local(df, .x, .y))
    ) %>%
    dplyr::select(-var) %>%                     # avoid duplicate 'var' on unnest
    tidyr::unnest(cols = c(res)) %>%
    dplyr::rename(exposure = temp, percentile = pc) %>%
    dplyr::mutate(sample = sample_label) %>%
    dplyr::select(sample, exposure, percentile, lag, var, n,
                  coef, se, OR, LCL, UCL, p, converged, error)
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



# ---- run analyses -----------------------------------------------------------
results_full <- run_models_for(data_full, "Full sample")
results_bfa  <- run_models_for(data_bfa,  "Burkina Faso")
results_gmb  <- run_models_for(data_gmb,  "Gambia")
results_ken  <- run_models_for(data_ken,  "Kenya")
results_urban  <- run_models_for(data_urban,  "Urban")
results_rural  <- run_models_for(data_rural,  "Rural")

# Combine for convenience
results_all <- dplyr::bind_rows(results_full, results_bfa, results_gmb, results_ken, results_urban, results_rural)

# ---- plotting dataset -------------------------------------------------------
plot_dat <- results_all %>%
  filter(!is.na(OR), !is.na(LCL), !is.na(UCL)) %>%
  filter(lag<=2) %>% 
  # apply requested exclusions
  #filter(!(percentile %in% c(80, 85)),
  #       !(lag %in% 3:5)) %>%
  mutate(
    exposure   = factor(exposure, levels = c("tmax", "tmean", "tmin")),
    lag        = factor(lag, levels = sort(unique(lag))),
    percentile = as.numeric(percentile),
    sample     = factor(sample, levels = c("Full sample", "Burkina Faso", "Gambia", "Kenya", "Urban", "Rural"))
  )


# ---- plot: separate panels per sample (rows) and exposure (columns) ---------
library(ggplot2)
library(dplyr)

# Set shared limits across panels
y_limits <- c(0.2, 5.0)


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
    breaks = c(0.2, 0.5, 1,2, 5)
  ) +
  labs(
    title = paste("APGAR score < ", 7, sep=""),
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
  filename = paste("figures/DBT_pct_lags.png", sep=""),       
  plot = gg_alt,                 # The plot object to save
  width = 15,                    # Width of the image in inches
  height = 12,                    # Height of the image in inches
  dpi = 500                      # Resolution in dots per inch
)



