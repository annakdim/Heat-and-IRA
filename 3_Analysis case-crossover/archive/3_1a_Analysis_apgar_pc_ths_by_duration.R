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


load("data/apgar_data_for_analysis.RData")

data <- apgar_temp %>% 
  mutate(country = as.character(country)) %>% 
  #generate low apgar score outcome
  mutate(del_apgar_n = as.numeric(del_apgar_n)) %>% 
  mutate(case = ifelse(del_apgar_n <= 7, 1, 0)) %>% 
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


# exposures parameters
ths  <- c(75, 90, 95, 98)
lags <- 0:7
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


# exposures parameters
ths    <- c(28, 29, 30, 31, 32, 33, 34, 35)
lags   <- 0:5
metric <- c("tmean")

for (v in metric) {
  for (L in lags) {
    in_col <- paste0(v, "_lag", L)
    for (th in ths) {
      out_col <- paste0(v, "_above", th, "_lag", L)
      data[[out_col]] <- as.integer(data[[in_col]] > th)
    }
  }
}


# exposures parameters
ths  <- c(36, 37, 38, 39, 40, 41, 42)
lags <- 0:5
metric <- c("tmax")

for (v in metric) {
  for (L in lags) {
    in_col <- paste0(v, "_lag", L)
    for (th in ths) {
      out_col <- paste0(v, "_above", th, "_lag", L)
      data[[out_col]] <- as.integer(data[[in_col]] > th)
    }
  }
}

# Generate extreme heat event by duration

data <- data %>% 
  filter(country=="Burkina Faso") %>% 
  #75th percentile
  mutate(HW1 = ifelse(tmean_above33_lag0==1, 1, 0)) %>% 
  mutate(HW2 = ifelse(tmean_above33_lag0==1 & tmean_above33_lag1==1, 1, 0)) %>% 
  mutate(HW3 = ifelse(tmean_above33_lag0==1 & tmean_above33_lag1==1 & tmean_above33_lag2==1, 1, 0)) %>% 
  mutate(HW4 = ifelse(tmean_above33_lag0==1 & tmean_above33_lag1==1 & tmean_above33_lag2==1 & tmean_above33_lag3==1, 1, 0)) %>% 
  mutate(HW5 = ifelse(tmean_above33_lag0==1 & tmean_above33_lag1==1 & tmean_above33_lag2==1 & tmean_above33_lag3==1 & tmean_above33_lag4==1, 1, 0)) %>% 
  #combine into single variable
  mutate(HW_duration = HW1 + HW2 + HW3 + HW4 + HW5) %>% 
  mutate(HW_duration = ifelse(HW_duration>2, 2, HW_duration))

unique(data$HW_duration)

data %>% group_by(HW_duration) %>% summarise(count = n())

model <- clogit(case ~ as.factor(HW_duration) + strata(strata_id), data = data)
summary(model) 

