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

### Prepare the data for the negative control test

load("data/PRECISE/Tmp_data_ERA5.RData")
load("data/PRECISE/Pre_data.RData")

temp <- temp %>% dplyr::select(-c(year, month, day))
pre <- pre %>% dplyr::select(-c(year, month, day))

load("data/PRECISE/apgar_data_for_analysis.RData")

apgar_full <- apgar_full %>% 
  mutate(date = date + 6)

## Link with the temperature data
apgar_temp <- apgar_full %>% 
  #day of birth
  left_join(temp, by = c("health_facility"="health_facility", "date"="date")) %>% 
  left_join(pre, by = c("health_facility"="health_facility", "date"="date")) %>% 
  rename_with(~ paste0(.x, "_lag0"), .cols = c(tmax:pre)) %>% 
  #lag 1
  mutate(date_lag1 = date-1) %>%
  left_join(temp, by = c("health_facility"="health_facility", "date_lag1"="date")) %>% 
  left_join(pre, by = c("health_facility"="health_facility", "date"="date")) %>% 
  rename_with(~ paste0(.x, "_lag1"), .cols = c(tmax:pre)) %>% 
  #lag 2
  mutate(date_lag2 = date-2) %>%
  left_join(temp, by = c("health_facility"="health_facility", "date_lag2"="date")) %>% 
  left_join(pre, by = c("health_facility"="health_facility", "date"="date")) %>% 
  rename_with(~ paste0(.x, "_lag2"), .cols = c(tmax:pre)) 

save(apgar_temp, file = "data/PRECISE/apgar_data_for_analysis_DBT_ERA5_negative_control.RData")



