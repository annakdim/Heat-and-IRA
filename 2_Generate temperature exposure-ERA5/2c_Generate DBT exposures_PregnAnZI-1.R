################################################################################
################################################################################
####                            PROJECT:                                    ####
####                                                                        ####  
####          STEP 2. Generate Dry-bulb temperature exposures               ####
################################################################################
################################################################################

rm(list =ls())

library(here)
library(haven)   
library(dplyr)
require(sf)
library(sp)
library(raster)
library(viridis)
library(rgdal)
library(tidyverse)
library(lubridate)
library(ggplot2)
#install.packages("rworldmap")
library(rworldmap)
library(zoo)


# Notes:
# Location for Pregnanzi-1 participants:
# Bundung Maternal and Neonatal Hospital, one of the two main Gambian Hospitals in 
# PregnAnZI-2. 

################################################################################
### 2. Generate exposures
### Case-crossover design
################################################################################

rm(list =ls())

# upload the temperature and pre data
load("data/PregnAnZI-2 trial/Tmp_data_ERA5.RData")
load("data/PregnAnZI-2 trial/Pre_data.RData")

temp <- temp %>% dplyr::select(-c(year, month, day))
pre <- pre %>% dplyr::select(-c(year, month, day))

load("data/Pregnanzi-1/apgar_data_for_analysis.RData")

## Link with the temperature data
apgar_temp <- apgar_full %>% 
  #day of birth
  left_join(temp, by = c("country"="country", "date"="date")) %>% 
  left_join(pre, by = c("country"="country", "date"="date")) %>% 
  rename_with(~ paste0(.x, "_lag0"), .cols = c(tmax:pre)) %>% 
  #lag 1
  mutate(date_lag1 = date-1) %>%
  left_join(temp, by = c("country"="country", "date_lag1"="date")) %>% 
  left_join(pre, by = c("country"="country", "date"="date")) %>% 
  rename_with(~ paste0(.x, "_lag1"), .cols = c(tmax:pre)) %>% 
  #lag 2
  mutate(date_lag2 = date-2) %>%
  left_join(temp, by = c("country"="country", "date_lag2"="date")) %>% 
  left_join(pre, by = c("country"="country", "date"="date")) %>% 
  rename_with(~ paste0(.x, "_lag2"), .cols = c(tmax:pre))   

save(apgar_temp, file = "data/Pregnanzi-1/apgar_data_for_analysis_DBT_ERA5.RData")


