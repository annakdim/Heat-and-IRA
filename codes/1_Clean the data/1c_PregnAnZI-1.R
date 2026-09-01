################################################################################
################################################################################
####                            PROJECT:                                    ####
####                                                                        ####  
####                     STEP 1. Clean the data                             ####
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

# Available variables
# Date of delivery: birthdate
# Apgar score: apgar
# Delivery outcome (alive or stillbirth): birthalv
# Congenital Malformation (we may not have this information though): SevCong (excel file)
# Birth weight: bwgt
# Singleton/twin: multpreg

# upload the survey data
data <- read_dta("data/Pregnanzi-1/baseline.dta")

# all labelled variables - factors with labels
data <- data %>% mutate(across(where(haven::is.labelled),
                               ~ haven::as_factor(., levels = "labels")))

# birth outcome - few cases of stillbirth
data %>% group_by(birthalv) %>% summarise(count = n()) 

# no information on severe congenital malformation
#data_visit0 <- read.csv("data/Pregnanzi-1/01_visit0.csv")
#unique(data_visit0$sevcong_e1_c1)

data  %>% 
  group_by(studyid, birthno) %>% 
  summarise(count = n()) %>% 
  filter(count>1)

unique(data$birthalv)
hist(data$bwgt)

var_list <- c("studyid","birthno", "birthdate", "season",
               'gesage', "bwgt", "multpreg", "birthalv",
               "delvmode", "datetimebirth",
               #vars for apgar score
               "apgar")

apgar <- data %>% 
  dplyr::select(all_of(var_list)) %>% 
  mutate(country = "Gambia") %>% 
  filter(!is.na(birthdate)) %>%                #date of birth is not missing
  #generate case date (delivery date): year, month and day of the week
  mutate(case_date = birthdate) %>%              
  mutate(year = as.numeric(format(case_date,'%Y'))) %>% 
  mutate(month = as.numeric(format(case_date,'%m'))) %>% 
  dplyr::mutate(week_day = lubridate::wday(case_date, label = TRUE)) %>% 
  # add time of delivery
  mutate(
    datetime_clean = trimws(datetimebirth),
    datetime_parsed = as.POSIXct(datetime_clean, format = "%Y-%m-%d %H:%M")
  ) %>% 
  mutate(delivery_hour = format(datetime_parsed, "%H:%M")) %>% 
  mutate(
    hour_num = as.numeric(format(datetime_parsed, "%H")),
    day_night = if_else(hour_num >= 7 & hour_num <= 18, "Day", "Night")
  ) %>% 
  #generate low apgar score outcome
  mutate(apgar_score = as.numeric(apgar)) %>% 
  mutate(twin = if_else(multpreg=="No", "No", "Yes")) %>%  #singleton birth (0) or twin (1)
  mutate(cs = if_else(delvmode=="Caesarean", "Yes", "No")) %>%    #mode of delivery: vaginal (0) or Caesarean (1)
  #subset(gesage >= 37) %>%                      #term births (between 37 and 41 weeks) - a lot of missing
  #subset(bwgt >= 2.5 & bwgt <= 4) %>%  #optimal birth weight (>=2.5 kg),
  mutate(birth_weight = bwgt) %>% 
  #subset(congen == "No") %>%           #absence of congenital anomalies
  #only keep live birth and fresh stillbirth (removing MSB)
  #subset(del_outc_n == "LB" | del_outc_n == "FSB") %>% 
  mutate(sample = "The Gambia - PregnAnZI-1") %>% 
  mutate(survey_site = sample) %>%
  mutate(stillbirth = ifelse(birthalv %in% c("No"), "Yes", "No")) %>% 
  mutate(season = ifelse(season=="Rainy", "Wet", "Dry")) %>% 
  #generate unique identifier variable for each case and control observation
  mutate(strata_id = as.numeric(factor(paste0(studyid, birthno)))) %>% 
  dplyr::select(country, sample, survey_site, strata_id, apgar_score, case_date, year, month, week_day, delivery_hour, day_night, stillbirth, birth_weight, twin, cs, season)

## Create a calendar (spanning the delivery period in the data)
min(unique(data$birthdate))
max(unique(data$birthdate))

cal <- as.data.frame(seq(as.Date("2013/04/01"), as.Date("2014/04/30"), by="day"))
colnames(cal) <- c("date")
cal <- cal %>% 
  dplyr::mutate(week_day = lubridate::wday(date, label = TRUE)) %>%  #add day of the week
  mutate(year = as.numeric(format(date,'%Y'))) %>% 
  mutate(month = as.numeric(format(date,'%m'))) 

apgar_full <- apgar %>% 
  #add control days matched by the year, month and week day  
  left_join(cal, by = c("year", "month", "week_day")) %>% 
  mutate(case = if_else(case_date==date, 1, 0)) %>% 
  dplyr::select(-c(year, month, week_day))

save(apgar_full, file = "data/Pregnanzi-1/apgar_data_for_analysis.RData")

