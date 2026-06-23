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

# upload the survey data

data <- read.csv("data/PRECISE/Gambia_updated_Dec_2025.csv")
names(data)

#data_kenya <- read.csv("data/PRECISE/KenyaPRECISE.csv")
#names(data_kenya)
#data_kenya <- data_kenya %>% rename(f2a_participant_id = ï..f2a_participant_id)
#data <- bind_rows(data, data_kenya) 

unique(data$country)
unique(data$health_facility)

data %>% group_by(country, health_facility) %>% summarise(count = n()) # 2 NA values for Gambia; 58 NA values for Kenya

data <- data %>% rename(participant_id = f2a_participant_id)

data$delivery_date <- as.Date(data$delivery_date, format = "%d-%b-%y")

var_list <- c("country", "health_facility","childID", 
               "delivery_date", "datetime_of_Delivery",
               "stillbirth", "Birthweight", "delivery_num_of_babies", 
               "f9_baby_malformation", "f9_baby_breathe_at_birth",
               "f9_difficulty_breathing_at_birth_resuscitation",
               "delivery_mode", 
               #vars for apgar score
               "apgarscore_5min", "apgarscore_1min")


unique(data$delivery_mode)

apgar <- data %>% 
  dplyr::select(all_of(var_list)) %>% 
  filter(!is.na(health_facility)) %>%      #location is not missing
  filter(!is.na(delivery_date)) %>%        #date of birth is not missing
  #generate case date (delivery date): year, month and day of the week
  mutate(case_date = delivery_date) %>%              
  mutate(year = as.numeric(format(case_date,'%Y'))) %>% 
  mutate(month = as.numeric(format(case_date,'%m'))) %>% 
  dplyr::mutate(week_day = lubridate::wday(case_date, label = TRUE)) %>% 
  mutate(season = if_else(month %in% c(6, 7, 8, 9, 10), "Wet", "Dry")) %>% 
  # add time of delivery  
  mutate(
    datetime_clean = trimws(datetime_of_Delivery),
    datetime_parsed = as.POSIXct(datetime_clean, format = "%d/%m/%Y %H:%M")
  ) %>% 
  mutate(delivery_hour = format(datetime_parsed, "%H:%M")) %>% 
  mutate(
    hour_num = as.numeric(format(datetime_parsed, "%H")),
    day_night = if_else(hour_num >= 7 & hour_num <= 18, "Day", "Night")
  ) %>% 
  #generate low apgar score outcome
  mutate(apgar_score = as.numeric(apgarscore_1min)) %>%
  mutate(apgar_score = ifelse(stillbirth == "Yes", 0 , apgar_score)) %>%
  filter(!is.na(apgar_score)) %>% 
  mutate(twin = if_else(delivery_num_of_babies=="Singleton", "No", "Yes")) %>%  #singleton birth (0) or twin (1)
  mutate(cs = if_else(delivery_mode=="Caesarean section", "Yes", "No")) %>%    #mode of delivery: vaginal or cesarean section
  mutate(birth_weight = Birthweight/1000) %>% 
  #restrict to low-risk births
  #subset(gesage >= 37) %>%                    #term births (between 37 and 41 weeks) - a lot of missing
  subset(f9_baby_malformation != "Yes") %>%    #absence of congenital anomalies
  mutate(sample = paste0(country, " - PRECISE")) %>% 
  mutate(survey_site = paste0(country, "  - ", health_facility)) %>%
  #generate unique identifier variable for each case and control observation
  mutate(strata_id = as.numeric(factor(paste0(health_facility, childID)))) %>% 
  dplyr::select(country, sample, survey_site, health_facility, strata_id, apgar_score, case_date, year, month, week_day, delivery_hour, day_night, stillbirth, birth_weight, twin, cs, season)


apgar <- apgar %>% 
  mutate(season = ifelse(country=="Kenya", NA, season))


## Create a calendar (spanning the delivery period in the data)
min(unique(data$delivery_date))
max(unique(data$delivery_date))

cal <- as.data.frame(seq(as.Date("2019-07-01"), as.Date("2023-04-30"), by="day"))
colnames(cal) <- c("date")
cal <- cal %>% 
  dplyr::mutate(week_day = lubridate::wday(date, label = TRUE)) %>%  #add day of the week
  mutate(year = as.numeric(format(date,'%Y'))) %>% 
  mutate(month = as.numeric(format(date,'%m'))) 

apgar_full <- apgar %>% 
  #add control days matched by the year, month and week day  
  left_join(cal, by = c("year", "month", "week_day"))   %>% 
  mutate(case = if_else(case_date==date, 1, 0)) %>% 
  dplyr::select(-c(year, month, week_day))

save(apgar_full, file = "data/PRECISE/apgar_data_for_analysis.RData")

