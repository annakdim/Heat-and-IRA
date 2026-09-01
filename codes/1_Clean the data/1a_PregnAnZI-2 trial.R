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
data <- read_dta("data/PregnAnZI-2 trial/PregnAnZI_2 Subgroups.dta")
var_info <- data.frame(
  variable = names(data),
  label = sapply(data, function(x) {
    lbl <- attr(x, "label")
    if (is.null(lbl)) "" else lbl   # replace NULL with empty string
  }),
  stringsAsFactors = FALSE
)

rownames(var_info) <- NULL

# all labelled variables - factors with labels
data <- data %>% mutate(across(where(haven::is.labelled),
                               ~ haven::as_factor(., levels = "labels")))

unique(data$season)

# add information on the time of delivery
df1 <- read.csv("data/PregnAnZI-2 trial/visit0_delivery_1.csv")
df1 <- df1 %>% dplyr::select("X1..Numéro.de.Randomisation", "X34..Heure.d.accouchement.Singleton.1er.Bébé")
df1 <- df1 %>% 
  rename("randno" = "X1..Numéro.de.Randomisation",
         "delivery_hour" = "X34..Heure.d.accouchement.Singleton.1er.Bébé")
  
df2 <- read.csv("data/PregnAnZI-2 trial/visit0_delivery_2.csv")
df2 <- df2 %>% dplyr::select("X1..Rand.No", "X34..Singleton.1St.Baby.Delivery.Time")
df2 <- df2 %>% 
  rename("randno" = "X1..Rand.No",
         "delivery_hour" = "X34..Singleton.1St.Baby.Delivery.Time")

df3 <- rbind(df1, df2) %>% 
  mutate(
    time_clean = trimws(delivery_hour),
    hour_num = as.numeric(substr(time_clean, 1, 2))   # extract HH as numeric
  ) %>% 
  mutate(
    day_night = if_else(hour_num >= 7 & hour_num <= 18, "Day", "Night")
  )

data <- data %>% 
  left_join(df3) 

# create a list of variable that will be used for the analysis
var_list <- c("country","randno", "idneonate", "del_dat_n", "season",
               "gesage", "del_wgt_n", "congen", "del_plural_n", "del_outc_n",
               "del_mod_n", "delivery_hour", "day_night",
               #vars for apgar score
               "del_apgar_n")


apgar <- data %>% 
  dplyr::select(all_of(var_list)) %>% 
  mutate(country = as.character(country)) %>% 
  # no missing
  filter(!is.na(del_dat_n)) %>% 
  filter(!is.na(del_apgar_n)) %>%     #remove if date of birth is missing
  #generate case date (delivery date): year, month and day of the week
  mutate(case_date = del_dat_n) %>%              
  mutate(year = as.numeric(format(case_date,'%Y'))) %>% 
  mutate(month = as.numeric(format(case_date,'%m'))) %>% 
  dplyr::mutate(week_day = lubridate::wday(case_date, label = TRUE)) %>% 
  #generate apgar score 
  mutate(apgar_score = as.numeric(del_apgar_n)) %>% 
  mutate(twin = if_else(del_plural_n==1, "No", "Yes")) %>%  #twin or singleton births
  mutate(cs = if_else(del_mod_n=="CS", "Yes", "No")) %>%    #mode of delivery: vaginal or cesarean section
  mutate(birth_weight = del_wgt_n) %>% 
  #restrict to low-risk births
  #subset(gesage >= 37) %>%                      #term births (between 37 and 41 weeks) - a lot of missing
  subset(congen == "No") %>%                     #absence of congenital anomalies
  #only keep live birth and fresh stillbirth (removing MSB)
  subset(del_outc_n == "LB" | del_outc_n == "FSB") %>% 
  mutate(sample = ifelse(country == "Gambia", "The Gambia - PregnAnZI-2 trial", "Burkina Faso - PregnAnZI-2 trial")) %>% 
  mutate(survey_site = sample) %>%
  mutate(stillbirth = ifelse(del_outc_n %in% c("FSB", "MSB"), "Yes", "No")) %>% 
  mutate(season = ifelse(season=="Wet (Jun-Oct)", "Wet", "Dry")) %>% 
  #generate unique identifier variable for each case and control observation
  mutate(strata_id = as.numeric(factor(paste0(country, randno, idneonate)))) %>% 
  dplyr::select(country, sample, survey_site, strata_id, apgar_score, case_date, year, month, week_day, delivery_hour, day_night, stillbirth, birth_weight, twin, cs, season)
 

## Create a calendar (spanning the delivery period in the data)
min(unique(data$del_dat_n))
max(unique(data$del_dat_n))

cal <- as.data.frame(seq(as.Date("2017/10/01"), as.Date("2021/04/30"), by="day"))
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

save(apgar_full, file = "data/PregnAnZI-2 trial/apgar_data_for_analysis.RData")

