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
#install.packages("rgdal")
library(rgdal)
library(tidyverse)
library(lubridate)
library(ggplot2)
#install.packages("rworldmap")
library(rworldmap)
library(zoo)


pts <- data.frame(country = c("Gambia", "Burkina Faso"), 
                   long = c(-16.67839, -2.191661), 
                   lat = c(13.42459, 12.67533))

# create a spatial points dataframe
sp <- pts 
coordinates(sp) <- c("long", "lat")
proj4string(sp) <- CRS("+proj=longlat +datum=WGS84")
plot(sp)

# Points in lon/lat (EPSG:4326)
pts_sf <- st_as_sf(sp) %>% st_set_crs(4326)
sf::sf_use_s2(TRUE)
buf_15km <- st_buffer(pts_sf, dist = 15000)     # geodesic, meters

# get the country boundaries
bound0 <- getMap(resolution = "low")
bound1 <- bound0[bound0@data$NAME %in% c("Burkina Faso", "Gambia"), ]
bound2 <- buffer(bound1, width=0.5) 

plot(bound1)
plot(bound2, add=T)
plot(buf_15km, add=T)
plot(sp, add=T)


################################################################################
### 1a. Link the locations with daily temperature data
### Temperature data source: NOAA
################################################################################

year = c(2000:2021) # select a longer time period to identify percentile distribution
y=2020

tmax = NULL

for (y in year){
  
  ### Import the gridded TMP data 
  rd0 <- brick(paste("D:/Anna/Data/ERA5/tmax/",y,".nc", sep="")) 
  plot(rd0,1)
  rd1 <- rotate(rd0) #convert the LONG from 0:360 to -180:180 degree 
  plot(rd1, 1)
  
  ## Restrict the spatial data to the country boundaries 
  cd0 <- crop(x = rd1, y = bound2)
  cd1 <- rasterize(x = bound2, y = cd0)
  cd2 <- mask(x = cd0, mask = cd1)
  #plot(cd2, 3)

  # Extract mean over each buffer polygon for all layers
  df3 <- raster::extract(cd2, as(buf_15km, "Spatial"),
                         fun = mean, na.rm = TRUE, df = TRUE)
  
  # Add the PSU information
  df3 <- cbind(sp@data, df3)

  # Convert to long format
  df4 <- df3 %>% 
    dplyr::select(-c(ID)) %>% 
    gather(key = date_num, value = tmax, -country)
  
  
  # Add day of the year
  df4 <- df4 %>% 
    mutate(
    doy = as.numeric(sub("X", "", date_num)),   # e.g. "X0" → 0, "X123" → 123
    date = as.Date(doy, origin = paste0(y, "-01-01")),
    year  = format(date, format = "%Y"),
    month = format(date, format = "%m"),
    day   =  format(date, format = "%d")) %>% 
  dplyr::select(-c(date_num, doy))
  
  tmax <- bind_rows(tmax, df4)
  
}

tmin = NULL

for (y in year){
  
  ### Import the gridded TMP data 
  rd0 <- brick(paste("D:/Anna/Data/ERA5/tmin/",y,".nc", sep="")) 
  #plot(rd0,1)
  rd1 <- rotate(rd0) #convert the LONG from 0:360 to -180:180 degree 
  #plot(rd1, 1)
  
  ## Restrict the spatial data to the country boundaries 
  cd0 <- crop(x = rd1, y = bound2)
  cd1 <- rasterize(x = bound2, y = cd0)
  cd2 <- mask(x = cd0, mask = cd1)
  #plot(cd2, 3)
  
  # Extract mean over each buffer polygon for all layers
  df3 <- raster::extract(cd2, as(buf_15km, "Spatial"),
                         fun = mean, na.rm = TRUE, df = TRUE)
  
  # Add the PSU information
  df3 <- cbind(sp@data, df3)

  # Convert to long format
  df4 <- df3 %>% 
    dplyr::select(-c(ID)) %>% 
    gather(key = date_num, value = tmin, -country)
  
  
  # Add day of the year
  df4 <- df4 %>% 
    mutate(
      doy = as.numeric(sub("X", "", date_num)),   # e.g. "X0" → 0, "X123" → 123
      date = as.Date(doy, origin = paste0(y, "-01-01")),
      year  = format(date, format = "%Y"),
      month = format(date, format = "%m"),
      day   =  format(date, format = "%d")) %>% 
    dplyr::select(-c(date_num, doy))
  
  tmin <- bind_rows(tmin, df4)
  
}


tmax <- tmax %>% filter(!is.na(tmax))
tmin <- tmin %>% filter(!is.na(tmin))

### Join the Tmax and Tmin data and calculate Tmean 
temp <- tmax %>% 
  left_join(tmin) %>% 
  mutate(tmean = (tmax+tmin)/2) %>% 
  dplyr::select(country, date, year, month, day, tmax, tmin, tmean)

# convert from kevin to celcius
temp <- temp %>% 
  mutate(tmax = tmax - 273.15,
         tmin = tmin - 273.15,
         tmean = tmean - 273.15) %>% 
  #add diurnal temperature range (DTR)
  mutate(dtr = tmax - tmin)

### Add temperature percentiles
temp <- temp %>% 
  group_by(country) %>% 
  mutate(tmean_pct = ntile(tmean, 100)) %>% 
  mutate(tmax_pct = ntile(tmax, 100)) %>% 
  mutate(tmin_pct = ntile(tmin, 100)) %>% 
  mutate(dtr_pct = ntile(dtr, 100)) %>%
  ungroup()


## Organize in ascending order by year, month, day
temp <- temp %>% 
  arrange(country, year, month, day) %>% 
  dplyr::select(country, date, year, month, day, tmax, tmin, tmean, dtr, tmax_pct, tmin_pct, tmean_pct, dtr_pct)


save(temp, file = "data/PregnAnZI-2 trial/Tmp_data_ERA5.RData")


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

load("data/PregnAnZI-2 trial/apgar_data_for_analysis.RData")

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
  
save(apgar_temp, file = "data/PregnAnZI-2 trial/apgar_data_for_analysis_DBT_ERA5.RData")


