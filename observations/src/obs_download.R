#!/usr/bin/env Rscript

# Download station observations from OGIMET by climate package:
# Usage: ./obs_download.R $station_country $year
#
# Variables:
#    station_ID 	- WMO station identifier
#    Lon 			- longitude
#    Lat 			- latitude
#    Date 			- date (and time) of observations
#    TC 			- air temperature at 2 metres above ground level. Values given in Celsius degrees
#    TdC 			- dew point temperature at 2 metres above ground level. Values given in Celsius degrees
#    TmaxC 			- maximum air temperature at 2 metres above ground level. Values given in Celsius degrees
#    TminC 			- minimum air temperature at 2 metres above ground level. Values given in Celsius degrees
#    ddd 			- wind direction
#    ffkmh 			- wind speed in km/h
#    Gustkmh 		- wind gust in km/h
#    P0hpa 			- air pressure at elevation of the station in hPa
#    PseahPa 		- sea level pressure in hPa
#    PTnd 			- pressure tendency in hPa
#    Nt 			- total cloud cover
#    Nh 			- cloud cover by high-level cloud fraction
#    HKm 			- height of cloud base
#    InsoD1			- insolation in hours
#    Viskm 			- visibility in kilometres
#    Snowcm			- depth of snow cover in centimetres
#    pr6 			- precicipitation totals in 6 hours
#    pr12 			- precicipitation totals in 12 hours
#    pr24 			- precicipitation totals in 24 hours
#    TemperatureCAvg - average air temperature at 2 metres above ground level. Values given in Celsius degrees
#    TemperatureCMax - maximum air temperature at 2 metres above ground level. Values given in Celsius degrees
#    TemperatureCMin - minimum air temperature at 2 metres above ground level. Values given in Celsius degrees
#    TdAvgC 		- average dew point temperature at 2 metres above ground level. Values given in Celsius degrees
#    HrAvg 			- average relative humidity. Values given in %
#    WindkmhDir 	- wind direction
#    WindkmhInt 	- wind speed in km/h
#    WindkmhGust 	- wind gust in km/h
#    PresslevHp 	- Sea level pressure in hPa
#    Precmm 		- precipitation totals in mm
#    TotClOct 		- total cloudiness in octants
#    lowClOct 		- cloudiness by low level clouds in octants
#    SunD1h 		- sunshine duration in hours
#    PreselevH 		- atmospheric pressure measured at altitude of station in hPa
#    SnowDepcm 		- depth of snow cover in centimetres

rm(list=ls()) 

# Input arguments
args = commandArgs(trailingOnly=TRUE)

# Read SYNOP observation from OGIMET
library(climate)
library(readr)

# Input arguments
station_country = argv[1] 
cyear = argv[2]

start = paste0(cyear,"-01-01")
ende  = paste0(cyear,"-12-31")

# Get list of CZ station
stations = stations_ogimet(country = station_country, add_map = TRUE)
# Save list of station
write_csv(stations, file.path("data", "metadata_stations.csv"))
lstat_ids = as.character(stations$wmo_id)

# Create directory where to save data
dir.create(file.path("data", cyear), showWarnings = FALSE)

# Start downloading stations
for (stat_id in lstat_ids){

    print(paste0("Read data for ", cyear, ": ", stat_id))
    file = file.path("data", cyear, paste0("obs_stat_", stat_id, "_", cyear,".rds"))

    skip_to_next <- FALSE
    if (!file.exists(file)){
        # Fetch observations
        tryCatch(
        expr = {
            obs_data <- meteo_ogimet(interval = "hourly", date = c(start, ende), station = stat_id);
            # Save observations
            saveRDS(obs_data, file = file)
            },
            error = function(e){
            message("Problem with fetching data.")     
            skip_to_next <<- TRUE 
            }
        )
    }else{
        print(paste("File:", file, "exists. Skip..."))
    }
    if (skip_to_next){ next }
}