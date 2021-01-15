#!/usr/bin/env Rscript

## bilinear interpolation of ensemble forecasts to station locations   
rm(list=ls()) 
path <- "/home/patrik/Work/czechglobe/TIGGE/observations/"

# Read SYNOP observation from OGIMET
library(climate)
library(readr)

# Get list of CZ station
stations_cz = stations_ogimet(country ="Czech+Republic", add_map = TRUE)
# Save list of station
write_csv(stations_cz, file.path(path, "data", "stations_cz.csv"))

lstat_ids = as.character(stations_cz$wmo_id)
  
# Downloading
years = c(2015:2019)
for (year in years){
    cyear = as.character(year)
    start = paste0(cyear,"-01-01")
    ende  = paste0(cyear,"-12-31")

    # Create directory where to save data
    dir.create(file.path(path, cyear), showWarnings = FALSE)

    for (stat_id in lstat_ids){

        print(paste0("Read data for ", cyear, ": ", stat_id))
        file = file.path(path, "data", cyear, paste0("obs_stat_", stat_id, "_", cyear,".rds"))

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
}
