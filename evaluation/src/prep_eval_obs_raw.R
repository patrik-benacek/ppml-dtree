#!/usr/bin/env Rscript

## Prepare CSV dataset from interpolated data
# Usage: 
# ./merge_interp_data.R $leadtime[24/240] $target[t2m/prec24]
# 
# Info:
#   * prec24: square-root transformation is applied

rm(list=ls())

# Input arguments
args = commandArgs(trailingOnly=TRUE)

## load data
DATADIR <- "../data_preproc/data/interp/"
METADIR <- "../observations/data/"
RESDIR  <- "results/"

target   = args[1] 
leadtime = paste0("ff",args[2],"h")

# debugging
#target   = 't2m'
#leadtime = 'ff24h'
#setwd('/home/patrik/Work/czechglobe/TIGGE/evaluation')

# Evaluation period
start_eval <- as.Date("2019-01-01 00:00", tz = "UTC")
end_eval <- as.Date("2019-12-31 00:00", tz = "UTC")

library(ncdf4)
library(tidyverse)

# Read metadata
meta <- read_csv(file.path(METADIR, "metadata_stations.csv"))[,c(1:2)]
colnames(meta) <- c('station_id', 'station')
meta$station = gsub(' / ', '-', meta$station)

print("Read target data ...")
# Target data
nc <- nc_open(file.path(DATADIR, leadtime, paste0("data_target_", target, "_interp.nc")))
  #if (target=='prec24'){
  #  # Square-root transformation
  #  print("--> Square root transformation for preciptation data!!")
  #  fcdata <- sqrt(ncvar_get(nc, paste0(target, "_fc")))
  #  obsdata <- sqrt(ncvar_get(nc, paste0(target, "_obs")))
  #}else{
  fcdata <- ncvar_get(nc, paste0(target, "_fc"))
  obsdata <- ncvar_get(nc, paste0(target, "_obs"))
  #}
  dates <- as.POSIXct(ncvar_get(nc, "time"), origin = "1970-01-01 00:00", tz = "UTC")
  stations <- ncvar_get(nc, "station")
  
  station_metadata <- list()
  station_metadata$altitudes <- ncvar_get(nc, "station_alt")
  station_metadata$latitudes <- ncvar_get(nc, "station_lat")
  station_metadata$longitudes <- ncvar_get(nc, "station_lon")
  station_metadata$locations <- ncvar_get(nc, "station_loc")
nc_close(nc)

dates_vec <- rep(dates, each = length(stations))
stations_vec <- rep(stations, length(dates))
# Prepare observations
obs <- data.frame(dates_vec, 
                  stations_vec,
                  c(obsdata)) 
names(obs) <- c("date", "station_id", "obs")
# Evaluation period only
obs <- subset(obs, date >= start_eval & date <= end_eval)
# Drop nan
obs <- obs[!is.na(obs$obs),] 
# Join metadata
obs <- obs %>% 
  right_join(meta) %>% 
  select(-station_id)
# Save obs
saveRDS(obs, file = file.path(RESDIR, paste0("eval_obs_", target, '_', leadtime, "_2019.RData")))

# Prepare forecast
data <- data.frame()
for (imem in seq(dim(fcdata)[2])){
  #print(paste0("Prepare for member: ", imem))
  d <- data.frame(dates_vec,
                  stations_vec,
                  round(c(fcdata[,imem,]), 4)
                  )
  names(d) <- c("date", "station_id", "fc")
  d$member = imem
  data = rbind(data, d)
}
# Broadcast members to wide format
data_wide <- data %>% 
  spread('member', 'fc') %>% 
  right_join(meta) %>% 
  select(-station_id)

# Evaluation period only
data_wide <- subset(data_wide, date >= start_eval & date <= end_eval)
  
# Save forecast
write.csv(data_wide, file = file.path(RESDIR, 'prediction', paste0("pred_Raw-Fcst_", target, "_", leadtime ,"_2019.csv")), row.names = FALSE)
