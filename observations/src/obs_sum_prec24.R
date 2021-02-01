#!/usr/bin/env Rscript
#---------------------------------------------------
# Prepare sum of 24-hour precipitation valid at 00utc
# Author: Patrik Benacek
#---------------------------------------------------
rm(list=ls()) 

library(tidyverse)
library(lubridate)
library(data.table)

# Input arguments
argv = commandArgs(trailingOnly=TRUE)

# Specify obs period
year = argv[1]

# Start script
files = list.files(path=file.path("data", year), pattern='^obs_stat_')
for (file in files){
  print(paste("Prepare data for", file))
  # Read ogimet data
  data = readRDS(file.path("data", year, file))
  # The 24h precipitation are measured at 06UTC. Therefore, we
  # shift them to 00UTC as model forecast.
  data$prec24 = shift(data$pr24, 6)
  # Transform to daily data
  data = data[hour(data$Date)==0,]
  # Save prec24h observations
  saveRDS(data[,c("station_ID", "Date", "prec24")], file = file.path("data", year, paste0('prec24_', file)))
}

print("Finish.")