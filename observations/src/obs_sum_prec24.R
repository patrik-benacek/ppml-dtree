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
args = commandArgs(trailingOnly=TRUE)

# Specify obs period
year = argv[1]

preproc_precip <- function(df){
  #
  # Prepare the 24h sum of precipitation at 00utc for input OGIMET data.
  #
  # Data subset
  prec <- df[which(hour(df$Date)%in%c(0,6,12,18)),c('Date','pr6','pr12','pr24')]
  # Fill sum of 6h precipitation
  prec$lg_pr6 = shift(prec$pr6, -1)
  prec$prr6 <- prec$pr12 - prec$lg_pr6
  prec[is.na(prec$pr6), 'pr6'] <- prec[is.na(prec$pr6), 'prr6']
  prec$prr6 <- NULL
  prec$lg_pr6 <- NULL
  
  # Sum precipitation by 6, 12, 18 and 00
  # Prepare sumday
  prec$sumday = ymd_hms(prec$Date) - hours(6)
  # Group by sumday
  totprec <- prec %>% 
    mutate(group_date = as.Date(sumday)) %>% 
    group_by(group_date) %>% 
    summarise(prec24 = sum(pr6, na.rm=TRUE)) %>% 
    rename(Date = group_date)
  # Associate to proper day
  #totprec$Date <- ymd(totprec$Date) + days(1)
  
  return(totprec)
}

# Start script
files = list.files(path=file.path("data", year), pattern='^obs_stat_')
for (file in files){
  print(paste("Prepare data for", file))
  # Read ogimet data
  data = readRDS(file.path("data", year, file))
  out = preproc_precip(data)
  # Save prec24h observations
  saveRDS(out, file = file.path("data", year, paste0('prec24_', file)))
}

print("Finish.")
