#!/usr/bin/env Rscript

#-----------------------------------------------------------------------
# Extract model orography in station location and save as CSV file
#-----------------------------------------------------------------------

rm(list=ls())

data_dir <- "/home/patrik/Work/czechglobe/TIGGE/data_preproc/generate_dataset/"
input_data_dir <- "/home/patrik/Work/czechglobe/TIGGE/data_preproc/interpolation/data"
metadata_id <- read.csv(file.path(data_dir, "data", "metadata_station.csv"))

leadtime="ff240h"
target='t2m'  # t2m/prec24

library(ncdf4)
nc <- nc_open(file.path(input_data_dir, leadtime, "data_aux_geo_interp.nc"))
metadata <- list()
metadata$station_id <- ncvar_get(nc, "station_id")
metadata$orog <- round(ncvar_get(nc, "orog"), 0)
nc_close(nc)

metadata <- as.data.frame(metadata)

# merge station and model metadata
metadata_id <- merge(metadata_id, metadata, by.x = "wmo_id", by.y = "station_id")

write.csv(metadata_id, file = file.path(data_dir, "data", "metadata_all.csv"), row.names = FALSE)

data <- read.csv(file.path(data_dir, "data", paste0("data_", target,"_", leadtime, "_2015_2019.csv")))
data$date = as.Date(data$date)

metadata_names <- c("date", "station", "station_names", "lon", "lat", "alt", "orog")
data_names <- colnames(data[,-which(colnames(data) %in% metadata_names)])

# Merge forecast/observation data with station/model metadata
data <- merge(data, metadata_id, by.x = "station", by.y = "wmo_id")[,c(metadata_names, data_names)]

# Sort data according to date
data <- data[order(data$date),]

# Write data
write.csv(data, file = file.path(data_dir, "data", paste0("data_wmeta_", target,"_", leadtime, "_2015_2019.csv")), row.names = FALSE)
