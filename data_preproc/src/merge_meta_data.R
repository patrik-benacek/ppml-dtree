#!/usr/bin/env Rscript

#-----------------------------------------------------------------------
# Add model and station metadata to the merged dataset
# Inputs:
#   * station-meta: observations/data/metadata_stations.csv
#   * model-meta:   data/interp/data_aux_geo_interp.nc      
#   * interp-data:  data/interp/$leadtime/data_$target_$leadtime.csv
# Output:
#   * dataset:
#-----------------------------------------------------------------------

rm(list=ls())

# Input arguments
args = commandArgs(trailingOnly=TRUE)

# Data location
input_data_dir <- "data/interp"
input_obs_dir <- "../observations/data"
out_data_dir <- "data/merged"

leadtime = paste0("ff",args[1],"h")
target   = args[2] 

# Read station metadata
metadata_id <- read.csv(file.path(input_obs_dir, "metadata_stations.csv"))

library(ncdf4)
# Read model orography metadata
nc <- nc_open(file.path(input_data_dir, "data_aux_geo_interp.nc"))
metadata <- list()
metadata$station_id <- ncvar_get(nc, "station_id")
metadata$orog <- round(ncvar_get(nc, "orog"), 0)
nc_close(nc)
metadata <- as.data.frame(metadata)

# merge station and model metadata
metadata_id <- merge(metadata_id, metadata, by.x = "wmo_id", by.y = "station_id")

# Output station+model metadata file
#write.csv(metadata_id, file = file.path(data_dir, "data", "metadata_all.csv"), row.names = FALSE)

data <- read.csv(file.path(out_data_dir, leadtime, paste0("data_", target,"_", leadtime, ".csv")))
data$date = as.Date(data$date)

metadata_names <- c("date", "station", "station_names", "lon", "lat", "alt", "orog")
data_names <- colnames(data[,-which(colnames(data) %in% metadata_names)])

# Merge forecast/observation data with station/model metadata
data <- merge(data, metadata_id, by.x = "station", by.y = "wmo_id")[,c(metadata_names, data_names)]

# Sort data according to date
data <- data[order(data$date),]

# Write data
write.csv(data, file = file.path(out_data_dir, leadtime, paste0("data_wmeta_", target,"_", leadtime, ".csv")), row.names = FALSE)