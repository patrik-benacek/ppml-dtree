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
input_data_dir <- "data/interp"
out_data_dir <- "data/merged"

leadtime = paste0("ff",args[1],"h")
target   = args[2] 

surface_attr = c('cape', 'sp', 'sshf', 'slhf', 'msl', 'u10', 'v10', 'd2m', 'ssr', 'str', 'skt', 'sm', 'tcc', 't2m', 'tp')
plev_attr = c('t', 'u', 'v', 'q', 'gh')

# Remove target from surface attributes
target_att <- ifelse(target=='prec24', 'tp', target)
surface_attr <- surface_attr[!surface_attr==target_att]

library(ncdf4)

print("Read target data ...")
# Target data
nc <- nc_open(file.path(input_data_dir, leadtime, paste0("data_target_", target, "_interp.nc")))
  if (target=='prec24'){
    # Square-root transformation
    print("--> Square root transformation for preciptation data!!")
    fcdata <- sqrt(ncvar_get(nc, paste0(target, "_fc")))
    obsdata <- sqrt(ncvar_get(nc, paste0(target, "_obs")))
  }else{
    fcdata <- ncvar_get(nc, paste0(target, "_fc"))
    obsdata <- ncvar_get(nc, paste0(target, "_obs"))
  }
  dates <- as.POSIXct(ncvar_get(nc, "time"), origin = "1970-01-01 00:00", tz = "UTC")
  stations <- ncvar_get(nc, "station")
  
  station_metadata <- list()
  station_metadata$altitudes <- ncvar_get(nc, "station_alt")
  station_metadata$latitudes <- ncvar_get(nc, "station_lat")
  station_metadata$longitudes <- ncvar_get(nc, "station_lon")
  station_metadata$locations <- ncvar_get(nc, "station_loc")
nc_close(nc)

# Covariate data 

# There is same number of surface and pressure level NetCDF files

# surface data [attr, stations, ensemble, date]
print("Read surface data ...")
nc <- nc_open(file.path(input_data_dir, leadtime, paste0("data_aux_srf_", target, "_interp.nc")))
dim_srf <- c(nc$dim$station$len, nc$dim$member$len, nc$dim$time$len)
data_srf <- array(dim=c(length(surface_attr), dim_srf))
date_srf <- as.POSIXct(ncvar_get(nc, "time"), origin = "1970-01-01 00:00", tz = "UTC")

for (i in seq_along(surface_attr)){
  feature = paste0(surface_attr[i], "_fc")
  data_srf[i,,,] <- ncvar_get(nc, feature)
}
nc_close(nc)

# data from pressure level 850 hPa [attr, stations, ensemble, date]
print("Read data at 850hPa pressure level ...")
nc <- nc_open(file.path(input_data_dir, leadtime, paste0("data_aux_pl850_", target, "_interp.nc")))
dim_plv <- c(nc$dim$station$len, nc$dim$member$len, nc$dim$time$len)
data_p500 <- data_p850 <- array(dim=c(length(plev_attr), dim_plv))
date_plv <- as.POSIXct(ncvar_get(nc, "time"), origin = "1970-01-01 00:00", tz = "UTC")
station_plv <- ncvar_get(nc, "station")

# Sanity check
# if (!all(dim_plv == dim_srf)){ stop("Inconsistency between number of surface and pressure NetCDF files.") }

for (i in seq_along(plev_attr)){
  feature = paste0(plev_attr[i], "_pl850_fc")
  data_p850[i,,,] <- ncvar_get(nc, feature)
}
nc_close(nc)

# data from pressure level 500 hPa [attr, stations, ensemble, date]
print("Read data at 500hPa pressure level ...")
nc <- nc_open(file.path(input_data_dir, leadtime, paste0("data_aux_pl500_", target, "_interp.nc")))
for (i in seq_along(plev_attr)){
  feature = paste0(plev_attr[i], "_pl500_fc")
  data_p500[i,,,] <- ncvar_get(nc, feature)
}
nc_close(nc)

# Mask only complete instances i.e. both surface and pressure level data exist
# Pressure level instances exist only if surface exist (see interpolation/interpolation_aux_pl850.R)
# find instances when srf exist and plv not
mask_date_plv <- date_srf %in% date_plv 
mask_station_plv <- stations %in% station_plv
mask_date_srf <- date_plv %in% date_srf 
mask_station_srf <- station_plv %in% stations
# mask fields determined from surface files
data_srf <- data_srf[,mask_station_plv,,mask_date_plv]
fcdata <- fcdata[mask_station_plv,,mask_date_plv]
obsdata <- obsdata[mask_station_plv,mask_date_plv]
dates <- dates[mask_date_plv]
stations <- stations[mask_station_plv]
data_p850 <- data_p850[,mask_station_srf,,mask_date_srf]
data_p500 <- data_p500[,mask_station_srf,,mask_date_srf]

# keep only mean and var
print("Calculate mean and variance features ...")
# observations
obs <- c(obsdata)
rm(obsdata)

# ens fc of target variable
target_mean <- apply(fcdata, c(1,3), mean)
target_var <- apply(fcdata, c(1,3), var)
rm(fcdata)

data_srf_mean <- apply(data_srf, c(1,2,4), mean)
data_srf_var <- apply(data_srf, c(1,2,4), var)
rm(data_srf)

data_p850_mean <- apply(data_p850, c(1,2,4), mean)
data_p850_var <- apply(data_p850, c(1,2,4), var)
rm(data_p850)

data_p500_mean <- apply(data_p500, c(1,2,4), mean)
data_p500_var <- apply(data_p500, c(1,2,4), var)
rm(data_p500)

##
## combine data into a data frame for easier handling 
##

# compute date and station vector
# NOTE: matrices have stations in row and date in column
# c(.) puts together columns into a vector
# here is a trivial example of matrix to vector transformations:
# mm <- matrix(paste0("st",rep(1:3, each = 2)), nrow = 3, byrow = T)
# mm[,1] <- paste0("day1-", mm[,1])
# mm[,2] <- paste0("day2-", mm[,2])
# mm
# c(mm)
# rep(c("day1", "day2"), each = 3)
# rep(paste0("st", 1:3), 2)
# rm(mm)

print("Prepare data frame ...")
dates_vec <- rep(dates, each = length(stations))
stations_vec <- rep(stations, length(dates))

data <- data.frame(dates_vec,
                   stations_vec,
                   obs,
                   c(target_mean),
                   c(target_var) 
                   )

names(data)[1:5] <- c("date", "station", "obs", paste0(target, "_mean"), paste0(target, "_var"))

# Surface features
# Mean
df_ncol <- ncol(data)
for (i in seq_along(surface_attr)){
  name_attr <- surface_attr[i]
  data <- cbind(data, c(data_srf_mean[i,,]))
  names(data)[df_ncol+i] <- paste0(name_attr, "_mean")
}
# Variance
df_ncol <- ncol(data)
for (i in seq_along(surface_attr)){
  name_attr <- surface_attr[i]
  data <- cbind(data, c(data_srf_var[i,,]))
  names(data)[df_ncol+i] <- paste0(name_attr, "_var")
}

# Features at pressure level 850hPa
# Mean
df_ncol <- ncol(data)
for (i in seq_along(plev_attr)){
  name_attr <- plev_attr[i]
  data <- cbind(data, c(data_p850_mean[i,,]))
  names(data)[df_ncol+i] <- paste0(name_attr, "_pl850_mean")
}
# Variance
df_ncol <- ncol(data)
for (i in seq_along(plev_attr)){
  name_attr <- plev_attr[i]
  data <- cbind(data, c(data_p850_var[i,,]))
  names(data)[df_ncol+i] <- paste0(name_attr, "_pl850_var")
}

# Features at pressure level 850hPa
# Mean
df_ncol <- ncol(data)
for (i in seq_along(plev_attr)){
  name_attr <- plev_attr[i]
  data <- cbind(data, c(data_p500_mean[i,,]))
  names(data)[df_ncol+i] <- paste0(name_attr, "_pl500_mean")
}
# Variance
df_ncol <- ncol(data)
for (i in seq_along(plev_attr)){
  name_attr <- plev_attr[i]
  data <- cbind(data, c(data_p500_var[i,,]))
  names(data)[df_ncol+i] <- paste0(name_attr, "_pl500_var")
}
# check output
#head(data)

# save output
#print("Save results as RDS ...")
#save(data, file = file.path(data_dir, "data", "data_all.Rdata"))
print("Save results as CSV ...")
write.csv(data, file = file.path(out_data_dir, paste0("data_", target,"_", leadtime,".csv")), row.names = FALSE)

