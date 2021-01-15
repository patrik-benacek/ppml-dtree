#!/usr/bin/env Rscript
#===============================================
# Deaccumulation script for TIGGE features.
# Arguments:
# * Input: 
#     - file1:   netcdf forecast             [ecmwf_srf_$yyyymmdd.nc]
#     - file2:   previous netcdf forecast    [ecmwf_tp_$yyyymmdd.nc]
#     - fctime:  leadtime of file1           [hours]
#     - accum:   accumulation time           [hours] 
# * Output: 
#     - outfile: deaccumulated forecast      [ecwmf_srf_$yyyymmdd.nc]
#
# Author: Patrik Benacek
#===============================================

rm(list=ls())

# Input arguments
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=5) {
  stop("Input arguments: [forec] [forec_prev] [leadtime [hours]] [accum_time [hours]] [outfile]", call.=FALSE)
}

file1   = args[1]
file2   = args[2]
fctime  = args[3]  
accum   = args[4] 
outfile = args[5] 

#file1    = '/home/patrik/Work/czechglobe/TIGGE/forecasts/data/2015/ecmwf_srf_20150101.nc'
#file2    = '/home/patrik/Work/czechglobe/TIGGE/forecasts/data/2015/ecmwf_tp_20150101.nc'
#outfile  = '/home/patrik/Work/czechglobe/TIGGE/forecasts/data/2015/ecmwf_deaccum_srf_20150101.nc'
#fctime = 240
#accum  = 24

# Variables that should be deaccumulated
attr2deacum =  c('sshf', 'slhf', 'ssr', 'str', 'tp')

# Time initialization
suppressMessages(library(lubridate))
init_date = as.Date(strsplit(strsplit(basename(file1),'_')[[1]][3], '\\.')[[1]][1], format='%Y%m%d')
fc_date = init_date + hours(fctime)
fc_date_prev = fc_date - hours(accum)

if (init_date > fc_date_prev) stop("The accumulation time is longer then the leadtime of forecast!")

if (file.info(file1)$size==0){
    stop(paste("Zero file:", file2, " --> Skip."))
}
if (file.info(file1)$size==0){
    stop(paste("Zero file:", file2, "--> Skip."))
}

# Read data
library(ncdf4)
nc1 <- nc_open(file1)
nc2 <- nc_open(file2)

# Get 
nc_vars1  <- attributes(nc1$var)$names
nc_vars2  <- attributes(nc2$var)$names
nc_fc_dates1 <- as.POSIXct(3600*ncvar_get(nc1, "time"), origin = "1900-01-01 00:00", tz = "UTC")
nc_fc_dates2 <- as.POSIXct(3600*ncvar_get(nc2, "time"), origin = "1900-01-01 00:00", tz = "UTC")

fc1_date_pos = which(nc_fc_dates1==fc_date)
fc2_date_pos = which(nc_fc_dates2==fc_date_prev)

nc_def <- list()
for (i in seq_along(nc_vars1)){
  # Read actual data
  att <- nc_vars1[i]
  data1 = ncvar_get(nc1, att)
  
  # Final dataset initialization
  if (i == 1){
    data <- array(dim = c(length(nc_vars1), dim(data1)))
  }
  data[i,,,,] = data1
  
  # Deaccumulate particular features
  if (att %in% attr2deacum){
    print(paste0("Start deaccumulation for: ", att)) 
    data2 = ncvar_get(nc2, att)
    if (fc1_date_pos>1){
        data_act = data1[,,,fc1_date_pos] 
    }else{
        data_act = data1
    }
    if (fc2_date_pos>1){
        data_prev = data2[,,,fc2_date_pos] 
    }else{
        data_prev = data2 
    }
    data[i,,,,fc1_date_pos] = data_act - data_prev
  }else{
    data[i,,,,] = data1
  }
  nc_def[[i]] <- ncvar_def(name = nc1$var[[att]]$name, units = nc1$var[[att]]$units,                                                                                                                                    
                     dim = nc1$dim, #compression = 1,                                                                                                                        
                     missval = NA, longname = nc1$var[[att]]$longname,                                                                                                                          
                     prec="single") 
}
# Create output files
ncout <- nc_create(outfile,
                   nc_def,                                                                                                      
                   force_v4=T)                                                                                                                                                             

# Write fileds to the output file
for (i in seq_along(nc_vars1)){
  att <- nc_vars1[i]
  ncvar_put(nc = ncout, varid = nc_def[[i]], vals = data[i,,,,])
}

# Close netcdf
nc_close(ncout)  
nc_close(nc1)
nc_close(nc2)

#nc <- nc_open(outfile)
#test <- ncvar_get(nc, 'tp')
#test[1,1,,3]
#nc_close(nc)
