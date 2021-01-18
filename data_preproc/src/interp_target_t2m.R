#!/usr/bin/env Rscript

##------------------------------------------------------------------
## Bilinear interpolation of ensemble forecasts to station locations
##------------------------------------------------------------------
rm(list=ls())
data_dir <- "../"

# Input arguments
args = commandArgs(trailingOnly=TRUE)

# Forecast lead times and position in nc files
fc_time = as.integer(args[1])

years = c(2015,2019)
years_regex = "201[5,6,7,8,9]"

# time position in NetCDF 1/3
if (fc_time==24) fc_time_pos = 1
if (fc_time==240) fc_time_pos = 3

print("Read data")
suppressMessages(library(lubridate))
suppressMessages(library(readr))
suppressMessages(library(ncdf4))  ## load nc files
suppressMessages(library(akima))  ## bilinear interpolation to station locations
#?bilinear
suppressMessages(library(abind))  ## append vectors / matrices / arrays 

#-------------------------
# Read observations 
#-------------------------
suffix = "srf"
if (fc_time > 24){
    print("Load deaccumulated surface data: ecmwf_srf_deacc_${yyyymmddnt}.nc")
    suffix = "srf_deacc"
}

# Use deaccumulated data if neccessary
fc_files <- Sys.glob(file.path(data_dir, "forecasts/data", years_regex, paste0("ecmwf_", suffix, "_", years_regex,"*.nc")))
if (fc_time>24) warnings(paste("Check if deaccumulated fileds are included in", fc_files[1]))

# Forecast NetCDF field initialization
nc <- nc_open(fc_files[1])
fc_raw_lat <- ncvar_get(nc, "latitude")                                                                                                                                                                   
fc_raw_lon <- ncvar_get(nc, "longitude")                                                                                                                                                                  
fc_raw_member <- ncvar_get(nc, "number")                                                                                                                                                                  
fc_raw_time_unconverted <- ncvar_get(nc, "time")                                                                                                                                                          
fc_raw_validtime0 <- as.POSIXct(3600*fc_raw_time_unconverted[fc_time_pos], origin = "1900-01-01 00:00", tz = "UTC")
rm(fc_raw_time_unconverted)   
nc_close(nc)

fc_time_hour <- hour(fc_raw_validtime0)

print("Read observations ...")

# Read station metadata
st_meta <- read_csv(file.path(data_dir, 'observations', 'data', 'stations_cz.csv'))
#str(st_meta)

# Read station data
st_temp <- data.frame()
for (stat_id in st_meta$wmo_id){
  filestats = Sys.glob(file.path(data_dir, 'observations/data/', years_regex, paste0("obs_stat_", stat_id, "_*.rds")))
  for (filestat in filestats){
    if (file.exists(filestat)){
      st_temp_ <- readRDS(filestat)
      # Subset of observations provided the initialization time is 00UTC
      st_temp_ = st_temp_[hour(st_temp_$Date)==fc_time_hour & minute(st_temp_$Date)==0,]
      st_temp <- rbind(st_temp, st_temp_) 
    }else{
      print(paste0("Station", stat_id, "does not exist."))
    }
  }
}

st_validtime <- st_temp$Date
st_stations <- st_temp$station_ID
stationlist <- st_meta$wmo_id

# array with interpolated forecasts [stationID, member, validtime]                                                                                                                                        
fc_interpolated_temp <- array(dim = c(length(st_meta$wmo_id),                                                                                                                                             
                                      length(fc_raw_member),                                                                                                                                              
                                      length(fc_raw_validtime0)))                                                                                        

# array with observed values [stationID, validtime]                                                                                                                                                       
obs_temp <- array(dim = c(length(st_meta$wmo_id),                                                                                                                                                         
                          length(fc_raw_validtime0)))  

# array with validtimes [validtime]                                                                                                                                                       
fc_raw_validtime <- array(dim = c(length(fc_raw_validtime0)))

print("Read forecast and interpolate to station locations:")

for(filename in fc_files){

  print(paste("-->", basename(filename)))
  
  # Check file size
  if (file.info(filename)$size == 0){ print("----> Not available."); next }

  # load forecast from NetCDF file
  nc <- nc_open(filename)
  
  # read forecast attributes
  fc_raw_temp <- ncvar_get(nc, "t2m")
  fc_raw_lat <- ncvar_get(nc, "latitude")                                                                                                                                                                   
  fc_raw_lon <- ncvar_get(nc, "longitude")                                                                                                                                                                  
  fc_raw_member <- ncvar_get(nc, "number")
  fc_raw_time_unconverted <- ncvar_get(nc, "time")
  fc_raw_validtime_append <- as.POSIXct(3600*fc_raw_time_unconverted[fc_time_pos], origin = "1900-01-01 00:00", tz = "UTC")
  rm(fc_raw_time_unconverted)
  
  nc_close(nc)
  
  # array with interpolated forecasts [stationID, member, validtime], to be appended later
  fc_interpolated_temp_append <- array(dim = c(length(st_meta$wmo_id),
                                               length(fc_raw_member),
                                               length(fc_raw_validtime_append)))
  
  # array with observed values [stationID, validtime], to be appended later
  obs_temp_append <- array(dim = c(length(st_meta$wmo_id),
                                   length(fc_raw_validtime_append)))
  
  for(thisst in stationlist){
    
    thisst_pos <- which(stationlist == thisst)

    if(thisst_pos %% 50 == 0){
      cat("station", thisst_pos, "of", length(stationlist), "starting at", paste(Sys.time()),"\n"); flush(stdout())
    }
    
    for(vtime in fc_raw_validtime_append){
      # last observation is recorded prior to last forecast valid time
      #if(vtime > "2015-12-31"){
      #  next
      #}
      vtime_pos <- which(fc_raw_validtime_append == vtime)
      
      st_temp_pos <- which(st_validtime == vtime & st_stations == thisst)
      # st_temp contains no NA values, in case of missing observations, st_temp_pos is "integer(0)"
      if(length(st_temp_pos) == 0){
        obs_temp_append[thisst_pos,vtime_pos] <- NA
      # There are duplicities in observations
      } else if(length(st_temp_pos)>1){
        st_temp_pos = st_temp_pos[1]
        obs_temp_append[thisst_pos,vtime_pos] <- st_temp$TC[st_temp_pos]
      }else{
        obs_temp_append[thisst_pos,vtime_pos] <- st_temp$TC[st_temp_pos]
      }
      
      #----------------------------------------------------------------
      # Bilinear interpolation of nc forecasts to station locations
      #----------------------------------------------------------------
      for(mem in 1:50){
        thisst_lat <- st_meta$lat[which(st_meta$wmo_id == thisst)]
        thisst_lon <- st_meta$lon[which(st_meta$wmo_id == thisst)]
        
        # find latitude and longitude of surrounding grib boxes
        lat_low <- max(fc_raw_lat[which(fc_raw_lat <= thisst_lat)])
        lat_high <- lat_low + 0.5
        lon_low <- max(fc_raw_lon[which(fc_raw_lon <= thisst_lon)])
        lon_high <- lon_low + 0.5
        x_coord <- c(lon_low, lon_high)
        y_coord <- c(lat_low, lat_high)
        
        # find corresponding positions in nc dimensions
        lat_pos <- c(which(fc_raw_lat == lat_low),which(fc_raw_lat == lat_high))
        lon_pos <- c(which(fc_raw_lon == lon_low),which(fc_raw_lon == lon_high))
        
        # bilinearly interpolate forecasts
        tmp <- bilinear(x = x_coord, y = y_coord, 
                        z = fc_raw_temp[lon_pos, lat_pos, mem, fc_time_pos],
                        x0 = thisst_lon, y0 = thisst_lat)$z      
        # ... and convert to degrees Celsius  
        fc_interpolated_temp_append[thisst_pos, mem, vtime_pos] <- tmp - 273.15
      }
    }
  }
  
  if (which(fc_files==filename)==1){
    # append new objects to first validate
    obs_temp[,1] <- obs_temp_append
    fc_interpolated_temp[,,1] <- fc_interpolated_temp_append
    fc_raw_validtime[1] <- fc_raw_validtime_append

  }else{
    # append new objects to existing ones
    obs_temp <- abind(obs_temp, 
                      obs_temp_append, 
                      along = 2)
    fc_interpolated_temp <- abind(fc_interpolated_temp,
                                  fc_interpolated_temp_append,
                                  along = 3)
    fc_raw_validtime <- c(fc_raw_validtime, fc_raw_validtime_append)
  }
  rm(obs_temp_append, fc_interpolated_temp_append, fc_raw_validtime_append)
}

#------------------------
# Write to netCDF file
#------------------------
print("Write to netCDF file ...")

# define dimensions
stationdim <- ncdim_def("station", "station_ID", as.integer(stationlist))
memberdim <- ncdim_def("member", "member_number", as.integer(1:50))
timedim <- ncdim_def(name = "time", vals = as.integer(fc_raw_validtime),
                     units = "seconds since 1970-01-01 00:00 UTC",
                     longname = "valid time of forecasts and observations, UTC")
# as.POSIXct(as.integer(fc_raw_validtime), tz = "UTC", origin = "1970-01-01 00:00")

# define variables
fillvalue <- NA

dlname <- "interpolated t2m ensemble forecast"
t2m_fc_def <- ncvar_def(name = "t2m_fc", units = "deg_C",
                        dim = list(stationdim,memberdim,timedim),
                        missval = fillvalue, longname = dlname,
                        prec="single")

dlname <- "t2m station observation"
t2m_obs_def <- ncvar_def(name = "t2m_obs", units = "deg_C",
                         dim = list(stationdim,timedim),
                         missval = fillvalue, longname = dlname,
                         prec="single")

dlname <- "altitude of station"
alt_def <- ncvar_def(name = "station_alt", units = "m",
                     dim = list(stationdim),
                     missval = fillvalue, longname = dlname,
                     prec = "single")

dlname <- "latitude of station"
lat_def <- ncvar_def(name = "station_lat", units = "degrees north",
                     dim = list(stationdim),
                     missval = fillvalue, longname = dlname,
                     prec = "single")

dlname <- "longitude of station"
lon_def <- ncvar_def(name = "station_lon", units = "degrees east",
                     dim = list(stationdim),
                     missval = fillvalue, longname = dlname,
                     prec = "single")

dlname <- "station ID"
id_def <- ncvar_def(name = "station_id", units = "",
                     dim = list(stationdim),
                     missval = fillvalue, longname = dlname,
                     prec = "single")

# character variables (location names) require special attention
dimnchar <- ncdim_def("nchar", "", 1:stationdim$len, create_dimvar=FALSE,
                      longname = "number of characters for locations")
location_def <- ncvar_def("station_loc", "", list(dimnchar,stationdim),
                          prec = "char", longname = "location of station")

## create nc file
ncfile_name <- file.path("data/interp", paste0("ff", fc_time, "h"), paste0("data_target_t2m_interp.nc"))
ncout <- nc_create(ncfile_name,
                   list(t2m_fc_def, t2m_obs_def, alt_def, lat_def, lon_def, id_def, location_def),
                   force_v4=T)

# put variables
ncvar_put(nc = ncout, varid = t2m_fc_def, vals = fc_interpolated_temp)
ncvar_put(nc = ncout, varid = t2m_obs_def, vals = obs_temp)
ncvar_put(nc = ncout, varid = alt_def, vals = st_meta$alt)
ncvar_put(nc = ncout, varid = lat_def, vals = st_meta$lat)
ncvar_put(nc = ncout, varid = lon_def, vals = st_meta$lon)
ncvar_put(nc = ncout, varid = id_def, vals = st_meta$wmo_id)
ncvar_put(nc = ncout, varid = location_def, vals = st_meta$station_names)

nc_close(ncout)

print("Finish successfully.")
