#!/usr/bin/env Rscript

##------------------------------------------------------------------
## Bilinear interpolation of ensemble forecasts to station locations
##------------------------------------------------------------------
rm(list=ls())
data_dir <- "/home/patrik/Work/czechglobe/TIGGE/"

years = c(2015)
years_regex = "2015"

# Forecast lead times and position in nc files
fc_time = 24
fc_time_pos = 1
fc_time_init = 0

print("Read data")

attributes = c('orog')

#-------------------------
# Read observations 
#-------------------------
suppressMessages(library(readr))

print("Read observations ...")

# Read station metadata
st_meta <- read_csv(file.path(data_dir, 'observations', 'stations_cz.csv'))
#str(st_meta)

#----------------------------------------------------------------
# Bilinear interpolation of nc forecasts to station locations
#----------------------------------------------------------------
suppressMessages(library(ncdf4))  ## load nc files
suppressMessages(library(akima))  ## bilinear interpolation to station locations
#?bilinear
suppressMessages(library(abind))  ## append vectors / matrices / arrays 

stationlist <- st_meta$wmo_id

fc_files <- Sys.glob(file.path(data_dir, "forecasts/data/ecmwf_geo_*.nc"))

# Forecast NetCDF field initialization
nc <- nc_open(fc_files[1])
fc_raw_lat <- ncvar_get(nc, "latitude")                                                                                                                                                                   
fc_raw_lon <- ncvar_get(nc, "longitude")                                                                                                                                                                  
fc_raw_member <- 1
fc_raw_time_unconverted <- ncvar_get(nc, "time")                                                                                                                                                          
fc_raw_validtime0 <- as.POSIXct(3600*fc_raw_time_unconverted[fc_time_pos], origin = "1900-01-01 00:00", tz = "UTC")
rm(fc_raw_time_unconverted)   
nc_close(nc)

# array with interpolated forecasts [stationID, member, validtime, attribute]
fc_interpolated <- array(dim = c(length(st_meta$wmo_id),
                                 length(fc_raw_member),
                                 length(fc_raw_validtime0),
                                 length(attributes)
                                 ))

# Cycle along auxilliary attributes
for (i in seq_along(attributes)){

    attr = attributes[i]

    # array with validtimes [validtime]
    fc_raw_validtime <- array(dim = c(length(fc_raw_validtime0)))
    
    # array with interpolated forecasts [stationID, member, validtime]
    fc_interpolated_attr <- array(dim = c(length(st_meta$wmo_id),
                                          length(fc_raw_member),
                                          length(fc_raw_validtime0)
    ))

    print("Read forecast and interpolate to station locations:")

    for(j in seq_along(fc_files)){

        filename = fc_files[j]

        print(paste("-->", basename(filename)))
  
        # load forecast from NetCDF file
        nc <- nc_open(filename)
  
        # read forecast attributes
        fc_raw_att <- ncvar_get(nc, attr)
        fc_raw_lat <- ncvar_get(nc, "latitude")
        fc_raw_lon <- ncvar_get(nc, "longitude")

        fc_raw_member <- 1
        fc_raw_time_unconverted <- ncvar_get(nc, "time")
        fc_raw_validtime_append <- as.POSIXct(3600*fc_raw_time_unconverted[fc_time_pos], origin = "1900-01-01 00:00", tz = "UTC")
        rm(fc_raw_time_unconverted)
  
        nc_close(nc)
  
        # array with interpolated forecasts [stationID, member, validtime], to be appended later
        fc_interpolated_append <- array(dim = c(length(st_meta$wmo_id),
                                                length(fc_raw_validtime_append)))
  
        for(thisst in stationlist){
    
            thisst_pos <- which(stationlist == thisst)
            thisst_lat <- st_meta$lat[which(st_meta$wmo_id == thisst)]
            thisst_lon <- st_meta$lon[which(st_meta$wmo_id == thisst)]
    
            # find latitude and longitude of surrounding grib boxes
            lat_pos_mindiff <- order(abs(fc_raw_lat - thisst_lat))[1]
            lon_pos_mindiff <- order(abs(fc_raw_lon - thisst_lon))[1]
    
            if(thisst_pos %% 50 == 0){
                cat("station", thisst_pos, "of", length(stationlist), "starting at", paste(Sys.time()),"\n"); flush(stdout())
            }
    
            for(vtime in fc_raw_validtime_append){

                vtime_pos <- which(fc_raw_validtime_append == vtime)
      
                # choose forecasts from nearest grid point
                fc_interpolated_append[thisst_pos, vtime_pos] <- fc_raw_att[lon_pos_mindiff, lat_pos_mindiff, vtime_pos]
            }
        }
  
        if (which(fc_files==filename)==1){
            # append new objects to first validate
            fc_interpolated_attr[,1,1] <- fc_interpolated_append
            fc_raw_validtime[1] <- fc_raw_validtime_append

        }else{
            # append new objects to existing ones
            fc_interpolated_attr  <- abind(fc_interpolated_attr,
                                           fc_interpolated_append,
                                           along = 3)
            fc_raw_validtime <- c(fc_raw_validtime, fc_raw_validtime_append)
        }
        rm(fc_interpolated_append)
    } # loop filenames
    
    # Reshape fc_interpolated_attr [stationID, member, validtime, attribute]
    fc_interpolated_attr = array(fc_interpolated_attr, , dim=c(dim(fc_interpolated_attr), 1))
    
    if (i==1){
        # append new objects
        fc_interpolated <- fc_interpolated_attr
    }else{
        # append new objects to existing ones
        fc_interpolated  <- abind(fc_interpolated,
                                  fc_interpolated_attr,
                                  along = 4)
    }
    rm(fc_interpolated_attr)
    
} # loop attributes
rm(fc_raw_att)

#------------------------
# Write to netCDF file
#------------------------
print("Write to netCDF file ...")

# define dimensions
stationdim <- ncdim_def("station", "station_ID", as.integer(stationlist))
timedim <- ncdim_def(name = "time", vals = as.integer(fc_raw_validtime),
                     units = "seconds since 1970-01-01 00:00 UTC",
                     longname = "valid time of forecasts and observations, UTC")
# as.POSIXct(as.integer(fc_raw_validtime), tz = "UTC", origin = "1970-01-01 00:00")

source(file.path(data_dir, "data_preproc/interpolation", "config_aux.R"))

# define variables
fillvalue <- NA
fc_def <- list()

for (i in seq_along(attributes)){

    attr = attributes[i]
    attr_info <- get_attribute_info(attr)
    sname = attr_info[1]
    lname = attr_info[2]
    units = attr_info[3]

    fc_def[[i]] <- ncvar_def(name = sname, units = units,
                             dim = list(stationdim,timedim),
                             missval = fillvalue, longname = lname,
                             prec="single")
}

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

meta_def <- list(alt_def, lat_def, lon_def, id_def, location_def)

## create nc file
ncfile_name <- file.path(data_dir, "data_preproc", "interpolation", "data", paste0("data_aux_geo_interp.nc"))
ncout <- nc_create(ncfile_name, 
                   c(fc_def, meta_def), 
                   force_v4=T)

# Put variables (fc_interpolated[statID, members=1, validtime, attribute])
for (i in seq_along(attributes)){
  ncvar_put(nc = ncout, varid = fc_def[[i]], vals = fc_interpolated[,1,,i])
}
ncvar_put(nc = ncout, varid = alt_def, vals = st_meta$alt)
ncvar_put(nc = ncout, varid = lat_def, vals = st_meta$lat)
ncvar_put(nc = ncout, varid = lon_def, vals = st_meta$lon)
ncvar_put(nc = ncout, varid = id_def, vals = st_meta$wmo_id)
ncvar_put(nc = ncout, varid = location_def, vals = st_meta$station_names)

nc_close(ncout)

print("Finish successfully.")
