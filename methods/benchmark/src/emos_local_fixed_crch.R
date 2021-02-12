#!/usr/bin/env Rscript                                                                                                                                                                       

##---------------------------------
## Run selected benchmark models
##---------------------------------
## * EMOS local, fixed window
## * chosen tuning parameters: none
## Precipitation case:
## square root transformation 
## we learn regression model in square-root space: x=4 --> predict_rain(sqrt(4))^2
## how use crch package: https://eeecon.uibk.ac.at/~zeileis/papers/KIT-2016.pdf

rm(list=ls())

pkgTest <- function(x)
  {
    if (!require(x,character.only = TRUE))
    {
      install.packages(x,dep=TRUE)
        if(!require(x,character.only = TRUE)) stop("Package not found")
    }
  }

# Load packages
suppressMessages(pkgTest("lubridate"))
suppressMessages(pkgTest("scoringRules"))
suppressMessages(pkgTest("crch"))
suppressMessages(pkgTest("readr"))

# Input arguments
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)<3) {
  stop("Input arguments: [target variable [t2m/prec24]] [start_train_period [year]] [forecast leadtime [hours]]", call.=FALSE)
}

target = args[1]
leadtime_hours = args[2]
start_train_year = args[3]

# Settings
exp_name <- "EMOS-local"
leadtime <- paste0("ff", leadtime_hours,"h")
train_date <- c(paste0(start_train_year, "-01-01"),"2019-01-01")
eval_date <- c("2019-01-01","2019-12-31")

# Summary
print("===========================================================")
print(paste0("Run benchmark experiment: ", exp_name))
print("===========================================================")
print("-------------------")
print("Setting")
print("-------------------")
print(paste0("Target: ", target))
print(paste0("Leadtime: ", leadtime_hours, "h"))
print(paste0("Train period: ", train_date[1], " to ", train_date[2]))
print(paste0("Test period : ", eval_date[1], " to ", eval_date[2]))

in_data_dir <- "../../data_preproc/data/processed"
out_data_dir <- "results"

source("src/config.R")

# Name of experiment
#exp_suffix = paste(target, leadtime, year(ymd(train_date[1])), year(ymd(train_date[2])), sep='_')
exp_suffix = paste(target, leadtime, year(ymd(train_date[1])), sep='_')

# Read station measurements
file_obs = paste0("data_", target, "_", leadtime, ".zip")
data <- read_csv(file.path(in_data_dir, file_obs))
data$date = as.Date(data$date)
data$station_names = gsub(' / ', '-', data$station_names)

print("")
print("-------------------")
print("Data:") 
print("-------------------")
print(paste("--> samples:", dim(data)[1]))
print(paste("--> columns:", dim(data)[2]))
print(paste("Num. of metadata:", ncol(data[,c(1:7)])))
print(paste("Num. of features:", ncol(data[,c(8:ncol(data))])))

data <- data[, -which(!(names(data) %in% c("obs", "date", "station_names", paste0(target,"_mean"), paste0(target,"_var"))))]
names(data)[4:5] <- c('fc_mean', 'fc_var')
  
# Data cleaning: square-root scaled data (negative->NAN), remove NAN values, omit 'perfect' ensemble prediction
data <- data[data$fc_var>0,]
data <- data[!is.na(data$obs),]
data <- data[!is.na(data$fc_mean),]
data$fc_std <- sqrt(data$fc_var)

train_start <- as.Date(paste(train_date[1], "00:00 UTC")) 
train_end <- as.Date(paste(train_date[2], "00:00 UTC")) - days(2)

data_train_all <- subset(data, date >= train_start & date <= train_end)
data_train_all <- data_train_all[complete.cases(data_train_all), ]

eval_start <- as.Date(paste(eval_date[1], "00:00 UTC"))
eval_end <- as.Date(paste(eval_date[2], "00:00 UTC"))
eval_dates <- seq(eval_start, eval_end, by = "1 day")

data_eval_all <- subset(data, date >= eval_start & date <= eval_end)

print("")
print("-------------------")
print("Forecast:") 
print("-------------------")
t1 <- Sys.time()
# List of stations
stations_list <- unique(data$station_names)
out_loc <- rep(NA, nrow(data_eval_all))
out_sc <- rep(NA, nrow(data_eval_all))
skip_stations <- c()
#par_out <- matrix(NA, ncol = 4, nrow = length(stations_list))

for(this_station in stations_list){
  
  data_train <- subset(data_train_all, station_names == this_station)
  
  # remove incomplete cases (= NA obs or fc)
  data_train <- data_train[complete.cases(data_train), ]
  
  # skip station if there are:
  # i) too few forecast cases (limit 10)
  # i) more than 75% singular values (e.g. zero precipitation) --> problem with hessian function
  if(nrow(data_train) < 10 | ((sum(data_train$obs==0)/nrow(data_train))>.75)){
    skip_stations <- append(skip_stations, this_station)
    print(paste0("Station:", this_station, ", zero:", sum(data_train$obs==0), " --> skip EMOS"))
    next
  }else{
    print(paste0("Station:", this_station, ", zero:", sum(data_train$obs==0)))
  }
  
  # determine optimal EMOS coefficients a,b,c,d using minimum CRPS estimation
  emos_crps <- get_crch_model(train=data_train, target=target)
  data_eval <- subset(data_eval_all, station_names == this_station)
  
  # save parameters
  ind_this_station <- which(data_eval_all$station_names == this_station)
  
  out_loc[ind_this_station] <- predict(emos_crps, type='location', newdata = data_eval[,c('fc_mean', 'fc_std')])
  out_sc[ind_this_station]  <- predict(emos_crps, type='scale'   , newdata = data_eval[,c('fc_mean', 'fc_std')])
}
t2 <- Sys.time()
print(paste0("Computational time: ", round(t2-t1, 1), "s"))

# Remove stations with singular data
if (target=='prec24'){
  # Make prediction
  data_eval_all$pp_mean <- out_loc
  data_eval_all$pp_std  <- out_sc
  data_eval_all[(data_eval_all$pp_mean<0 & !is.na(data_eval_all$pp_mean)), 'pp_mean'] <- 0
  # Transform back from square-root precipitation 
  data_eval_all$obs     <- data_eval_all$obs^2
  data_eval_all$fc_mean <- data_eval_all$fc_mean^2
  data_eval_all$pp_mean <- data_eval_all$pp_mean^2
}else{
  data_eval_all$pp_mean <- out_loc
  data_eval_all$pp_std  <- out_sc 
}
data_eval_all <- data_eval_all[!(data_eval_all$station_names%in%skip_stations), ]

# Model evaluation by scoringrule: 
print("-------------------")
print("CRPS scores:")
print("-------------------")
scores <- get_crps_scores(eval=data_eval_all, target=target)

print(paste0("--> Raw Forecast: "  , round(mean(scores[[1]], na.rm=TRUE), 3)))
print(paste0("--> ", exp_name, ": ", round(mean(scores[[2]], na.rm=TRUE), 3)))
print("")

print("CRPS score summary")
print("--> Raw Forecast:") 
print(summary(scores[[1]]))
print("")
print(paste0("--> ", exp_name, ":"))
print(summary(scores[[2]]))
print("")

df_out <- data.frame(date = data_eval_all$date,
                     station = data_eval_all$station_names,
                     mean = data_eval_all$pp_mean,
                     std  = data_eval_all$pp_std
                     )

# Output the NGR forecast of square-root precipitation
write.csv(df_out,
          file = file.path(out_data_dir, paste0("pred_", exp_name, "_", exp_suffix, ".csv")),
          row.names = FALSE)

print("Finished Successfully")
