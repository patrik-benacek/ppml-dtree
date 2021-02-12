#!/usr/bin/env Rscript                                                                                                                                                                       

pkgTest <- function(x)
  {
    if (!require(x,character.only = TRUE))
    {
      install.packages(x,dep=TRUE)
        if(!require(x,character.only = TRUE)) stop("Package not found")
    }
  }

suppressMessages(pkgTest("readr"))

# Input arguments
args = commandArgs(trailingOnly=TRUE)

exp_name = 'Raw-Fcst'
target   = args[1] 
leadtime = args[2]

in_data_dir = "../data_preproc/data/processed"

# Read station measurements                                                                                                                                                                  
file_obs = paste0("data_", target, "_ff", leadtime, "h.zip")                                                                                                                    
data <- read_csv(file.path(in_data_dir, file_obs))                                                                                                                                           
data$date = as.Date(data$date)                                                                                                                                                               
data$station_names = gsub(' / ', '-', data$station_names)
  
data <- data[, -which(!(names(data) %in% c("obs", "date", "station_names", paste0(target,"_mean"), paste0(target,"_var"))))]                                                                       
names(data)[2:5] <- c('station', 'obs', 'mean', 'var')

# Evaluation period
start_eval <- as.Date("2019-01-01 00:00", tz = "UTC")
end_eval <- as.Date("2019-12-31 00:00", tz = "UTC") 
data <- subset(data, date >= start_eval & date <= end_eval)

# Data cleaning: square-root scaled data (negative->NAN), remove NAN values, omit 'perfect' ensemble prediction                                                                              
data <- data[data$var>0,]                                                                                                                                                                 
data <- data[!is.na(data$obs),]                                                                                                                                                              
data <- data[!is.na(data$mean),]                                                                                                                                                          
data$std <- sqrt(data$var) 
data$var <- NULL

if (target=="prec24"){
    print("Data are transformed back from square-root transformation (applied before running ML models.)")
    data$obs  = data$obs**2
    data$mean = data$mean**2
}

saveRDS(data[,c('date', 'station', 'obs')], file = file.path('results', paste0("eval_obs_", target, '_ff', leadtime, "h_2019.RData")))

write.csv(data[,c('date', 'station', 'mean', 'std')], file = file.path('results', 'prediction', paste0('pred_', exp_name, "_", target, "_ff", leadtime ,"h_2019.csv")), row.names = FALSE)
