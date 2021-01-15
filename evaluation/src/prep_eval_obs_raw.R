#!/usr/bin/env Rscript                                                                                                                                                                       

exp_name <- 'Raw-Fcst'
target   <- 't2m'
leadtime <- 240

in_data_dir = "/home/patrik/Work/czechglobe/TIGGE/data_preproc/generate_dataset/data"
out_data_dir = "/home/patrik/Work/czechglobe/TIGGE/evaluation/results"

# Read station measurements                                                                                                                                                                  
file_obs = paste0("data_wmeta_", target, "_ff", leadtime, "h_2015_2019.csv")                                                                                                                    
data <- read.csv(file.path(in_data_dir, file_obs))                                                                                                                                           
data$date = as.Date(data$date)                                                                                                                                                               
  
data <- data[, -which(!(names(data) %in% c("obs", "date", "station", paste0(target,"_mean"), paste0(target,"_var"))))]                                                                       
names(data)[2:5] <- c('station_id', 'obs', 'mean', 'var')

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

saveRDS(data[,c('date', 'station_id', 'obs')], file = file.path(out_data_dir, paste0("eval_obs_", target, '_ff', leadtime, "h_2019.RData")))

write.csv(data[,c('date', 'station_id', 'mean', 'std')], file = file.path(out_data_dir, 'prediction', paste0('pred_', exp_name, "_", target, "_ff", leadtime ,"h_2019_2019.csv")), row.names = FALSE)
