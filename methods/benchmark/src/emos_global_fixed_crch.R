#!/usr/bin/env Rscript                                                                                                                                                                       

##---------------------------------
## Run selected benchmark models
##---------------------------------
## * EMOS global, fixed window
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
exp_name <- "EMOS-global"
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
print(paste0("All cases:", nrow(data_train_all)))
print(paste0("Zero precipitation:", sum(data_train_all$obs==0)))
t1 <- Sys.time()
# Heteroscedastic censored regression with CRPS minimization. Left censored by 0. 
emos_crps <- get_crch_model(train=data_train_all, target=target)
t2 <- Sys.time()
print(paste0("Computational time: ", round(t2-t1, 1), "s"))

if (target=='prec24'){
  # Make prediction
  data_eval_all$pp_mean <- predict(emos_crps, type='location', newdata = data_eval_all[,c('fc_mean', 'fc_std')])
  data_eval_all$pp_std  <- predict(emos_crps, type='scale'   , newdata = data_eval_all[,c('fc_mean', 'fc_std')])
  data_eval_all[(data_eval_all$pp_mean<0 & !is.na(data_eval_all$pp_mean)), 'pp_mean'] <- 0
  # Transform back from square-root precipitation 
  data_eval_all$obs     <- data_eval_all$obs^2
  data_eval_all$fc_mean <- data_eval_all$fc_mean^2
  data_eval_all$pp_mean <- data_eval_all$pp_mean^2
}else{
  data_eval_all$pp_mean <- predict(emos_crps, type='location', newdata = data_eval_all[,c('fc_mean', 'fc_std')])
  data_eval_all$pp_std  <- predict(emos_crps, type='scale'   , newdata = data_eval_all[,c('fc_mean', 'fc_std')])
}

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

# Scatter plot
#png(file.path(out_data_dir, paste0("plot_scatt_", exp_name, "_", exp_suffix, '.png')))
#m <- lm(obs~fc_mean, data = data)
#plot(obs~fc_mean, data = data[data$obs<max(data$fc_mean, na.rm=TRUE),], pch=20, col='lightblue')
#abline(0,1, col='black')
#abline(m, col='red')
#abline(emos_crps, col='blue')
#legend('topright', c('1-to-1', 'OLS', 'NGR'), lty=1, col = c('black', 'red', 'blue'), cex=.7)
#dev.off()

print("Finished Successfully")
