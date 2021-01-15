#!/usr/bin/env Rscript

DATADIR=file.path("/home/patrik/Work/czechglobe/TIGGE/data_preproc/generate_dataset/data")

data = read.csv(file.path(DATADIR, "data_wmeta_prec24_ff240h_2015_2019.csv"))
data$date = as.Date(data$date)

libus = subset(data, station=11520)

pdf(file.path(DATADIR, "plot_control_deacc.pdf"))
plot(libus$date, libus$prec24_mean, type="p", pch=20, cex=.2)
plot(libus$date, libus$slhf_mean, type="p", pch=20, cex=.2)
dev.off()
