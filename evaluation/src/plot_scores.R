#!/usr/bin/env Rscript                                                                                                                                                                       
#--------------------------------------------------------------------                                                                                                                        
# Visualisation input score files
# Usage: ./plot_scores.R ScoreFiles[regex]                                                                                                                         
# Author: Patrik Benacek                                                                                                                                                                     
#----------------------------------------------                                                                                                                                              
rm(list=ls())  

pkgTest <- function(x)
  {
    if (!require(x,character.only = TRUE))
    {
      install.packages(x,dep=TRUE)
        if(!require(x,character.only = TRUE)) stop("Package not found")
    }
  }

suppressMessages(pkgTest("tidyverse"))
suppressMessages(pkgTest("ggmap"))

# Input arguments                                                                                                                                                                            
args = commandArgs(trailingOnly=TRUE)                                                                                                                                                        

files = Sys.glob(args)

# Paths                                                                                                                                                                                    
outpath <- "results/plots"
obspath <- "../observations/data"

is_reference_exp = FALSE
allscores = data.frame()
for (i in seq_along(files)){
  file = files[i]
  print(paste0("--> Read file: ", file))
  # Check if file exist
  if (!file.exists(file)){
    stop(paste0("Input file", file, "does not exist."))
  }
  # Read exp metadata
  filename = gsub('.Rdata', '', basename(file))
  meta = strsplit(filename, "_")[[1]]
  # Check the name of the ScoreFile
  if (length(meta)!=6){
    stop("Input filename has not valid name. Use e.g.: eval_scores_EMOS-global_t2m_ff24h_2015.Rdata")
  }                                                                                                                                                                                            
  # Check the consistency of input ScoreFiles
  if (i==1){
    expname     = meta[3]  
    target      = meta[4] 
    fc_time     = meta[5]
    start_train = meta[6]
  }else{
    expname     = meta[3]
    if(meta[4] != target) stop(paste("ScoreFiles Target Inconsistency:", target, 'and', meta[4]))
    if(meta[5] != fc_time) stop(paste("ScoreFiles Leadtime Inconsistency:", fc_time, 'and', meta[5]))
    if(meta[6] != start_train) warning(paste("ScoreFiles TrainPeriod Inconsistency:", start_train, 'and', meta[6]))
  }
  # Bind expdata
  allscores = rbind(allscores, readRDS(file)) 
}

# Experiment summary
expnames <- unique(allscores$exp)
print(paste("The following experiments are compared:"))
print(expnames)

if ('Raw-Fcst' %in% expnames){
  is_reference_exp = TRUE
  print("Reference exist.")
}

# Read stations' metedata
obs_meta = read.csv(file.path(obspath, 'metadata_stations.csv'), 
                    col.names = c('station_id', 'station', 'lon', 'lat', 'alt'))
# trimws(lstrip+rstrip) + substitute
obs_meta$station = trimws(gsub(' / ', '-', obs_meta$station))

# Merge station with scorefiles
allscores = merge(allscores, obs_meta, by = 'station')

#--------------------------
# Visualisation
#--------------------------
# Save plot function
savefig <- function(plot, name, type='png', width=2.2*1400, height=2.2*960, res=2*300){
  fig_suffix <- paste0(paste(target, fc_time, start_train, sep="_"), '.', type)
  if (type=='png'){
    png(file.path(outpath, paste('plot', name, fig_suffix, sep='_')), width=width, height=height, res=res)
    print(plot)
    dev.off()
  }else if(type=='pdf'){
    pdf(file.path(outpath, paste('plot', name, fig_suffix, sep='_')), width=width, height=height)
    print(plot)
    dev.off()
  }
}

# PIT scores
pl <- allscores %>% 
  ggplot(aes(pit)) +
  geom_histogram(aes(y=..density..), color='black', fill='white', bins=20) + 
  ylim(c(0,1.5)) + 
  geom_hline(yintercept = 1, lty=2) +
  xlab('PIT') + #ggtitle("PIT histograms") +
  facet_wrap(~exp) +
  theme_bw()
savefig(pl, name="PIT_overall")

# CRPS scores
pl <- allscores %>% 
  group_by(exp) %>% 
  summarise(
    mCRPS = mean(crps)
  ) %>% 
  ggplot(aes(x=factor(exp), y=mCRPS, fill=factor(exp))) +
  geom_col(color='black', fill='grey') + 
  #geom_col(color='black') + # use to colorize columns
  geom_text(aes(label = round(mCRPS, 2), y=0), vjust = -0.5, size=3) +
  xlab('Model') + ylab('Mean CRPS') +
  theme_bw()
savefig(pl, name="CRPS_overall")

# Absolute Error
pl <- allscores %>% 
  group_by(exp) %>% 
  summarise(
    mAE = mean(ae)
  ) %>% 
  ggplot(aes(x=factor(exp), y=mAE, fill=factor(exp))) +
  geom_col(color='black', fill='grey') + 
  #geom_col(color='black') + # use to colorize columns
  geom_text(aes(label = round(mAE, 2), y=0), vjust = -0.5, size=3) +
  xlab('Model') + ylab('AE of mean forecast') +
  theme_bw()
savefig(pl, name="AE_overall")
  
# CRPS score values in domain
register_google(key="AIzaSyCM_lvOP7nXiLHFuV96DL4mrXBitJWRHEc")

map <- get_googlemap(                                                                                                                                                                        
  center = c(15.55, 49.8), zoom = 7, maptype = "terrain", scale = 1,
  style = 'feature:road|element:all|visibility:off&style=feature:all|element:labels|visibility:off&sensor=false'
  #color = "bw"
)                                                                                                                                                                                            

# Prepare domain crps scores
crps_scores <- allscores %>% 
  group_by(exp, station) %>% 
  summarise(
    mCRPS = mean(crps)
  ) %>% 
  left_join(obs_meta, by = 'station')

pl <- ggmap(map) + 
  geom_point(data = crps_scores, aes(x = lon, y = lat), colour = "black", size=2.5) + 
  geom_point(data = crps_scores, aes(x = lon, y = lat, colour = mCRPS), size=2) + 
  scale_colour_gradient(low = "blue", high = "red") + 
  labs(x = "Longitude", y = "Latitude", col="Mean CRPS") + 
  theme(legend.position="bottom") +
  facet_wrap(~exp)
savefig(pl, name="CRPS_domain", height = 2.2*1400, width=2.2*1400)

# Relative crps scores wrt RAW forecast: 1-(exp/RAW) where 0=RAWquality, 1=100%improvement)
if (is_reference_exp){
  spread_crps <- spread(crps_scores, exp, mCRPS)
  for (exp in expnames){
    if (exp=='Raw-Fcst') next()
    spread_crps[,exp] = 100*(1-(spread_crps[,exp]/spread_crps[,"Raw-Fcst"]))
  }
  crps_scores_rel <- spread_crps %>% select(-alt, -station_id,-`Raw-Fcst`) %>% 
    gather(exp, rCRPS, expnames[expnames!="Raw-Fcst"]) 

  pl <- ggmap(map) + 
    geom_point(data = crps_scores_rel, aes(x = lon, y = lat), colour = "black", size=2.5) + 
    geom_point(data = crps_scores_rel, aes(x = lon, y = lat, colour = rCRPS), size=2) + 
    scale_colour_gradient2(low = "red", mid="white", high = "blue", midpoint=0) + 
    labs(x = "Longitude", y = "Latitude", col="CRPS improvement \n wrt RawFc [%]") + 
    theme(legend.position="bottom") +
    facet_wrap(~exp)
  savefig(pl, name="rCRPS_domain", height = 2.2*1400, width=2.2*1400)
}
