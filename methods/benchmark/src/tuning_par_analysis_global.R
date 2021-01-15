#!/usr/bin/env Rscript                                                                                                                                                                       

## -----------------------------------------------------------------
## Analysis code for qrf results / tuning parameter influence
## -----------------------------------------------------------------

rm(list=ls())

in_data_dir <- "/home/patrik/Work/czechglobe/TIGGE/postproc_benchmark/results/qrf_tuning"

# memTest version [+ pre-existing files]
ntree.try <- c(125,250,500,1000)
nodesize.try <- c(5,10,15,20,50)
mtry.try <- c(10,15,20,25,30,40,50)
replace.try <- c(1)

pars <- expand.grid(as.factor(ntree.try), 
                    as.factor(nodesize.try), 
                    as.factor(mtry.try), 
                    as.factor(replace.try))
names(pars) <- c("ntree", "nodesize", "mtry", "replace")

savenames <- paste0("qrf_global_ntree", pars[,1], "_nodesize", pars[,2], "_mtry", pars[,3], "_repl", pars[,4], ".Rdata")

res <- rep(NA, length(savenames))
for(filename in savenames){
  file = file.path(in_data_dir, filename)
  if(file.exists(file)){
    load(file)
    res[which(savenames == filename)] <- mean(qrf_crps)
  }
}

pars$crps <- res
head(pars)

pars[with(pars, order(crps)), ]
