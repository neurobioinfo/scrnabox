#!/usr/bin/env Rscript

####################
# step4 -- msd (pre-step for doublet detection and demultiplexing)
####################

## load parameters
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
.libPaths(r_lib_path)

## set seed for replicability
set.seed(1234)

## load library
packages<-c('Seurat','ggplot2', 'dplyr','foreach', 'doParallel')
lapply(packages, library, character.only = TRUE)

## create list of existing Seurat objects from step 3
sample_name<-list.files(path = paste(output_dir, "/step3/objs3",sep=""),pattern = "*.rds")
print("list of sample")
print(sample_name)

if(length(sample_name)<1) {
   print("You do not have any object from step 3 ")
}

for (i in 1:length(sample_name)) {
  if (!grepl(".rds",tolower(sample_name[i]), fixed = TRUE)){
     print(c(sample_name[i],"is not R rds"))
  }
}  

## load parameters text file
source(paste(output_dir,'/job_info/parameters/step4_par.txt',sep=""))

## detect number of available cores
numCores <- detectCores()
cl <- makeCluster(numCores-1)
registerDoParallel(cl) 

## retrieve and print antibody labels
foreach (i_s=1:length(sample_name)) %do% {   
    seu<-readRDS(paste(output_dir,'/step3/objs3/',sample_name[i_s], sep=""))
    DefaultAssay(seu) <- "HTO"
    seu <- MULTIseqDemux(seu, assay = "HTO", quantile = par_quantile, autoThresh = par_autoThresh, maxiter = par_maxiter) 
    write.csv(unique(seu$MULTI_ID), paste(output_dir,'/step4/info4/',sample_name[i_s],"_old_antibody_label_MULTIseqDemuxHTOcounts.csv",sep=""))
}
