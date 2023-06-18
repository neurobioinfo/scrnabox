#!/usr/bin/env Rscript

####################
# step4
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
.libPaths(r_lib_path)

set.seed(1234)
packages<-c('Seurat','ggplot2', 'dplyr')
lapply(packages, library, character.only = TRUE)

sample_name<-list.files(path = paste(output_dir, "/step3/objs3",sep=""),pattern = "*.rds")

print("list of sample")
print(sample_name)

if(length(sample_name)<1) {
   print("You do not have any object from step 2 ")
}

for (i in 1:length(sample_name)) {
  if (!grepl(".rds",tolower(sample_name[i]), fixed = TRUE)){
     print(c(sample_name[i],"is not R rds"))
  }
}   


library(foreach)
library(doParallel)
numCores <- detectCores()
cl <- makeCluster(numCores-1)
registerDoParallel(cl) 

foreach (i_s=1:length(sample_name)) %do% {   
    set.seed(1234)
    seu1<-readRDS(paste(output_dir,'/step3/objs3/',sample_name[i_s], sep=""))
    DefaultAssay(seu1) <- "HTO"
    seu1 <- NormalizeData(seu1, assay = "HTO", normalization.method = "CLR")
    seu1 <- FindVariableFeatures(seu1, selection.method = "mean.var.plot")
    seu1 <- ScaleData(seu1, features = VariableFeatures(seu1))
    seu1 <- MULTIseqDemux(seu1, assay = "HTO", quantile = 0.9, autoThresh = TRUE, maxiter = 5) 
    write.csv(table(seu1$MULTI_ID), paste(output_dir,'/step4/info4/seu',i_s,"MULTIseqDemuxHTOcounts.csv",sep=""))
}
