#!/usr/bin/env Rscript

####################
# step4
## Create seurat object 
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
.libPaths(r_lib_path)
# source(args[3])
set.seed(1234)
packages<-c('Seurat','ggplot2', 'dplyr')
lapply(packages, library, character.only = TRUE)

# output_dir="/home/samamiri/scratch/Darkgenome_run/pipeline/scrnabox.svn_run/des"
list<-read.csv(paste(output_dir, "/job_output/sample_dir.list",sep=""),header=FALSE)
sample_name<-read.csv(paste(output_dir, "/job_output/sample.list",sep=""),header=FALSE)

library(foreach)
library(doParallel)
numCores <- detectCores()
cl <- makeCluster(numCores-1)
registerDoParallel(cl) 


foreach (i_s=1:dim(sample_name)[1]) %do% {   
    set.seed(1234)
    seu1<-readRDS(paste(output_dir,'/step3/objs',"/seu",i_s,".rds", sep=""))
    DefaultAssay(seu1) <- "HTO"
    seu1 <- NormalizeData(seu1, assay = "HTO", normalization.method = "CLR")
    seu1 <- FindVariableFeatures(seu1, selection.method = "mean.var.plot")
    seu1 <- ScaleData(seu1, features = VariableFeatures(seu1))
    seu1 <- MULTIseqDemux(seu1, assay = "HTO", quantile = 0.9, autoThresh = TRUE, maxiter = 5) 
    write.csv(table(seu1$MULTI_ID), paste(output_dir,'/step4/objs/seu',i_s,"MULTIseqDemuxHTOcounts.csv",sep=""))
}