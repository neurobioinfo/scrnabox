#!/usr/bin/env Rscript

####################
# step2 
## Create seurat object 

args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
.libPaths(r_lib_path)


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

seu_list<-list()

seu_list<-foreach (i_s=1:dim(sample_name)[1]) %do% {  
    seu<-readRDS(paste(output_dir,'/step4/objs',"/seu",i_s,"dem.rds", sep=""))
    DefaultAssay(seu) <- "RNA"
    seu <- Seurat::NormalizeData(seu)
    seu <- Seurat::FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
   }  

seu_anchors <- Seurat::FindIntegrationAnchors(object.list = seu_list)
seu_int <- Seurat::IntegrateData(anchorset = seu_anchors)
Seurat::DefaultAssay(seu_int) <- "integrated"
saveRDS(seu_int, paste(output_dir,'/step5/objs',"/seu_int.rds", sep=""))


