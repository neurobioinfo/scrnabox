#!/usr/bin/env Rscript

####################
# step3 
## Create seurat object 

args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
NFRNAL=as.numeric(args[3])
NFRNAU=as.numeric(args[4])
PMT=as.numeric(args[5])


print(output_dir)
print(NFRNAL)
print(NFRNAU)
print(PMT)

.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr')
# output_dir="/home/samamiri/scratch/Darkgenome_run/pipeline/scrnabox.svn_run/des"
list<-read.csv(paste(output_dir, "/job_output/sample_dir.list",sep=""),header=FALSE)
sample_name<-read.csv(paste(output_dir, "/job_output/sample.list",sep=""),header=FALSE)
set.seed(1234)
library(foreach)
library(doParallel)

numCores <- detectCores()
cl <- makeCluster(numCores-1)
registerDoParallel(cl) 

foreach (i=1:dim(sample_name)[1]) %do% {    
    set.seed(1234)
    seu<-readRDS(paste(output_dir,'/step2/objs',"/seurat_object.",sample_name$V1[i],".rds", sep=""))
    DA.hash<-seu
    hashtags <- rownames(DA.hash[["HTO"]])
    seu <-  DA.hash
    seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")
    seu <- subset(seu, subset = nFeature_RNA > NFRNAL & nFeature_RNA < NFRNAU & percent.mt < PMT)
    # seu <- subset(seu, subset = nFeature_RNA > 300 & nFeature_RNA < 6500 & percent.mt < 25)
    saveRDS(seu, paste(output_dir,'/step3/objs',"/seu",i,".rds", sep=""))
}

