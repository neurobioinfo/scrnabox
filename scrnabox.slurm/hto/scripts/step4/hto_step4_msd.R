#!/usr/bin/env Rscript

####################
# step4 -- msd (pre-step for doublet detection and demultiplexing)
####################
stepp0="Step IV msd"
cat("##########################################################################\n")
start_time0 <- Sys.time()
cat(stepp0,"has commenced.\n")
cat("##########################################################################\n")

stepp="Loading libraries and configuring parameters"
cat("#####################################\n")
cat(stepp, "started\n")
start_time <- Sys.time()

## load parameters
# rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
.libPaths(r_lib_path)

## set seed for replicability
set.seed(1234)

## load library
packages<-c('Seurat','ggplot2', 'dplyr','foreach', 'doParallel')
invisible(lapply(packages, library, character.only = TRUE))

## load Seurat objects
if (exists("par_seurat_object")) {                                                  
    stepp="Step I"
    cat("#####################################\n")
    cat("The step",stepp, "started\n")
    start_time <- Sys.time()
    sample_name<-list.files(path = par_seurat_object)
    sample_nameb<-gsub(".rds","",sample_name)
    if(length(sample_name)<1) {
        print("You do not have any existing Seurat object")
    }
    for (i in 1:length(sample_name)) {
        if (!grepl(".rds",tolower(sample_name[i]), fixed = TRUE)){
            print(c(sample_name[i],"is not R rds"))
        }
        }
    cat("The step",stepp,"finished. Total time to achive:",as.numeric (Sys.time() - start_time, units = "mins"),"minutes\n")
    cat("#####################################\n")
} else {
    stepp="Step I"
    cat("#####################################\n")
    cat("The step",stepp, "started\n")
    start_time <- Sys.time()
    sample_name<-list.files(path = paste(output_dir, "/step3/objs3",sep=""))
    sample_nameb<-gsub(".rds","",sample_name)
    if(length(sample_name)<1) {
        print("You do not have any object from step 3 ")
    }
    for (i in 1:length(sample_name)) {
        if (!grepl(".rds",tolower(sample_name[i]), fixed = TRUE)){
            print(c(sample_name[i],"is not R rds"))
        }
        }
    cat("The step",stepp,"finished. Total time to achive:",as.numeric (Sys.time() - start_time, units = "mins"),"minutes\n")
    cat("#####################################\n")
}

## load parameters text file
source(paste(output_dir,'/job_info/parameters/step4_par.txt',sep=""))

cat(stepp,"has been achieved. Total time:",as.numeric (Sys.time() - start_time, units = "mins"),"minutes\n")
cat("#####################################\n")

## detect number of available cores
numCores <- detectCores()
cl <- makeCluster(numCores-1)
registerDoParallel(cl) 

stepp="retrieve and print antibody labels"
cat("#####################################\n")
cat(stepp, "started\n")
start_time <- Sys.time()
## retrieve and print antibody labels
foreach (i_s=1:length(sample_name)) %do% {  

    seu<-readRDS(paste(output_dir,'/step3/objs3/',sample_name[i_s], sep=""))
    DefaultAssay(seu) <- "HTO"
    seu <- MULTIseqDemux(seu, assay = "HTO", quantile = par_quantile, autoThresh = par_autoThresh, maxiter = par_maxiter) 
    write.csv(unique(seu$MULTI_ID), paste(output_dir,'/step4/info4/',sample_name[i_s],"_old_antibody_label_MULTIseqDemuxHTOcounts.csv",sep=""))
    cat("The step",stepp,"finished. Total time to achive:",as.numeric (Sys.time() - start_time, units = "mins"),"minutes\n")
    cat("#####################################\n")
}
cat(stepp,"has been achieved. Total time:",as.numeric (Sys.time() - start_time, units = "mins"),"minutes\n")
cat("#####################################\n")

cat("##########################################################################\n")
cat(stepp0,"successfully completed. Total time:",as.numeric (Sys.time() - start_time0, units = "mins"),"minutes\n")
cat("##########################################################################\n")

