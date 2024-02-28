#!/usr/bin/env Rscript

####################
# step7 -- module score
####################

stepp0="Module score"
cat("##########################################################################\n")
start_time0 <- Sys.time()
cat(stepp0,"has commenced.\n")
cat("##########################################################################\n")

stepp="Loading libraries and configuring parameters"
cat("#####################################\n")
cat(stepp, "started\n")
start_time <- Sys.time()

## set sample ID metadata column -- this is standard and does not require parameter modification
par_level_genotype <- "Sample_ID"

## load parameters
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
pipeline_home=args[3]

## load library
.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr', 'xlsx')
invisible(lapply(packages, library, character.only = TRUE))

## load existing Seurat objects
sample_name<-list.files(path = paste(output_dir, "/step6/objs6",sep=""),pattern = "*.rds")
seu_int<-readRDS(paste(output_dir,'/step6/objs6/',sample_name, sep=''))

## load parameters
source(paste(output_dir,'/job_info/parameters/step7_par.txt',sep=""))

cat(stepp,"has been achieved. Total time:",as.numeric (Sys.time() - start_time, units = "mins"),"minutes\n")
cat("#####################################\n")

stepp="Computing module score"
cat("#####################################\n")
cat(stepp, "started\n")
start_time <- Sys.time()


## set cell identity to the clustering resolution defined by the user
Idents(seu_int) <- par_level_cluster

## create directories for module score
## figures 
OUT_DIR_figs <- paste(output_dir,"/step7/figs7",sep='') 
OUT_dir_figs_module_score <- paste(OUT_DIR_figs,"/module_score/",sep='') 
dir.create(OUT_dir_figs_module_score)
## info
OUT_DIR_info <- paste(output_dir,"/step7/info7",sep='') 
OUT_dir_info_module_score <- paste(OUT_DIR_info,"/module_score/",sep='') 
dir.create(OUT_dir_info_module_score)

## set output directory
PWD=OUT_dir_info_module_score
setwd(PWD)

## compute module score
source(paste(pipeline_home,'/tools/module_score.R',sep=''))

cat(stepp,"has been achieved. Total time:",as.numeric (Sys.time() - start_time, units = "mins"),"minutes\n")
cat("#####################################\n")

## save session info
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step7/info7/sessionInfo_module_score.txt', sep=""))
if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}

cat("##########################################################################\n")
cat(stepp0,"successfully completed. Total time:",as.numeric (Sys.time() - start_time0, units = "mins"),"minutes\n")
cat("##########################################################################\n")

