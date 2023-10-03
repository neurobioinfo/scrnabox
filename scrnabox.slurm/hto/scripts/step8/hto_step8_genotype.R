#!/usr/bin/env Rscript

####################
# v1.38
# step8 -- add sample-sample contrast
####################

## load parameters
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]

## load library
.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr','stringi','tidyverse','Matrix')
lapply(packages, library, character.only = TRUE)

## load step7 seurat object
sample_name<-list.files(path = paste(output_dir, "/step7/objs7",sep=""),pattern = "*.rds")
seu_int<-readRDS(paste(output_dir,'/step7/objs7/',sample_name, sep=''))

## load parameters
source(paste(output_dir,'/job_info/parameters/step8_par.txt',sep=""))

################### import the right Seurat object ###################
## load name of existing Seurat objects
sample_name<-list.files(path = paste(output_dir, "/step7/objs7",sep=""),pattern = "*.rds")

if(file.exists(paste(output_dir,'/step8/objs8/','seu_step8.rds', sep = ""))){
    seu_int<-readRDS(paste(output_dir,'/step8/objs8/','seu_step8.rds', sep=''))
}else{
    seu_int<-readRDS(paste(output_dir,'/step7/objs7/',sample_name, sep=''))
}
################### ############################## ###################

## add info directory
OUT_DIR_info <- paste(output_dir,"/step8/info8",sep='') 
OUT_dir_info_sample <- paste(OUT_DIR_info,"/sample_sample_contrasts/",sep='') 
dir.create(OUT_dir_info_sample)


################### sample-sample dge analysis ###################
dd<-read.csv(paste(output_dir,'/job_info/parameters/step8_contrast_genotype.txt',sep='/'), sep="")

for(i in 1:nrow(dd)){
Idents(seu_int) <- dd[i,2]    
DGE <- FindMarkers(seu_int, ident.1 = dd[i,3], ident.2 = dd[i,4],  logfc.threshold = 0)
#write dge
write.csv(DGE, file = paste(OUT_dir_info_sample, dd[i,1],'_DEG.csv', sep=""), quote = FALSE, sep = ",")
}
################### ########################## ###################

## print RDS object
saveRDS(seu_int, paste(output_dir,'/step8/objs8',"/seu_step8.rds", sep=""))

## save rna expression matrix
if (tolower(par_save_RNA)=='yes') {
    mat <- GetAssayData(object = seu_int, assay = "RNA", slot = "data")
    #  write.csv(mat, paste(output_dir,'/step2/info2/',sample_name[i],"_RNA.csv", sep=""))
    writeMM(mat,file= paste(output_dir,'/step8/info8/seu',"_RNA.txt", sep=""))
}

## save metadata dataframe
if (tolower(par_save_metadata)=='yes') {
    write.csv(seu_int[[]], file = paste(output_dir,'/step8/info8/seu_MetaData.txt', sep=""), quote = TRUE, sep = ",")
}

## write session info
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step8/info8/sessionInfo_sample_sample_cont.txt', sep=""))
file.remove("Rplots.pdf")

