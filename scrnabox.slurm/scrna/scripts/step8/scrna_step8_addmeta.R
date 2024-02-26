#!/usr/bin/env Rscript

####################
# step8 -- add Metadata
####################
stepp0="Step VIII addmeta"
cat("##########################################################################\n")
start_time0 <- Sys.time()
cat(stepp0,"has commenced.\n")
cat("##########################################################################\n")

stepp="Loading libraies and parameters"
cat("#####################################\n")
cat(stepp, "started\n")
start_time <- Sys.time()
## load parameters
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]

## load library
.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr','stringi','tidyverse','Matrix', 'ggrepel','DESeq2','EnhancedVolcano')
invisible(lapply(packages, library, character.only = TRUE))

## load parameters
source(paste(output_dir,'/job_info/parameters/step8_par.txt',sep=""))

if (exists("par_seurat_object")) {                                                  
    seu_int<-readRDS(par_seurat_object)
} else {
    sample_name<-list.files(path = paste(output_dir, "/step7/objs7",sep=""),pattern = "*.rds")
    seu_int<-readRDS(paste(output_dir,'/step7/objs7/',sample_name, sep=''))
}  

cat(stepp,"has been achieved. Total time:",as.numeric (Sys.time() - start_time, units = "mins"),"minutes\n")
cat("#####################################\n")

stepp="Adding metadata"
cat("#####################################\n")
cat(stepp, "started\n")
start_time <- Sys.time()

## load meatadata
meta_data <- read.delim(par_metadata, header = T, sep = ",") 
new_meta <- colnames(meta_data)

## create existing metdata dataframe
metadata_df <- data.frame(seu_int@meta.data)

## merge existing metdata and new metdata
merge_metadata_df <- merge(metadata_df,meta_data, by = par_merge_meta )
df<- merge_metadata_df[,new_meta]
nrow(metadata_df) == nrow(df)

## set colnames into a list
colnames <- colnames(df)
nrow(df)

## add metadata
for (i in 1:ncol(df)){
seu_int <- AddMetaData(seu_int, metadata = df[,i], col.name = colnames[i])
}

## print RDS object
saveRDS(seu_int, paste(output_dir,'/step8/objs8',"/seu_step8.rds", sep=""))


cat(stepp,"has been achieved. Total time:",as.numeric (Sys.time() - start_time, units = "mins"),"minutes\n")
cat("#####################################\n")

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
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step8/info8/sessionInfo_add_metadata.txt', sep=""))

if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}


cat("##########################################################################\n")
cat(stepp0,"successfully completed. Total time:",as.numeric (Sys.time() - start_time0, units = "mins"),"minutes\n")
cat("##########################################################################\n")
