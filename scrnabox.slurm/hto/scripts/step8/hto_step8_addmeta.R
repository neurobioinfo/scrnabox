#!/usr/bin/env Rscript

####################
# step8 -- add Metadata
####################

## load parameters
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]

## load library
.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr','stringi','tidyverse','Matrix', 'ggrepel','DESeq2','EnhancedVolcano')
lapply(packages, library, character.only = TRUE)

## load parameters
source(paste(output_dir,'/job_info/parameters/step8_par.txt',sep=""))

## load step7 seurat object
if (exists("par_seurat_object")) {                                                  
    seu_int<-readRDS(par_seurat_object)
} else {
    sample_name<-list.files(path = paste(output_dir, "/step7/objs7",sep=""),pattern = "*.rds")
    seu_int<-readRDS(paste(output_dir,'/step7/objs7/',sample_name, sep=''))
}  

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

