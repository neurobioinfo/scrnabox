#!/usr/bin/env Rscript

####################
# step7
## Annotation

args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr')
lapply(packages, library, character.only = TRUE)

sample_name<-list.files(path = paste(output_dir, "/step6/objs6",sep=""),pattern = "*.rds")
seu_int<-readRDS(paste(output_dir,'/step6/objs6/',sample_name, sep=''))

source(paste(output_dir,'/job_info/parameters/step7_par.txt',sep=""))
PWD=paste(output_dir,'/step7/objs7/', sep='')
setwd(PWD)
seu_int <- SetIdent(seu_int, value = level_cluster)
reference0 <-readRDS(reference)
options(future.globals.maxSize = futureglobalsmaxSize)
transfer.anchors <- FindTransferAnchors(reference = reference0, query = seu_int, dims = 1:30)
# predictions <- TransferData(anchorset = transfer.anchors, refdata = reference0[[level_celltype]],dims = 1:30)
eval(parse(text = paste('predictions <- TransferData(anchorset = transfer.anchors, refdata = reference0$',level_celltype ,',dims = 1:30)', sep='')))
seu_int <- AddMetaData(object = seu_int, metadata = predictions)
saveRDS(seu_int,paste(output_dir,'/step7/objs7','/seu_step7.rds', sep=''))
write.csv(colnames(seu_int[[]]), file= paste(output_dir,'/step7/info7/meta_info_seu_step7',".txt", sep=""))

if (tolower(Save_RNA)=='yes') {
    mat <- GetAssayData(object = seu_int, assay = "RNA", slot = "data")
    #  write.csv(mat, paste(output_dir,'/step2/info2/',sample_name[i],"_RNA.csv", sep=""))
    writeMM(mat,file= paste(output_dir,'/step7/info7/seu',"_RNA.txt", sep=""))
}

if (tolower(Save_metadata)=='yes') {
    write.csv(seu_int[[]], file = paste(output_dir,'/step7/info7/seu_MetaData.txt', sep=""), quote = TRUE, sep = ",")
}

writeLines(capture.output(sessionInfo()), paste(output_dir,'/step7/info7/sessionInfo_fta.txt', sep=""))

