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

seu_int<-readRDS(paste(output_dir,'/step6/objs','/seu_int_clu.rds', sep=''))
source(paste(output_dir,'/job_output/step7_par.txt',sep=""))
PWD=paste(output_dir,'/step7/objs/', sep='')
setwd(PWD)
seu_int <- SetIdent(seu_int, value = level_cluster)
reference0 <-readRDS(reference)
options(future.globals.maxSize = futureglobalsmaxSize)
transfer.anchors <- FindTransferAnchors(reference = reference0, query = seu_int, dims = 1:30)
# predictions <- TransferData(anchorset = transfer.anchors, refdata = reference0[[level_celltype]],dims = 1:30)
eval(parse(text = paste('predictions <- TransferData(anchorset = transfer.anchors, refdata = reference0$',level_celltype ,',dims = 1:30)', sep='')))
seu_int <- AddMetaData(object = seu_int, metadata = predictions)
saveRDS(seu_int,paste(output_dir,'/step7/objs','/seu_annot.rds', sep=''))
