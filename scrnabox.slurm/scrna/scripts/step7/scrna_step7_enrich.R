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
# seu_int<-readRDS(paste(output_dir,'/step6/objs/',sample_name, sep=''))

# seu_int<-readRDS(paste(output_dir,'/step6/objs','/seu_step6.rds', sep=''))
source(paste(output_dir,'/job_info/parameters/step7_par.txt',sep=""))
PWD=paste(output_dir,'/step7/info7/', sep='')
setwd(PWD)
# seu_int <- SetIdent(seu_int, value = level_cluster)
# reference0 <-readRDS(reference)
options(future.globals.maxSize = futureglobalsmaxSize)
# transfer.anchors <- FindTransferAnchors(reference = reference0, query = seu_int, dims = 1:30)
# predictions <- TransferData(anchorset = transfer.anchors, refdata = reference0[[level_celltype]],dims = 1:30)
# eval(parse(text = paste('predictions <- TransferData(anchorset = transfer.anchors, refdata = reference0$',level_celltype ,',dims = 1:30)', sep='')))
# seu_int <- AddMetaData(object = seu_int, metadata = predictions)
# saveRDS(seu_int,paste(output_dir,'/step7/objs','/seu_step7.rds', sep=''))


##################
# level_cluster='integrated_snn_res.0.7'
# PWD='~/Desktop/annot/'
PSUE=paste(output_dir,'/step6/objs6/',sample_name, sep='')
# top_sel=5
# db <- c('Descartes_Cell_Types_and_Tissue_2021','CellMarker_Augmented_2021','Azimuth_Cell_Types_2021')
scrnaboxR::annotation(level_cluster,PWD,PSUE,top_sel,db)

writeLines(capture.output(sessionInfo()), paste(output_dir,'/step7/info7/sessionInfo_enrich.txt', sep=""))

