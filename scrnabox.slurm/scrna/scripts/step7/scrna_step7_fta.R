#!/usr/bin/env Rscript

####################
# step7 -- reference-based annotation
####################

## load parameters
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
.libPaths(r_lib_path)

## load libraries
packages<-c('Seurat','ggplot2', 'dplyr')
lapply(packages, library, character.only = TRUE)

## load existing Seurat objects
sample_name<-list.files(path = paste(output_dir, "/step6/objs6",sep=""),pattern = "*.rds")
seu_int<-readRDS(paste(output_dir,'/step6/objs6/',sample_name, sep=''))

## load parameters
source(paste(output_dir,'/job_info/parameters/step7_par.txt',sep=""))

## set output directory
PWD=paste(output_dir,'/step7/objs7/', sep='')
setwd(PWD)

## create reference-based annotation-specific directory
OUT_DIR_figs <- paste(output_dir,"/step7/figs7",sep='') 
OUT_dir_figs_reference <- paste(OUT_DIR_figs,"/reference_based_annotation/",sep='') 
dir.create(OUT_dir_figs_reference)

## set user defined clustering resolution
seu_int <- SetIdent(seu_int, value = par_level_cluster)

## load reference Seurat object
reference0 <-readRDS(par_reference)

## load parallelization parameters
options(future.globals.maxSize = par_futureglobalsmaxSize)

## find transfer anchors between reference and query Seurat objects
transfer.anchors <- FindTransferAnchors(reference = reference0, query = seu_int, dims = 1:par_FindTransferAnchors_dim, reference.reduction = "pca")

## add reference-based annotations to the qeury object
eval(parse(text = paste('predictions <- TransferData(anchorset = transfer.anchors, refdata = reference0$',par_level_celltype ,',dims = 1:',par_FindTransferAnchors_dim,')', sep='')))
seu_int <- AddMetaData(object = seu_int, metadata = predictions)

## save query Seurat object with reference annotation predicitions
saveRDS(seu_int,paste(output_dir,'/step7/objs7','/seu_step7.rds', sep=''))

## save metadata information
write.csv(colnames(seu_int[[]]), file= paste(output_dir,'/step7/info7/meta_info_seu_step7',".txt", sep=""))

## Print a umap projection showing the predicted cell types on the query object 
reference0 <- RunUMAP(reference0, dims = 1:30, reduction = "pca", return.model = TRUE)
seu_int <- MapQuery(anchorset = transfer.anchors, reference = reference0, query = seu_int,
    refdata = list(celltype = par_level_celltype), reference.reduction = "pca", reduction.model = "umap")
p1 <- DimPlot(reference0, reduction = "umap", group.by = par_level_celltype, label = TRUE, label.size = 3,repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(seu_int, reduction = "umap", group.by = "predicted.id", label = TRUE, label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2
ggsave(file = paste(OUT_dir_figs_reference,'UMAP_transferred_labels.pdf', sep=''))

## save RNA expression matrix
if (tolower(par_save_RNA)=='yes') {
    mat <- GetAssayData(object = seu_int, assay = "RNA", slot = "data")
    writeMM(mat,file= paste(output_dir,'/step7/info7/seu',"_RNA.txt", sep=""))
}

## save metadata 
if (tolower(par_save_metadata)=='yes') {
    write.csv(seu_int[[]], file = paste(output_dir,'/step7/info7/seu_MetaData.txt', sep=""), quote = TRUE, sep = ",")
}

## save session info
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step7/info7/sessionInfo_fta.txt', sep=""))
if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}
