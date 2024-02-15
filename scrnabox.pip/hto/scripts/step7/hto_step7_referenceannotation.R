#!/usr/bin/env Rscript

###############################################################################
# step7 -- refrence-based annotation
###############################################################################

## load parameters
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
pipeline_home=args[3]

## load library
.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr', 'xlsx', 'Matrix')
lapply(packages, library, character.only = TRUE)

## load parameters
source(paste(output_dir,'/job_info/parameters/step7_par.txt',sep=""))

################### import the right Seurat object ###################
## load name of existing Seurat objects
sample_name<-list.files(path = paste(output_dir, "/step6/objs6",sep=""),pattern = "*.rds")

if (exists("par_seurat_object")) {                                                   
    seu_int<-readRDS(par_seurat_object)
}else{
    if(file.exists(paste(output_dir,'/step7/objs7/','seu_step7.rds', sep = ""))){
        seu_int<-readRDS(paste(output_dir,'/step7/objs7/','seu_step7.rds', sep=''))
    }else{
        seu_int<-readRDS(paste(output_dir,'/step6/objs6/',sample_name, sep=''))
    }
}
################### ############################## ###################


## set output directory
PWD=paste(output_dir,'/step7/objs7/', sep='')
setwd(PWD)

## create reference-based directory
#figures
OUT_DIR_figs <- paste(output_dir,"/step7/figs7",sep='') 
OUT_dir_figs_reference <- paste(OUT_DIR_figs,"/reference_based_annotation/",sep='') 
dir.create(OUT_dir_figs_reference)
#info
OUT_DIR_info <- paste(output_dir,"/step7/info7",sep='') 
OUT_dir_info_reference <- paste(OUT_DIR_info,"/reference_based_annotation/",sep='') 
dir.create(OUT_dir_info_reference)

## set user defined clustering resolution
seu_int <- SetIdent(seu_int, value = par_level_cluster)

## load reference Seurat object
reference0 <-readRDS(par_reference)
DefaultAssay(reference0) <- "RNA" ## new code

## load parallelization parameters
options(future.globals.maxSize = par_futureglobalsmaxSize)

# perform standard preprocessing on reference object
reference0<- NormalizeData(reference0)
reference0 <- FindVariableFeatures(reference0)
reference0<- ScaleData(reference0)
reference0 <- RunPCA(object = reference0, assay = "RNA", npcs = par_FindTransferAnchors_dim)

## find transfer anchors between reference and query Seurat objects
transfer.anchors <- FindTransferAnchors(reference = reference0, query = seu_int, dims = 1:par_FindTransferAnchors_dim, reference.reduction = "pca")

## add reference-based annotations to the qeury object
eval(parse(text = paste('predictions <- TransferData(anchorset = transfer.anchors, refdata = reference0$',par_level_celltype ,',dims = 1:',par_FindTransferAnchors_dim,')', sep='')))
seu_int <- AddMetaData(object = seu_int, metadata = predictions)

# Add metadata column for reference object
seu_int$temp_temp_2 <- seu_int@meta.data$predicted.id
name_meta <- names(seu_int@meta.data) 
length <- length(name_meta)
name_meta[length] <- paste(par_reference_name, "_predictions", sep = "")
names(seu_int@meta.data) <- name_meta

## save query Seurat object with reference annotation predicitions
saveRDS(seu_int,paste(output_dir,'/step7/objs7','/seu_step7.rds', sep=''))

## save metadata information
write.csv(colnames(seu_int[[]]), file= paste(output_dir,'/step7/info7/meta_info_seu_step7',".txt", sep=""))

## Print a umap projection showing the predicted cell types on the query object 
reference0 <- RunUMAP(reference0, dims = 1:par_FindTransferAnchors_dim, reduction = "pca", return.model = TRUE)
seu_int <- MapQuery(anchorset = transfer.anchors, reference = reference0, query = seu_int,
    refdata = list(celltype = par_level_celltype), reference.reduction = "pca", reduction.model = "umap")
p1 <- DimPlot(reference0, reduction = "umap", group.by = par_level_celltype, label = TRUE, label.size = 3,repel = TRUE, raster=FALSE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(seu_int, reduction = "umap", group.by = "predicted.id", label = TRUE, label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2
ggsave(file = paste(OUT_dir_figs_reference,par_reference_name,'_UMAP_transferred_labels.pdf', sep=''), dpi = 300, height = 7, width = 14, unit = 'in')

## print summary table
df <- data.frame(seu_int@meta.data)
cluster_list <- list(par_level_cluster)
keep_list <- append(cluster_list, "predicted.id")
df <- df[, (colnames(df) %in% keep_list)]
df_summary <- table(df[,1], df[,2])
df_summary <- data.frame(df_summary)
colnames(df_summary) <- c("cluster", "cell_type", "number_of_cells")
    
if(file.exists(paste(output_dir,'/step7/info7/reference_based_annotation/',par_reference_name, '_prediction_summary.xlsx', sep=""))){
file.remove(paste(output_dir,'/step7/info7/reference_based_annotation/',par_reference_name, '_prediction_summary.xlsx', sep=""))
}
    
write.xlsx(df_summary, file=paste(output_dir,'/step7/info7/reference_based_annotation/',par_reference_name, '_prediction_summary.xlsx', sep=""),row.names=FALSE)

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
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step7/info7/sessionInfo_reference.txt', sep=""))
if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}
