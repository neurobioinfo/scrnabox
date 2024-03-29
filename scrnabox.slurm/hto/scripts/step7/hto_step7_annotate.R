#!/usr/bin/env Rscript

###############################################################################
# step7 -- Annotate clusters
###############################################################################
## load parameters
stepp0="Step VII annotate"
cat("##########################################################################\n")
start_time0 <- Sys.time()
cat(stepp0,"has commenced.\n")
cat("##########################################################################\n")

stepp="Loading libraries and configuring parameters"
cat("#####################################\n")
cat(stepp, "started\n")
start_time <- Sys.time()

args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
pipeline_home=args[3]

## load library
.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr', 'xlsx', 'Matrix')
invisible(lapply(packages, library, character.only = TRUE))

## load parameters
source(paste(output_dir,'/job_info/parameters/step7_par.txt',sep=""))

################### import the right Seurat object ###################
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

cat(stepp,"has been achieved. Total time:",as.numeric (Sys.time() - start_time, units = "mins"),"minutes\n")
cat("#####################################\n")
################### ############################## ###################

stepp="Cell identity"
cat("#####################################\n")
cat(stepp, "started\n")
start_time <- Sys.time()

## set cell identity to the clustering resolution defined by the user
Idents(seu_int) <- par_annotate_resolution

## create directories for annotation cluster 
## figures 
OUT_DIR_figs <- paste(output_dir,"/step7/figs7",sep='') 
OUT_dir_figs_annotate <- paste(OUT_DIR_figs,"/annotate/",sep='') 
dir.create(OUT_dir_figs_annotate)

## add cluster annotation
cluster.ids<-par_annotate_labels

## set cluster resolution and rename cluster identities
seu_int <- SetIdent(seu_int, value = par_annotate_resolution)
names(cluster.ids) <- levels(seu_int)    
seu_int <- RenameIdents(seu_int, cluster.ids) 
seu_int$temp_temp_1 <- Idents(seu_int)

name_meta <- names(seu_int@meta.data) 
length <- length(name_meta)
name_meta[length] <- par_name_metadata
names(seu_int@meta.data) <- name_meta

##save RDS object
saveRDS(seu_int,paste(output_dir,'/step7/objs7','/seu_step7.rds', sep=''))

cat(stepp,"has been achieved. Total time:",as.numeric (Sys.time() - start_time, units = "mins"),"minutes\n")
cat("#####################################\n")

stepp="DimPlot"
cat("#####################################\n")
cat(stepp, "started\n")
start_time <- Sys.time()
## print UMAP with cluster annotation
DimPlot(seu_int, reduction = "umap", label = TRUE, pt.size = 0.5, raster = FALSE) + NoLegend()
ggsave(file = paste(OUT_dir_figs_annotate,par_name_metadata,'_cluster_annotation.pdf', sep=''))

## print UMAP splitted
DimPlot(seu_int, reduction = "umap", split.by = "Sample_ID", label = TRUE, pt.size = 0.5, ncol = 3, raster = FALSE) + NoLegend()
ggsave(file = paste(OUT_dir_figs_annotate,par_name_metadata,'_split_cluster_annotation.pdf', sep=''), dpi = 300, height = 15, width = 15, unit = 'in')

cat(stepp,"has been achieved. Total time:",as.numeric (Sys.time() - start_time, units = "mins"),"minutes\n")
cat("#####################################\n")

## save metadata information
write.csv(colnames(seu_int[[]]), file= paste(output_dir,'/step7/info7/meta_info_seu_step7',".txt", sep=""))

## save rna expression matrix
if (tolower(par_save_RNA)=='yes') {
    mat <- GetAssayData(object = seu_int, assay = "RNA", slot = "data")
    writeMM(mat,file= paste(output_dir,'/step7/info7/seu',"_RNA.txt", sep=""))
}

## save metadata dataframe
if (tolower(par_save_metadata)=='yes') {
    write.csv(seu_int[[]], file = paste(output_dir,'/step7/info7/seu_MetaData.txt', sep=""), quote = TRUE, sep = ",")
}

## write session info
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step7/info7/sessionInfo_annotate.txt', sep=""))

if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}

cat("##########################################################################\n")
cat(stepp0,"successfully completed. Total time:",as.numeric (Sys.time() - start_time0, units = "mins"),"minutes\n")
cat("##########################################################################\n")

