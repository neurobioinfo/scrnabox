#!/usr/bin/env Rscript

####################
# v1.38
# step6 -- clustering
####################

## load parameters
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
pipeline_home=args[3]

## load library
.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr','cowplot','clustree','Matrix')
lapply(packages, library, character.only = TRUE)

## load parameters
source(paste(output_dir,'/job_info/parameters/step6_par.txt',sep=""))

## determine the default assay depending on whether users performed integration or not
if (tolower(par_skip_integration)=='yes') {
    sample_name<-list.files(path = paste(output_dir, "/step5/objs5",sep=""),pattern = "*.rds")
    if(length(sample_name)<1) {
    print("You do not have any object from step 5")
    }
    for (i in 1:length(sample_name)) {
    if (!grepl(".rds",tolower(sample_name[i]), fixed = TRUE)){
        print(c(sample_name[i],"is not R rds"))
    }
    }   
    seu_int<-readRDS(paste(output_dir,'/step5/objs5/',sample_name, sep=''))
    par_whatAssay<-'RNA'
} else {
    sample_name<-list.files(path = paste(output_dir, "/step5/objs5",sep=""),pattern = "*.rds")
    if(length(sample_name)<1) {
    print("You do not have any object from step 5")
    }
    for (i in 1:length(sample_name)) {
    if (!grepl(".rds",tolower(sample_name[i]), fixed = TRUE)){
        print(c(sample_name[i],"is not R rds"))
    }
    }   
    seu_int<-readRDS(paste(output_dir,'/step5/objs5/',sample_name, sep=''))
    par_whatAssay<-'integrated'
}

## set default assay
par_whatAssay
Seurat::DefaultAssay(seu_int) <-  par_whatAssay 

## cluster at each user-defined clustering resolution
seu_int <- FindNeighbors(seu_int,  dims = 1:par_FindNeighbors_dims, k.param = par_FindNeighbors_k.param, prune.SNN = par_FindNeighbors_prune.SNN)
seu_int <- Seurat::FindClusters(seu_int, resolution = par_FindClusters_resolution)
seu_int <- RunUMAP(seu_int, dims = 1:par_RunUMAP_dims)

## print clustree plot 
clustree(seu_int@meta.data, prefix = paste0(par_whatAssay,"_snn_res."))
ggsave(paste(output_dir,"/step6/figs6/clustree_int.pdf", sep=""), dpi = 300, height = 9, width = 5, unit = 'in')

## print UMAPs at different clustering resolutions
for (i in par_FindClusters_resolution){
    Seurat::DimPlot(seu_int, group.by = paste(par_whatAssay,"_snn_res.",i,sep=''))
    ggsave(paste(output_dir,'/step6/figs6/',par_whatAssay,"_snn_res.",i,".pdf",sep=""))
}

## save Seurat object
saveRDS(seu_int, paste(output_dir,'/step6/objs6',"/seu_step6.rds", sep=""))

## save meta information
write.csv(colnames(seu_int[[]]), paste(output_dir,'/step6/info6',"/meta_info.csv", sep=""))

## print sample labels
write.table(unique(seu_int$MULTI_ID_Lables),file=paste(output_dir,'/step6/info6',"/class_label_phenotypes.txt", sep=""),col.names=FALSE)

## perform ARI computation to evaluate the stability of clusters at different clustering resolutions
if (tolower(par_compute_ARI)=='yes') {
source(paste(pipeline_home,'/general_codes/Rand_index.R',sep=''))
}

## save RNA expression matrix
if (tolower(par_save_RNA)=='yes') {
    mat <- GetAssayData(object = seu_int, assay = "RNA", slot = "data")
    writeMM(mat,file= paste(output_dir,'/step6/info6/seu',"_RNA.txt", sep=""))
}

## save metadata dataframe
if (tolower(par_save_metadata)=='yes') {
    write.csv(seu_int[[]], file = paste(output_dir,'/step6/info6/seu_MetaData.txt', sep=""), quote = TRUE, sep = ",")
}

## write session information
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step6/info6/sessionInfo.txt', sep=""))
if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}