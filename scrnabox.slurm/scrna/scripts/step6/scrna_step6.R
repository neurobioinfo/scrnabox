#!/usr/bin/env Rscript

####################
# step 6: clustering
####################
stepp0="Step IV"
cat("##########################################################################\n")
start_time0 <- Sys.time()
cat(stepp0,"has commenced.\n")
cat("##########################################################################\n")

stepp="Loading libraries and configuring parameters"
cat("#####################################\n")
cat(stepp, "started\n")
start_time <- Sys.time()
## load parameters
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
pipeline_home=args[3]

## load libraries
.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr','cowplot','clustree', 'Matrix')
invisible(lapply(packages, library, character.only = TRUE))

## load parameters
source(paste(output_dir,'/job_info/parameters/step6_par.txt',sep=""))

## load Seurat object
if (exists("par_seurat_object")) {                                                  
    sample_name<-list.files(path = par_seurat_object)
    if(length(sample_name)<1) {
    print("You do not have any existing Seurat object")
    }
    for (i in 1:length(sample_name)) {
    if (!grepl(".rds",tolower(sample_name[i]), fixed = TRUE)){
        print(c(sample_name[i],"is not R rds"))
    }
    }  
    seu_int<-readRDS(par_seurat_object)
} else {
    sample_name<-list.files(path = paste(output_dir, "/step5/objs5",sep=""),pattern = "*.rds")
    if(length(sample_name)<1) {
    print("You do not have any object from step 3 ")
    }
    for (i in 1:length(sample_name)) {
    if (!grepl(".rds",tolower(sample_name[i]), fixed = TRUE)){
        print(c(sample_name[i],"is not R rds"))
    }
    }
    seu_int<-readRDS(paste(output_dir,'/step5/objs5/',sample_name, sep=''))
}  

## determine the default assay
if (tolower(par_skip_integration)=='yes') { 
    par_whatAssay<-'RNA'
} else {
    par_whatAssay<-'integrated'
}
cat(stepp,"has been achieved. Total time:",as.numeric (Sys.time() - start_time, units = "mins"),"minutes\n")
cat("#####################################\n")
stepp="FindNeighbors"
cat("#####################################\n")
cat(stepp, "started\n")
start_time <- Sys.time()
## set default assay depending on whether integration was performed or not
    Seurat::DefaultAssay(seu_int) <- par_whatAssay
    ## cluster at each user-defined clustering resolution
    seu_int <- FindNeighbors(seu_int,  dims = 1:par_FindNeighbors_dims, k.param = par_FindNeighbors_k.param, prune.SNN = par_FindNeighbors_prune.SNN)
cat(stepp,"has been achieved. Total time:",as.numeric (Sys.time() - start_time, units = "mins"),"minutes\n")
cat("#####################################\n")
stepp="FindClusters"
cat("#####################################\n")
cat(stepp, "started\n")
start_time <- Sys.time()
    seu_int <- Seurat::FindClusters(seu_int, resolution = par_FindClusters_resolution)
cat(stepp,"has been achieved. Total time:",as.numeric (Sys.time() - start_time, units = "mins"),"minutes\n")
cat("#####################################\n")
stepp="RunUMAP"
cat("#####################################\n")
cat(stepp, "started\n")
start_time <- Sys.time()
    seu_int <- RunUMAP(seu_int, dims = 1:par_RunUMAP_dims)
cat(stepp,"has been achieved. Total time:",as.numeric (Sys.time() - start_time, units = "mins"),"minutes\n")
cat("#####################################\n")
stepp="clustree"
cat("#####################################\n")
cat(stepp, "started\n")
start_time <- Sys.time()
## print clustree plot 
    clustree(seu_int@meta.data, prefix = paste0(par_whatAssay,"_snn_res."))
    ggsave(paste(output_dir,"/step6/figs6/clustree_int.pdf", sep=""), dpi = 300, height = 9, width = 8, unit = 'in')
    ## print UMAPs at different resolutions
    for (i in par_FindClusters_resolution){
        Seurat::DimPlot(seu_int, group.by = paste(par_whatAssay,"_snn_res.",i,sep=''))
        ggsave(paste(output_dir,'/step6/figs6/',par_whatAssay,"_snn_res.",i,".pdf",sep=""))
    }
cat(stepp,"has been achieved. Total time:",as.numeric (Sys.time() - start_time, units = "mins"),"minutes\n")
cat("#####################################\n")
## save Seurat object as RDS 
saveRDS(seu_int, paste(output_dir,'/step6/objs6',"/seu_step6.rds", sep=""))

## print metadata information
write.csv(colnames(seu_int[[]]), paste(output_dir,'/step6/info6',"/meta_info.csv", sep=""))

## perform ARI computation to evaluate the stability of clusters at different clustering resolutions
if (tolower(par_compute_ARI)=='yes') {
source(paste(pipeline_home,'/tools/Rand_index.R',sep=''))
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

## write session info
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step6/info6/sessionInfo.txt', sep=""))
if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}

cat("##########################################################################\n")
cat(stepp0,"successfully completed. Total time:",as.numeric (Sys.time() - start_time0, units = "mins"),"minutes\n")
cat("##########################################################################\n")
