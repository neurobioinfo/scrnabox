#!/usr/bin/env Rscript

##########################################
# v1.38
# step5 -- integration and linear dimensional reduction
##########################################

## load parameters
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
.libPaths(r_lib_path)

## load library
packages<-c('Seurat','ggplot2', 'dplyr', 'foreach', 'doParallel','Matrix')
lapply(packages, library, character.only = TRUE)

## load a list of existsing Seurat objects
sample_name<-list.files(path = paste(output_dir, "/step4/objs4",sep=""),pattern = "*.rds")
if(length(sample_name)<1) {
   print("You do not have any object from step 2 ")
}
for (i in 1:length(sample_name)) {
  if (!grepl(".rds",tolower(sample_name[i]), fixed = TRUE)){
     print(c(sample_name[i],"is not R rds"))
  }
}   

## load parameters
source(paste(output_dir,'/job_info/parameters/step5_par.txt',sep=""))

## detect available cores for parallel processing
numCores <- detectCores()
cl <- makeCluster(numCores-1)
registerDoParallel(cl) 

## create empty list to be populated by existing Seurat objects
seu_list<-list()

################################################################################################################################################
## if users want to integrate Seurat objects
################################################################################################################################################

if (tolower(par_integrate_seurat)=='yes') {
## load exisiting Seurat objects, find variable features, and scale
seu_list<-foreach (i_s=1:length(sample_name)) %do% {  
    seu<-readRDS(paste(output_dir,'/step4/objs4/',sample_name[i_s], sep=""))
    DefaultAssay(seu) <- par_DefaultAssay
    seu<- NormalizeData(seu,normalization.method = par_normalization.method, scale.factor =par_scale.factor) ##new code
    seu<- FindVariableFeatures(seu, selection.method = par_selection.method, nfeatures = par_nfeatures)
    seu<- ScaleData(seu, verbose = FALSE)
}  

## select inetgration features
features <- Seurat::SelectIntegrationFeatures(object.list = seu_list, nfeatures = par_nfeatures)

## find integration anchors
seu_anchors <- Seurat::FindIntegrationAnchors(object.list = seu_list, dims = 1:par_FindIntegrationAnchors_dim, anchor.features = features)

## integrate Seurat objects
seu_int <- Seurat::IntegrateData(anchorset = seu_anchors,dims = 1:par_FindIntegrationAnchors_dim)

## set default assay to integrated
Seurat::DefaultAssay(seu_int) <- "integrated"

## scale integrated Seurat object
seu_int <- ScaleData(seu_int, verbose = FALSE)

## run PCA and UMAP on integrated Seurat object
seu_int <- RunPCA(seu_int, npcs = par_RunPCA_npcs, verbose = FALSE)
seu_int <- RunUMAP(seu_int, dims = 1:par_RunUMAP_dims, n.neighbors =par_RunUMAP_n.neighbors)

## print PCA
DimPlot(seu_int, reduction = "pca", group.by="orig.ident", raster = FALSE )
ggsave(paste(output_dir,'/step5/figs5',"/integrated_DimPlot_pca.pdf", sep=""))

## print elbow plot
ElbowPlot(seu_int, ndims = par_RunPCA_npcs)  
ggsave(paste(output_dir,'/step5/figs5',"/integrated_elbow.pdf", sep=""))

## print UMAP
DimPlot(seu_int, reduction = "umap", group.by="orig.ident", raster = FALSE )
ggsave(paste(output_dir,'/step5/figs5',"/integrated_DimPlot_umap.pdf", sep=""))

## print Jack straw plot
if (tolower(par_compute_jackstraw)=='yes') {
seu_int <- JackStraw(seu_int, num.replicate = 100,dims = par_RunPCA_npcs)
seu_int <- ScoreJackStraw(seu_int, dims = 1:par_RunPCA_npcs)
JackStrawPlot(seu_int, dims = 1:par_RunPCA_npcs)
ggsave(paste(output_dir,'/step5/figs5',"/integrated_Jackstraw_plot.pdf", sep=""))
}

## save integrated Seurat object as RDS
saveRDS(seu_int, paste(output_dir,'/step5/objs5',"/seu_step5.rds", sep=""))

## save metadata info
write.csv(colnames(seu_int[[]]), paste(output_dir,'/step5/info5',"/integrated_meta_info_seu_step5.csv", sep=""))
}


################################################################################################################################################
## if users do not want to integrate Seurat objects (i.e. only have one Seurat object)
################################################################################################################################################

if (tolower(par_one_seurat)=='yes') {
## load exisiting Seurat objects, find variable features, and scale
seu_list<-foreach (i_s=1:length(sample_name)) %do% {  
    seu_int<-readRDS(paste(output_dir,'/step4/objs4/',sample_name[i_s], sep=""))
    DefaultAssay(seu_int) <- par_DefaultAssay
    seu_int<- NormalizeData(seu_int,normalization.method = par_normalization.method, scale.factor =par_scale.factor) ##new code
    seu_int<- FindVariableFeatures(seu_int, selection.method = par_selection.method, nfeatures = par_nfeatures)
    seu_int<- ScaleData(seu_int, verbose = FALSE)
}  

## run PCA and UMAP on  Seurat object
seu_int <- RunPCA(seu_int, npcs = par_RunPCA_npcs, verbose = FALSE)
seu_int <- RunUMAP(seu_int, dims = 1:par_RunUMAP_dims, n.neighbors =par_RunUMAP_n.neighbors)

## print PCA
DimPlot(seu_int, reduction = "pca", group.by="orig.ident", raster = FALSE )
ggsave(paste(output_dir,'/step5/figs5',"/DimPlot_pca.pdf", sep=""))

## print elbow plot
ElbowPlot(seu_int, ndims = par_RunPCA_npcs)   
ggsave(paste(output_dir,'/step5/figs5',"/elbow.pdf", sep=""))

## print UMAP
DimPlot(seu_int, reduction = "umap", group.by="orig.ident", raster = FALSE )
ggsave(paste(output_dir,'/step5/figs5',"/DimPlot_umap.pdf", sep=""))

## print Jack straw plot (optional)
if (tolower(par_compute_jackstraw)=='yes') {
seu_int <- JackStraw(seu_int, num.replicate = 100,dims = par_RunPCA_npcs) 
seu_int <- ScoreJackStraw(seu_int, dims = 1:par_RunPCA_npcs)
JackStrawPlot(seu_int, dims = 1:par_RunPCA_npcs)
ggsave(paste(output_dir,'/step5/figs5',"/Jackstraw_plot.pdf", sep=""))
}

## save Seurat object as RDS
saveRDS(seu_int, paste(output_dir,'/step5/objs5',"/seu_step5.rds", sep=""))

## save metadata info
write.csv(colnames(seu_int[[]]), paste(output_dir,'/step5/info5',"/meta_info_seu_step5.csv", sep=""))
}


################################################################################################################################################
## if users just want to merge their Seurat objects
################################################################################################################################################
if (tolower(par_merge_seurat)=='yes') {
## normalize exisiting Seurat objects
seu_list<-foreach (i_s=1:length(sample_name)) %do% {  
    seu<-readRDS(paste(output_dir,'/step4/objs4/',sample_name[i_s], sep=""))
    DefaultAssay(seu) <- par_DefaultAssay
    seu <- Seurat::NormalizeData(seu,normalization.method = par_normalization.method,scale.factor =par_scale.factor)
}  

## merge Seurat objects
length_list <- length(seu_list)

vars <- seu_list[2:length_list]

seu_int <- merge(seu_list[1], dput(as.character(vars)), add.cell.ids = dput(as.character(seu_list)), project = "merge",
                         merge.data = TRUE)


## set default assay to integrated
Seurat::DefaultAssay(seu_int) <- "RNA"

## find variable features
seu_int<- FindVariableFeatures(seu_int, selection.method = par_selection.method, nfeatures = par_nfeatures)

## scale integrated assay
seu_int <- ScaleData(seu_int, verbose = FALSE)

## run PCA and UMAP on integrated Seurat object
seu_int <- RunPCA(seu_int, npcs = par_RunPCA_npcs, verbose = FALSE)
seu_int <- RunUMAP(seu_int, dims = 1:par_RunUMAP_dims, n.neighbors =par_RunUMAP_n.neighbors)

## print PCA
DimPlot(seu_int, reduction = "pca", group.by="orig.ident", raster = FALSE )
ggsave(paste(output_dir,'/step5/figs5',"/merge_DimPlot_pca.pdf", sep=""))

## print elbow plot
ElbowPlot(seu_int)
ggsave(paste(output_dir,'/step5/figs5',"/merge_elbow.pdf", sep=""))

## print UMAP
DimPlot(seu_int, reduction = "umap", group.by="orig.ident", raster = FALSE)
ggsave(paste(output_dir,'/step5/figs5',"/merge_DimPlot_umap.pdf", sep=""))

## print Jack straw plot (optional)
if (tolower(par_compute_jackstraw)=='yes') {
seu <- JackStraw(seu, num.replicate = 100,dims = par_RunPCA_npcs) #added this figure
seu <- ScoreJackStraw(seu, dims = 1:par_RunPCA_npcs)
JackStrawPlot(seu, dims = 1:par_RunPCA_npcs)
ggsave(paste(output_dir,'/step5/figs5',"/merge_Jackstraw_plot.pdf", sep=""))
}

## save Seurat object as RDS
saveRDS(seu_int, paste(output_dir,'/step5/objs5',"/seu_step5.rds", sep=""))

## save metadata information
write.csv(colnames(seu_int[[]]), paste(output_dir,'/step5/info5',"/merge_meta_info_seu_step5.csv", sep=""))

}



## save RNA expression matrix
if (tolower(par_save_RNA)=='yes') {
    mat <- GetAssayData(object = seu_int, assay = "RNA", slot = "data")
    #  write.csv(mat, paste(output_dir,'/step2/info2/',sample_name[i],"_RNA.csv", sep=""))
    writeMM(mat,file= paste(output_dir,'/step5/info5/seu_int',"_RNA.txt", sep=""))
}

## save metadata dataframe
if (tolower(par_save_metadata)=='yes') {
    write.csv(seu_int[[]], file = paste(output_dir,'/step5/info5/seu_int_MetaData.txt', sep=""), quote = TRUE, sep = ",")
}

## save session info
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step5/info5/sessionInfo.txt', sep=""))
if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}
