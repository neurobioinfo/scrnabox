#!/usr/bin/env Rscript

####################
# step6

args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
# scrna_method=args[3]
.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr','cowplot','clustree')
lapply(packages, library, character.only = TRUE)
# library('cowplot')
# library('clustree')

sample_name<-list.files(path = paste(output_dir, "/step5/objs5",sep=""),pattern = "*.rds")
if(length(sample_name)<1) {
   print("You do not have any object from step 2 ")
}

for (i in 1:length(sample_name)) {
  if (!grepl(".rds",tolower(sample_name[i]), fixed = TRUE)){
     print(c(sample_name[i],"is not R rds"))
  }
}   

seu_int<-readRDS(paste(output_dir,'/step5/objs5/',sample_name, sep=''))


source(paste(output_dir,'/job_info/parameters/step6_par.txt',sep=""))


DefaultAssay(seu_int) <- "RNA"
seu_int <- NormalizeData(seu_int)
seu_int <- ScaleData(seu_int)
# Do the PCA
seu_int <- FindVariableFeatures(seu_int, selection.method = par_selection.method, nfeatures = par_nfeatures)
seu_int <- RunPCA(seu_int)  # should use variable ft by default

png(paste(output_dir,'/step6/figs6',"/elbow_pca.png", sep=""))
fig1<-DimPlot(seu_int, reduction = "pca")
fig2<-ElbowPlot(seu_int)
print(fig1 + fig2)
dev.off()


#######
Seurat::DefaultAssay(seu_int) <- "integrated"
seu_int <- ScaleData(seu_int, verbose = FALSE)
seu_int <- RunPCA(seu_int, npcs = par_npcs_pca, verbose = FALSE)
seu_int <- RunUMAP(seu_int, dims = 1:par_dims_umap, n.neighbors =par_n.neighbors)
seu_int <- FindNeighbors(seu_int,  dims = 1:par_dims_fn, k.param = par_k.param, prune.SNN = par_prune.SNN)
seu_int <- Seurat::FindClusters(seu_int, resolution = par_resolution)


clustree(seu_int@meta.data, prefix = "integrated_snn_res.")
ggsave(paste(output_dir,'/step6/figs6',"/clustree.png", sep=""))
for (i in par_resolution){
Seurat::DimPlot(seu_int, group.by = paste("integrated_snn_res.",i,sep=''))
ggsave(paste(output_dir,'/step6/figs6',"/integrated_snn_res.",i,".png",sep=""))
}

########


saveRDS(seu_int, paste(output_dir,'/step6/objs6',"/seu_step6.rds", sep=""))
write.csv(colnames(seu_int[[]]), paste(output_dir,'/step6/info6',"/meta_info.csv", sep=""))


# if (scrna_method=='HTO') {
    # write.table(unique(seu_int$MULTI_ID_Lables),file=paste(output_dir,'/step6/info6',"/class_label_phenotypes.txt", sep=""),col.names=FALSE)
# }



if (tolower(Save_RNA)=='yes') {
    mat <- GetAssayData(object = seu_int, assay = "RNA", slot = "data")
    writeMM(mat,file= paste(output_dir,'/step6/info6/seu',"_RNA.txt", sep=""))
}

if (tolower(Save_metadata)=='yes') {
    write.csv(seu_int[[]], file = paste(output_dir,'/step6/info6/seu_MetaData.txt', sep=""), quote = TRUE, sep = ",")
}

writeLines(capture.output(sessionInfo()), paste(output_dir,'/step6/info6/sessionInfo.txt', sep=""))

