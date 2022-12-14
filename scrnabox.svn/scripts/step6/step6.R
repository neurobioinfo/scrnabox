#!/usr/bin/env Rscript

####################
# step6
## Create seurat object 

args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr')
lapply(packages, library, character.only = TRUE)
library('cowplot')
library('clustree')

seu_int<-readRDS(paste(output_dir,'/step5/objs','/seu_int.rds', sep=''))

source(paste(output_dir,'/job_output/step6_par.txt',sep=""))


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
# seu_int <- FindNeighbors(seu_int, dims = 1:30)
seu_int <- FindNeighbors(seu_int,  dims = 1:par_dims_fn, k.param = par_k.param, prune.SNN = par_prune.SNN)
seu_int <- Seurat::FindClusters(seu_int, resolution = par_resolution)

# png(paste(output_dir,'/step6/figs',"/clustree.png", sep=""))
# clustree(seu_int@meta.data, prefix = "integrated_snn_res.")
# dev.off()

# png(paste(output_dir,'/step6/figs',"/dim_10.png", sep=""))
# Seurat::DimPlot(seu_int, group.by = "integrated_snn_res.0.1")
# dev.off()

# png(paste(output_dir,'/step6/figs',"/dim_12.png", sep=""))
# Seurat::DimPlot(seu_int, group.by = "integrated_snn_res.0.2")
# dev.off()

# png(paste(output_dir,'/step6/figs',"/dim_15.png", sep=""))
# Seurat::DimPlot(seu_int, group.by = "integrated_snn_res.0.3")
# dev.off()

# png(paste(output_dir,'/step6/figs',"/dim_17.png", sep=""))
# Seurat::DimPlot(seu_int, group.by = "integrated_snn_res.0.7")
# dev.off()

# png(paste(output_dir,'/step6/figs',"/dim.png", sep=""))
# Seurat::DimPlot(seu_int)
# dev.off()

clustree(seu_int@meta.data, prefix = "integrated_snn_res.")
ggsave(paste(output_dir,'/step6/figs6',"/clustree.png", sep=""))
for (i in par_resolution){
Seurat::DimPlot(seu_int, group.by = paste("integrated_snn_res.",i,sep=''))
ggsave((paste(output_dir,'/step6/figs6','integrated_snn_res.',i, sep=''))
}

########

write.table(unique(seu_int$Lables),file=paste(output_dir,'/step6/objs',"/class_label_phenotypes.txt", sep=""),col.names=FALSE)

saveRDS(seu_int, paste(output_dir,'/step6/objs',"/seu_int_clu.rds", sep=""))
write.csv(colnames(seu_int[[]]), paste(output_dir,'/step6/objs',"/mata_info.csv", sep=""))


