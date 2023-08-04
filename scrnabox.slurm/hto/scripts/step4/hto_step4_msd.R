#!/usr/bin/env Rscript

####################
# step4
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
.libPaths(r_lib_path)

set.seed(1234)
packages<-c('Seurat','ggplot2', 'dplyr','foreach', 'doParallel')
lapply(packages, library, character.only = TRUE)

sample_name<-list.files(path = paste(output_dir, "/step3/objs3",sep=""),pattern = "*.rds")

print("list of sample")
print(sample_name)

if(length(sample_name)<1) {
   print("You do not have any object from step 2 ")
}

for (i in 1:length(sample_name)) {
  if (!grepl(".rds",tolower(sample_name[i]), fixed = TRUE)){
     print(c(sample_name[i],"is not R rds"))
  }
}   


numCores <- detectCores()
cl <- makeCluster(numCores-1)
registerDoParallel(cl) 

foreach (i_s=1:length(sample_name)) %do% {   
    seu<-readRDS(paste(output_dir,'/step3/objs3/',sample_name[i_s], sep=""))
    DefaultAssay(seu) <- "HTO"
    if (tolower(normlazation_and_scalaing)=='yes'){
        seu <- NormalizeData(seu, assay = "HTO", normalization.method = par_normalization.method,scale.factor =par_scale.factor)
        seu <- FindVariableFeatures(seu, selection.method = par_selection.method, nfeatures = par_nfeatures)
        topsel <- head(Seurat::VariableFeatures(seu), par_top)
        write.csv(topsel, file = paste(output_dir,'/step4/info4/most_variable_genes_',sample_name[i],'.txt', sep=""), quote = TRUE, sep = ",")
        vf_plot <- Seurat::VariableFeaturePlot(seu)
        Seurat::LabelPoints(plot = vf_plot,points = topsel, repel = TRUE)
        ggsave(paste(output_dir,'/step4/figs4/VariableFeaturePlot',sample_name[i_s],".png",sep=""))
        seu <- ScaleData(seu)
    }
    ## dimensionality reduction
    if (tolower(dimensionality_reduction)=='yes'){
        seu <- RunPCA(seu, npcs = par_npcs_pca, verbose = FALSE)
        DimPlot(seu, reduction = "pca")
        ggsave(paste(output_dir,'/step4/info4/',"dimplot_pca",sample_name[i_s],".png",sep=""))
        ElbowPlot(seu)
        ggsave(paste(output_dir,'/step4/info4/',"elbowplot",sample_name[i_s],".png",sep=""))
        Seurat::DimHeatmap(seu, dims = 1:par_dims, cells = par_cells, balanced = TRUE)
        ggsave(paste(output_dir,'/step4/info4/',"dimheatplot.",sample_name[i_s],".png",sep=""))
        FeaturePlot(seu, reduction = par_reduction, features = par_features)
        ggsave(paste(output_dir,'/step4/info4/',"featureplot",sample_name[i_s],".png",sep=""))
        seu <- RunUMAP(seu, dims = 1:par_dims_umap, n.neighbors =par_n.neighbors)
        Seurat::DimPlot(seu, reduction = "umap")
        ggsave(paste(output_dir,'/step4/figs4',"dimplot_umap",sample_name[i_s],".png",sep=""))
    }
    seu <- MULTIseqDemux(seu, assay = "HTO", quantile = par_quantile, autoThresh = par_autoThresh, maxiter = par_maxiter) 
    write.csv(table(seu$MULTI_ID), paste(output_dir,'/step4/info4/',sample_name[i_s],"MULTIseqDemuxHTOcounts.csv",sep=""))
}


