#!/usr/bin/env Rscript

####################
# step4
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
.libPaths(r_lib_path)

set.seed(1234)
packages<-c('Seurat','ggplot2', 'dplyr','foreach', 'doParallel','Matrix','tidyverse')
lapply(packages, library, character.only = TRUE)


sample_name<-list.files(path = paste(output_dir, "/step3/objs3",sep=""),pattern = "*.rds")
sample_nameb<-gsub(".rds","",sample_name)

if(length(sample_name)<1) {
   print("You do not have any object from step 2 ")
}

for (i in 1:length(sample_name)) {
  if (!grepl(".rds",tolower(sample_name[i]), fixed = TRUE)){
     print(c(sample_name[i],"is not R rds"))
  }
}   


source(paste(output_dir,'/job_info/parameters/step4_par.txt',sep=""))
# library('tidyverse')
# dd<-read.csv(paste(output_dir,'/job_info/parameters/step4_antibody_label.txt',sep=''), header=FALSE, sep="")
# new.names<-str_split(dd[2,],",")[[1]]
# old.names<-str_split(dd[1,],",")[[1]]

old.names<-old_antibody_label
new.names<-new_antibody_label


numCores <- detectCores()
cl <- makeCluster(numCores-1)
registerDoParallel(cl) 


foreach (i_s=1:length(sample_name)) %do% {  
    set.seed(1234)
    seu<-readRDS(paste(output_dir,'/step3/objs3/',sample_name[i_s], sep=""))
    DefaultAssay(seu) <- "HTO"
    if (tolower(normlazation_and_scalaing)=='yes'){
        seu <- NormalizeData(seu, assay = "HTO", normalization.method = par_normalization.method,scale.factor =par_scale.factor)
        seu <- FindVariableFeatures(seu, selection.method = par_selection.method, nfeatures = par_nfeatures)
        topsel <- head(Seurat::VariableFeatures(seu), par_top)
        write.csv(topsel, file = paste(output_dir,'/step4/info4/most_variable_genes_',sample_nameb[i_s],'.txt', sep=""), quote = TRUE, sep = ",")
        vf_plot <- Seurat::VariableFeaturePlot(seu)
        Seurat::LabelPoints(plot = vf_plot,points = topsel, repel = TRUE)
        ggsave(paste(output_dir,'/step4/figs4/VariableFeaturePlot',sample_nameb[i_s],".png",sep=""))
        seu <- ScaleData(seu)
    }
    ## dimensionality reduction
    if (tolower(dimensionality_reduction)=='yes'){
        seu <- RunPCA(seu, npcs = par_npcs_pca, verbose = FALSE)
        DimPlot(seu, reduction = "pca")
        ggsave(paste(output_dir,'/step4/figs4/',"dimplot_pca",sample_nameb[i_s],".png",sep=""))
        ElbowPlot(seu)
        ggsave(paste(output_dir,'/step4/figs4/',"elbowplot",sample_nameb[i_s],".png",sep=""))
        # Seurat::DimHeatmap(seu, dims = 1:par_dims, cells = par_cells, balanced = TRUE)
        # ggsave(paste(output_dir,'/step4/info4/',"dimheatplot.",sample_name[i_s],".png",sep=""))
        # FeaturePlot(seu, reduction = par_reduction, features = par_features)
        # ggsave(paste(output_dir,'/step4/info4/',"featureplot",sample_name[i_s],".png",sep=""))
        seu <- RunUMAP(seu, dims = 1:par_dims_umap, n.neighbors =par_n.neighbors)
        Seurat::DimPlot(seu, reduction = "umap")
        ggsave(paste(output_dir,'/step4/figs4/',"dimplot_umap",sample_nameb[i_s],".png",sep=""))
    }
    seu <- MULTIseqDemux(seu, assay = "HTO", quantile = par_quantile, autoThresh = par_autoThresh, maxiter = par_maxiter) 
    write.csv(table(seu$MULTI_ID), paste(output_dir,'/step4/info4/',sample_nameb[i_s],"MULTIseqDemuxHTOcounts.csv",sep=""))
    Idents(seu) <- "MULTI_ID"
    RidgePlot(seu, assay = "HTO", features = rownames(seu[["HTO"]]), group.by = "MULTI_ID", ncol =par_RidgePlot_ncol) 
    ggsave(paste(output_dir,'/step4/figs4/',sample_nameb[i_s],"RidgeplotHTOMSD.png",sep=""))
    Idents(seu) <- "HTO_classification.global"
    VlnPlot(seu, features = "nCount_RNA", pt.size = 0.01, log = TRUE)
    ggsave(paste(output_dir,'/step4/figs4/seu',i_s,"nCounts_demult_groupsMSD.png",sep=""))
    DoHeatmap(seu, features = rownames(seu[["HTO"]]), group.by = "MULTI_ID")
    ggsave(paste(output_dir,'/step4/figs4/seu',i_s,"HeatmapHTOIDMSD.png",sep=""))
    # DotPlot(seu, group.by = "MULTI_ID", features = rownames(seu[["HTO"]])) + theme(axis.text.x = element_text(angle = 90)) # Done 
    DotPlot(seu, features = rownames(seu1[["HTO"]]), group.by = "MULTI_ID")  # sug by MF
    ggsave(paste(output_dir,'/step4/figs4/seu',i_s,"DotPlotHTOIDMSD.png",sep=""))
    multi.names <- unique(seu@meta.data$MULTI_ID)
    Idents(seu)  <- "MULTI_ID"
    levels(seu@meta.data$MULTI_ID)
    for (i in 1:length(old.names)){
        newIdent <- new.names[i]
        names(newIdent) <- old.names[i]
        seu <- RenameIdents(object = seu, newIdent)
    }
    seu[["MULTI_ID_Lables"]] <- Idents(seu)
    seu[['MULTI_classification']] <- NULL
    if (tolower(dropDN)=='yes') {
        print('The following are deleted')
        print(label_dropDN)
            DefaultAssay(seu) <- "RNA"
            aa<-as.character(unique(seu$MULTI_ID))[!as.character(unique(seu$MULTI_ID)) %in% label_dropDN]
            Idents(seu) <- "MULTI_ID"
            seu=subset(seu,idents=aa)
    }  
    if (tolower(Save_RNA)=='yes') {
       mat <- GetAssayData(object = seu, assay = "RNA", slot = "data")
      writeMM(mat,file= paste(output_dir,'/step4/info4/seu',i_s,"_RNA.txt", sep=""))
    }
    if (tolower(Save_metadata)=='yes') {
      write.csv(seu[[]], file = paste(output_dir,'/step4/info4/seu_MetaData',i_s,'.txt', sep=""), quote = TRUE, sep = ",")
    }
    write.csv(colnames(seu[[]]), file= paste(output_dir,'/step4/info4/meta_info',sample_nameb[i_s],".txt", sep=""))
    saveRDS(seu, paste(output_dir,'/step4/objs4/seu',i_s,".rds", sep=""))
    write.csv(table(seu$MULTI_ID), paste(output_dir,'/step4/info4/seu',i_s,"MULTIseqDemuxHTOcounts.csv",sep=""))
}

writeLines(capture.output(sessionInfo()), paste(output_dir,'/step4/info4/sessionInfo.txt', sep=""))
