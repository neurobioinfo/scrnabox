#!/usr/bin/env Rscript

####################
# step4
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
.libPaths(r_lib_path)

set.seed(1234)
packages<-c('Seurat','ggplot2', 'dplyr')
lapply(packages, library, character.only = TRUE)
library('Matrix')

sample_name<-list.files(path = paste(output_dir, "/step3/objs3",sep=""),pattern = "*.rds")

if(length(sample_name)<1) {
   print("You do not have any object from step 2 ")
}

for (i in 1:length(sample_name)) {
  if (!grepl(".rds",tolower(sample_name[i]), fixed = TRUE)){
     print(c(sample_name[i],"is not R rds"))
  }
}   


source(paste(output_dir,'/job_info/parameters/step4_par.txt',sep=""))
library('tidyverse')
# dd<-read.csv(paste(output_dir,'/job_info/parameters/step4_antibody_label.txt',sep=''), header=FALSE, sep="")
# new.names<-str_split(dd[2,],",")[[1]]
# old.names<-str_split(dd[1,],",")[[1]]

old.names<-old_antibody_label
new.names<-new_antibody_label


library(foreach)
library(doParallel)
numCores <- detectCores()
cl <- makeCluster(numCores-1)
registerDoParallel(cl) 


foreach (i_s=1:length(sample_name)) %do% {  
    set.seed(1234)
    seu1<-readRDS(paste(output_dir,'/step3/objs3/',sample_name[i_s], sep=""))
    DefaultAssay(seu1) <- "HTO"
    seu1 <- NormalizeData(seu1, assay = "HTO", normalization.method = "CLR")
    seu1 <- FindVariableFeatures(seu1, selection.method = "mean.var.plot")
    seu1 <- ScaleData(seu1, features = VariableFeatures(seu1))
    seu1 <- MULTIseqDemux(seu1, assay = "HTO", quantile = 0.9, autoThresh = TRUE, maxiter = 5) 
    write.csv(table(seu1$MULTI_ID), paste(output_dir,'/step4/info4/seu',i_s,"MULTIseqDemuxHTOcounts.csv",sep=""))
    Idents(seu1) <- "MULTI_ID"
    RidgePlot(seu1, assay = "HTO", features = rownames(seu1[["HTO"]])[1:6], group.by = "MULTI_ID", ncol =2)
    ggsave(paste(output_dir,'/step4/figs4/seu',i_s,"RidgeplotHTOMSD.png",sep=""))
    Idents(seu1) <- "HTO_classification.global"
    VlnPlot(seu1, features = "nCount_RNA", pt.size = 0.01, log = TRUE)
    ggsave(paste(output_dir,'/step4/figs4/seu',i_s,"nCounts_demult_groupsMSD.png",sep=""))
    DoHeatmap(seu1, features = rownames(seu1[["HTO"]])[1:6], group.by = "MULTI_ID")
    ggsave(paste(output_dir,'/step4/figs4/seu',i_s,"HeatmapHTOIDMSD.png",sep=""))
    ggsave(paste(output_dir,'/step4/figs4/seu',i_s,"DotPlotHTOIDMSD.png",sep=""))
    multi.names <- unique(seu1@meta.data$MULTI_ID)
    Idents(seu1)  <- "MULTI_ID"
    levels(seu1@meta.data$MULTI_ID)
    for (i in 1:length(old.names)){
        newIdent <- new.names[i]
        names(newIdent) <- old.names[i]
        seu1 <- RenameIdents(object = seu1, newIdent)
    }
    seu1[["MULTI_ID_Lables"]] <- Idents(seu1)
    seu1[['MULTI_classification']] <- NULL
    if (tolower(dropDN)=='yes') {
        print('The following are deleted')
        print(label_dropDN)
            DefaultAssay(seu1) <- "RNA"
            aa<-as.character(unique(seu1$MULTI_ID))[!as.character(unique(seu1$MULTI_ID)) %in% label_dropDN]
            Idents(seu1) <- "MULTI_ID"
            seu1=subset(seu1,idents=aa)
    }  
    if (tolower(Save_RNA)=='yes') {
       mat <- GetAssayData(object = seu1, assay = "RNA", slot = "data")
      writeMM(mat,file= paste(output_dir,'/step4/info4/seu',i_s,"_RNA.txt", sep=""))
    }
    if (tolower(Save_metadata)=='yes') {
      write.csv(seu1[[]], file = paste(output_dir,'/step4/info4/seu_MetaData',i_s,'.txt', sep=""), quote = TRUE, sep = ",")
    }
    write.csv(colnames(seu1[[]]), file= paste(output_dir,'/step4/info4/meta_info_seu',i_s,".txt", sep=""))
    saveRDS(seu1, paste(output_dir,'/step4/objs4/seu',i_s,"dem.rds", sep=""))
    write.csv(table(seu1$MULTI_ID), paste(output_dir,'/step4/info4/seu',i_s,"MULTIseqDemuxHTOcounts.csv",sep=""))
}

writeLines(capture.output(sessionInfo()), paste(output_dir,'/step4/info4/sessionInfo.txt', sep=""))
