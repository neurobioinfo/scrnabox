#!/usr/bin/env Rscript

####################
# step4
## Create seurat object 
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
.libPaths(r_lib_path)
# source(args[3])
set.seed(1234)
packages<-c('Seurat','ggplot2', 'dplyr')
lapply(packages, library, character.only = TRUE)

# output_dir="/home/samamiri/scratch/Darkgenome_run/pipeline/scrnabox.svn_run/des"
list<-read.csv(paste(output_dir, "/job_output/sample_dir.list",sep=""),header=FALSE)
sample_name<-read.csv(paste(output_dir, "/job_output/sample.list",sep=""),header=FALSE)

library('tidyverse')
dd<-read.csv(paste(output_dir,'/job_output/step4_par.txt',sep=''), header=FALSE, sep="")
new.names<-str_split(dd[2,],",")[[1]]
old.names<-str_split(dd[1,],",")[[1]]

library(foreach)
library(doParallel)
numCores <- detectCores()
cl <- makeCluster(numCores-1)
registerDoParallel(cl) 


foreach (i_s=1:dim(sample_name)[1]) %do% {  
    set.seed(1234)
    seu1<-readRDS(paste(output_dir,'/step3/objs',"/seu",i_s,".rds", sep=""))
    DefaultAssay(seu1) <- "HTO"
    seu1 <- NormalizeData(seu1, assay = "HTO", normalization.method = "CLR")
    seu1 <- FindVariableFeatures(seu1, selection.method = "mean.var.plot")
    seu1 <- ScaleData(seu1, features = VariableFeatures(seu1))
    seu1 <- MULTIseqDemux(seu1, assay = "HTO", quantile = 0.9, autoThresh = TRUE, maxiter = 5) 
    write.csv(table(seu1$MULTI_ID), paste(output_dir,'/step4/objs/seu',i_s,"MULTIseqDemuxHTOcounts.csv",sep=""))
    Idents(seu1) <- "MULTI_ID"
    # png(paste(output_dir,'/step4/figs/seu',i_s,"RidgeplotHTOMSD.png",sep=""))
    RidgePlot(seu1, assay = "HTO", features = rownames(seu1[["HTO"]])[1:6], group.by = "MULTI_ID", ncol =2)
    ggsave(paste(output_dir,'/step4/figs4/seu',i_s,"RidgeplotHTOMSD.png",sep=""))
    # dev.off()
    Idents(seu1) <- "HTO_classification.global"
    # png(paste(output_dir,'/step4/figs/seu',i_s,"nCounts_demult_groupsMSD.png",sep=""))
    VlnPlot(seu1, features = "nCount_RNA", pt.size = 0.01, log = TRUE)
    ggsave(paste(output_dir,'/step4/figs4/seu',i_s,"nCounts_demult_groupsMSD.png",sep=""))
    # dev.off()
    # png(paste(output_dir,'/step4/figs/seu',i_s,"HeatmapHashIDMSD.png",sep=""))
    DoHeatmap(seu1, features = rownames(seu1[["HTO"]])[1:6], group.by = "MULTI_ID")# done
    ggsave(paste(output_dir,'/step4/figs4/seu',i_s,"HeatmapHashIDMSD.png",sep=""))
    # dev.off()
    # png(paste(output_dir,'/step4/figs/seu',i_s,"DotPlotHashIDMSD.png",sep=""))
    DotPlot(seu1, group.by = "MULTI_ID", features = rownames(seu1[["HTO"]])[1:6]) + theme(axis.text.x = element_text(angle = 90)) # Done 
    ggsave(paste(output_dir,'/step4/figs4/seu',i_s,"DotPlotHashIDMSD.png",sep=""))
    # dev.off()
    multi.names <- unique(seu1@meta.data$MULTI_ID)
    Idents(seu1)  <- "MULTI_ID"
    # old.names <- c("B0251-TotalSeqB","B0252-TotalSeqB","B0253-TotalSeqB","B0254-TotalSeqB","B0255-TotalSeqB","B0256-TotalSeqB","Doublet","Negative")
    # new.names <- c("AIW002","SNCA-A53T","GBA-KO","Parkin-KO","PINK1-KO","SNCA-KO","Doublet","Negative")
    levels(seu1@meta.data$MULTI_ID)
    for (i in 1:length(old.names)){
    newIdent <- new.names[i]
    names(newIdent) <- old.names[i]
    seu1 <- RenameIdents(object = seu1, newIdent)
    }
    seu1[["Lables"]] <- Idents(seu1)
    # table(seu1$MULTI_ID)
    # table(seu1$MULTI_classification) # show all the matches doublets
    # write.csv(table(seu1$MULTI_ID), paste(output_dir,'/step4/objs/seu',i_s,"MULTIseqDemuxHTOcounts.csv",sep=""))
    saveRDS(seu1, paste(output_dir,'/step4/objs/seu',i_s,"dem.rds", sep=""))
}