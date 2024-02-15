#!/usr/bin/env Rscript

####################

args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
seus=args[3]

.libPaths(r_lib_path)

packages<-c('Seurat','ggplot2', 'dplyr')
lapply(packages, library, character.only = TRUE)

loption<-read.table(paste(output_dir,'/job_info/parameters/stepint_par.txt',sep=""),header = FALSE)
lseus<-read.table(seus,header = FALSE)

if(dim(loption)[1]<dim(lseus)[1]){
    print('Check step_inteseu_par.txt, it has smaller length than number of seurats')
      loption[(dim(loption)[1]+1):dim(lseus)[1],1]='yes'
}

print('You are integrating seurat subject')
print('files are')
lseus

seu_list<-list()
library(foreach)
library(doParallel)

numCores <- detectCores()
cl <- makeCluster(numCores-1)
registerDoParallel(cl) 
# i<-1
foreach (i=1:dim(lseus)[1]) %do% {    
    seu<-readRDS(lseus[i,1])
    Seurat::DefaultAssay(seu) <- "RNA"
    if (loption[i,1]=='yes') {
        seu <- Seurat::NormalizeData(seu)
        seu <- Seurat::FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
    }
    seu_list[[i]]<-seu
    rm(seu)
}

seu_anchors <- Seurat::FindIntegrationAnchors(object.list = seu_list)
seu_int <- Seurat::IntegrateData(anchorset = seu_anchors)
Seurat::DefaultAssay(seu_int) <- "integrated"
saveRDS(seu_int, paste(output_dir,"/seu_inetgrated.rds", sep=""))
writeLines(capture.output(sessionInfo()), paste(output_dir,'/sessionInfo_inetgrated.txt', sep=""))
