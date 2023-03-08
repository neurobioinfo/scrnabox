#!/usr/bin/env Rscript

####################
# step2 
## Create seurat object 
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]


.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr')
lapply(packages, library, character.only = TRUE)
# output_dir="/home/samamiri/scratch/Darkgenome_run/pipeline/scrnabox.svn_run/des"
list<-read.csv(paste(output_dir, "/job_output/sample_dir.list",sep=""),header=FALSE)
sample_name<-read.csv(paste(output_dir, "/job_output/sample.list",sep=""),header=FALSE)

library(foreach)
library(doParallel)
numCores <- detectCores()
cl <- makeCluster(numCores-1)
registerDoParallel(cl) 

foreach (i=1:dim(sample_name)[1]) %do% {    
    datadirs <- file.path(list$V1[i],   "ouput_folder","outs","raw_feature_bc_matrix")
    names(datadirs)=sample_name$V1[i]
    sparse_matrix <- Seurat::Read10X(data.dir = datadirs)
    seurat_object <- Seurat::CreateSeuratObject(counts = sparse_matrix$`Gene Expression`)
    seurat_object[['HTO']] = Seurat::CreateAssayObject(counts = sparse_matrix$`Antibody Capture`)
    nam <- paste("seurat_object", sample_name$V1[i], sep = ".")
    assign(nam, seurat_object)
    saveRDS(get(nam),paste(output_dir,'/step2/objs',"/seurat_object.",sample_name$V1[i],".rds", sep=""),compress=TRUE)
    seu<-get(nam)
    seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")
    seu <- subset(seu, subset = percent.mt < 100)
    print(i)
    # png(file = paste(output_dir,'/step2/figs/vioplot',sample_name$V1[i],".png", sep=""))
    Seurat::VlnPlot(seu, group.by= "orig.ident", features = c("nFeature_RNA","nCount_RNA","percent.mt"), pt.size = 0.1,ncol = 3) + NoLegend()
    # dev.off()
    ggsave(paste(output_dir,'/step2/figs2/vioplot',sample_name$V1[i],".png", sep=""))
}


# for (i in 1:dim(sample_name)[1]) {    
#     set.seed(1234)
#     seu<-readRDS(paste(output_dir,'/step2/objs',"/seurat_object.",sample_name$V1[i],".rds", sep=""))
#     seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")
#     seu <- subset(seu, subset = percent.mt < 100)
#     print(i)
#     # png(file = paste(output_dir,'/step2/figs/vioplot',sample_name$V1[i],".png", sep=""))
#     Seurat::VlnPlot(seu, group.by= "orig.ident", features = c("nFeature_RNA","nCount_RNA","percent.mt"), pt.size = 0.1,ncol = 3) + NoLegend()
#     # dev.off()
#     ggsave(paste(output_dir,'/step2/figs/vioplot',sample_name$V1[i],".png", sep=""))
# }

