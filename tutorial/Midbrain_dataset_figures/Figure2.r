##########################################################################################
## Figure 2A
##########################################################################################
# Output directly from pipeline

##########################################################################################
## Figure 2B-D
##########################################################################################
## load parameters
output_dir = "/path/to/output/directory"
r_lib_path = "/path/to/R/library"

## load library
.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr', 'foreach', 'doParallel','Matrix', 'ggpubr')
lapply(packages, library, character.only = TRUE)

## load parameters
source(paste(output_dir,'/job_info/parameters/step2_par.txt',sep=""))
           
### import Seurat  object 
seu<-readRDS('/path/to/step2/objs2/Control1.rds')

## PCA cell cycle score
seu <- RunPCA(seu, features = c(cc.genes$g2m.genes, cc.genes$s.genes),raster = FALSE)
DimPlot(seu, group.by = "Phase", pt.size =1, raster =FALSE) +
theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) +
        xlab("PC1") + ylab("PC2") +ggtitle("")
ggsave(paste(output_dir,'/manuscript_figures/','ccpca.pdf', sep=""))

## print violin plot for QC metrics
## nefeature
n_feature <- Seurat::VlnPlot(seu, group.by= "orig.ident", features = "nFeature_RNA", pt.size = 0.001,raster = FALSE) + NoLegend()  + 
ggtitle("") + 
xlab("Control1") + 
ylab("Unique RNA transcripts") +
theme(axis.text.x = element_blank(),
axis.text.y = element_text(size = 22),
axis.title = element_text(size = 22, face = "bold")) +
scale_y_continuous(limits = c(0,12000), breaks = c(0,2500,5000,7500,10000,12500,15000))
ggsave(paste(output_dir,'/manuscript_figures/','nfeature_',".pdf", sep=""))

## ncount
n_count <- Seurat::VlnPlot(seu, group.by= "orig.ident", features = "nCount_RNA", pt.size = 0.01,raster = FALSE) + NoLegend() + 
ggtitle("") + 
xlab("Control1") + 
ylab("Total RNA transcripts") +
theme(axis.text.x = element_blank(),
axis.text.y = element_text(size = 22),
axis.title = element_text(size = 22, face = "bold")) +
scale_y_continuous(limits = c(0,250000),breaks = c(0,25000, 50000, 75000, 100000,125000, 150000, 175000, 200000, 225000, 250000))
ggsave(paste(output_dir,'/manuscript_figures/','total_',".pdf", sep=""))

## percent mt
mito <- Seurat::VlnPlot(seu, group.by= "orig.ident", features = "percent.mt", pt.size = 0.01,raster = FALSE) + NoLegend() + 
ggtitle("") + 
xlab("Control1") + 
ylab("Percent mitochondrial RNA") +
theme(axis.text.x = element_blank(),
axis.text.y = element_text(size = 22),
axis.title = element_text(size = 22, face = "bold")) +
scale_y_continuous(limits = c(0,100), breaks = c(0,12.5,25,37.5,50,62.5,75,87.5,100))
ggsave(paste(output_dir,'/manuscript_figures/','mito_',".pdf", sep=""))

## percent ribo
library(ggbreak)
ribo <- Seurat::VlnPlot(seu, group.by= "orig.ident", features = "percent.ribo", pt.size = 0.01,raster = FALSE) + NoLegend() + 
ggtitle("") + 
xlab("Control1") + 
ylab("Percent ribosomal RNA") +
theme(axis.text.x = element_blank(),
axis.text.y = element_text(size = 22),
axis.title = element_text(size = 22, face = "bold")) +
scale_y_continuous(limits = c(0,100), breaks = c(0,2.5,5,7.5,10,11,100)) +
scale_y_break(c(11, 100))
ggsave(paste(output_dir,'/manuscript_figures/','ribo_' ,".pdf", sep=""))

