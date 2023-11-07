###########################################################################################
## Figure 6A-D
##########################################################################################
## load parameters
output_dir = "/path/to/output/directory"
r_lib_path = "/path/to/R/library"

## load library
.libPaths(r_lib_path)

## set seed for reproducibility
set.seed(1234)

## load library
packages<-c('Seurat','ggplot2', 'dplyr','foreach', 'doParallel','Matrix','tidyverse')
lapply(packages, library, character.only = TRUE)

## load list of existing Seurat objects from Step 3
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

## load parameters
source(paste(output_dir,'/job_info/parameters/step4_par.txt',sep=""))

## load old and new antibody labels
old.names<-par_old_antibody_label
new.names<-par_new_antibody_label

## detect number of available cores
numCores <- detectCores()
cl <- makeCluster(numCores-1)
registerDoParallel(cl) 

## perform doublet detection and demultiplexing for each Seurat object
i_s=1
#foreach (i_s=1:length(sample_name)) %do% {  
set.seed(1234)
seu<-readRDS(paste(output_dir,'/step3/objs3/',sample_name[i_s], sep=""))
DefaultAssay(seu) <- "HTO"

## normalize and scale HTO assay
seu <- NormalizeData(seu, assay = "HTO", normalization.method = par_normalization.method,scale.factor =par_scale.factor)
seu <- FindVariableFeatures(seu, selection.method = par_selection.method, nfeatures = par_nfeatures)
seu <- ScaleData(seu)

## linear dimensional reduction on HTO assay (optional)
if (tolower(par_dimensionality_reduction)=='yes'){
seu <- RunPCA(seu, npcs = par_npcs_pca, verbose = FALSE)
DimPlot(seu, reduction = "pca")
ggsave(paste(output_dir,'/step4/figs4/',sample_nameb[i_s],"_HTO_dimplot_pca",".pdf",sep=""))
seu <- RunUMAP(seu, dims = 1:par_dims_umap, n.neighbors =par_n.neighbors)
Seurat::DimPlot(seu, reduction = "umap")
ggsave(paste(output_dir,"_HTO_dimplot_umap", ".pdf",sep=""))
}

## doublet detection and demultiplexing with MULTIseqDemux
seu <- MULTIseqDemux(seu, assay = "HTO", quantile = par_quantile, autoThresh = par_autoThresh, maxiter = par_maxiter) 
write.csv(table(seu$MULTI_ID), paste(output_dir,'/step4/info4/',sample_nameb[i_s],"_MULTIseqDemuxHTOcounts.csv",sep=""))

## rename droplet identities
Idents(seu) <- "MULTI_ID"

## rename identies with new antibody label
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

## print ridge plot
RidgePlot(seu, assay = "HTO", features = rownames(seu[["HTO"]]), group.by = "MULTI_ID_Lables", ncol =4, cols = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F","#FF7F00", "#CAB2D6", "#6A3D9A"))
ggsave(paste(output_dir,"/Ridgeplot_HTO_MSD.pdf",sep=""),dpi = 300, height = 8, width = 16, unit = 'in' )

## print violin plot for ncount_RNA
Idents(seu) <- "HTO_classification.global"
VlnPlot(seu, features = "nCount_RNA", pt.size = 0.01, log = TRUE, group.by = "MULTI_ID_Lables") +
theme(axis.text.x = element_text(size =16),
axis.text.y = element_text(size =16),
axis.title = element_text(face="bold", size =16),
legend.text = element_text( size =16),
plot.title = element_text( size =16),
legend.position = "none") + 
scale_fill_manual(values = c( "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F","#FF7F00", "#CAB2D6", "#6A3D9A")) +
xlab("Droplet assignment") + ylab("Total RNA transcripts") +
ggtitle("")
ggsave(paste(output_dir,"/Violin.pdf",sep=""), height = 7, width = 6)

## print heatmap
DoHeatmap(seu, features = rownames(seu[["HTO"]]), group.by = "MULTI_ID_Lables", size = 5, group.colors = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F","#FF7F00", "#CAB2D6", "#6A3D9A")) +
theme(axis.text.y = element_text(size = 14),
axis.text.x = element_blank(),
axis.title = element_text(size = 14, face = "bold"),
legend.text = element_text(size = 14),
legend.title = element_blank(),
legend.position = "right") +
ggtitle("") +
guides(colour = "none")
ggsave(paste(output_dir,"/heat.pdf",sep=""), width = 6, height = 6)

## print dotplot
DotPlot(seu, features = rownames(seu[["HTO"]]), group.by = "MULTI_ID_Lables") + 
theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1, size = 16),
axis.text.y = element_text(size =16),
axis.title = element_text(face="bold", size =16),
legend.text = element_text( size =16),
plot.title = element_text( size =16)) + 
xlab("Barcodes") + ylab("Droplet assignment") +
ggtitle("") 
ggsave(paste(output_dir,"/dot.pdf",sep=""), width = 7, height = 6) 
