###########################################################################################
## Figure 4A
##########################################################################################
## load parameters
output_dir = "/path/to/output/directory"
r_lib_path = "/path/to/R/library"

## load library
.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr', 'foreach', 'doParallel','Matrix', 'ggpubr')
lapply(packages, library, character.only = TRUE)
library(RColorBrewer)

## load integrated Seurat object
seu_int<-readRDS('/path/to/step6/objs6/seu_step6.rds')

## plot
Seurat::DimPlot(seu_int, group.by = 'integrated_snn_res.0.2') +
theme(axis.text.x = element_text(size = 16),
axis.text.y = element_text(size = 16),
axis.title = element_text(size = 16, face = "bold"),
legend.text = element_text(size = 16),
legend.title = element_blank()) +
xlab("UMAP1") + ylab("UMAP2") + 
scale_color_manual(values = rev(c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFD92F", "#B15928", "lightgrey",  "#666666"))) +
ggtitle("")
ggsave(paste(output_dir,'/manuscript_figures/clustered_UMAP.pdf',sep="")) 


##########################################################################################
## Figure 4B
##########################################################################################
gene_sets <- read.delim('/path/to/step7/info7/marker/top_sel.csv', header = T, sep = ",", na.strings=c("","NA"))
seu_int<-readRDS('/path/to/step6/objs6/seu_step6.rds')

heat_map<-DoHeatmap(seu_int, features = gene_sets$gene, size=0, angle =90, group.bar.height = 0.02, group.by = 'integrated_snn_res.0.2', 
group.colors = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFD92F", "#B15928", "lightgrey",  "#666666")) +
theme(axis.text.y = element_text(size = 20),
axis.text.x = element_blank(),
axis.title = element_text(size = 20, face = "bold"),
legend.text = element_text(size = 20),
legend.title = element_blank()) +
ggtitle("") 
ggsave(file = paste(output_dir,'/manuscript_figures/heatclust.pdf',sep=""), dpi = 300, height = 11, width = 11, unit = 'in' )


##########################################################################################
## Figure 4C
##########################################################################################
## output directly from pipeline

##########################################################################################
## Figure 4D-E
##########################################################################################
## load Seurat object
seu_int<-readRDS('/path/to/step7/objs7/seu_step7.rds')

## load parameters
par_select_features_list= c("MOBP", "VCAN", "AQP4","FOXJ1", "CD74", "CLDN5", "GFRB", "SLC17A6", "GAD2", "GRIK1", "TH", "CADPS2")
DefaultAssay(seu_int) <- "RNA"

## dotplot
dot_plt <- DotPlot(seu_int, features = par_select_features_list, group.by = 'clustering_3') + 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size =16),
axis.text.y = element_text(size =16),
axis.title = element_text(face="bold", size =16),
legend.text = element_text( size =16),
legend.title = element_text( size =16))
ggsave(file = paste(output_dir,'/manuscript_figures/annotated_dotplot.pdf',sep=""))

## load parameters
par_select_features_list= c("CD74")

# violin plot
vln_plt <- VlnPlot(seu_int, features = par_select_features_list, group.by = 'integrated_snn_res.0.2', pt.size = 0) +
theme(axis.text.x = element_text(size =16),
axis.text.y = element_text(size =16),
axis.title = element_text(face="bold", size =16),
legend.text = element_text( size =16),
legend.title = element_text( size =16)) + 
scale_fill_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFD92F", "#B15928", "lightgrey",  "#666666")) +
ggtitle("") +
ylab("CD74 Expression Level")
ggsave(file = paste(output_dir,'/manuscript_figures/annotated_list_violin_plot.pdf', sep=''))


##########################################################################################
## Figure 4F
##########################################################################################
## output directly from pipeline


##########################################################################################
## Figure 4G
##########################################################################################
## load Seurat objects
seu_int<-readRDS('/path/to/step7/objs7/seu_step7.rds')
reference0<-readRDS("/path/to/kamath/reference.rds")

reference0<- NormalizeData(reference0)
reference0 <- FindVariableFeatures(reference0)
reference0<- ScaleData(reference0)
reference0 <- RunPCA(object = reference0, assay = "RNA", npcs = 10)
reference0 <- RunUMAP(reference0, dims = 1:10, reduction = "pca", return.model = TRUE)

## plot reference
p1 <- DimPlot(reference0, reduction = "umap", group.by = "Cell_Type", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) + NoLegend() + ggtitle("Reference annotations") +
theme(axis.text.x = element_text(size =10),
axis.text.y = element_text(size =10),
axis.title = element_text(face="bold", size =10),
legend.text = element_text( size =10),
plot.title = element_text( size =10)) + 
scale_color_manual(values = c("astro" = "#33A02C",
"endo"="#FF7F00",
"olig"="#A6CEE3",
"mg"="#E31A1C",
"da"="#1F78B4",
"nonda"="#FB9A99",
"opc"="#FDBF6F"))+
xlab("UMAP1") + ylab("UMAP2")

## plot query
p2 <- DimPlot(seu_int, reduction = "umap", group.by = "predicted.id", label = FALSE, label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels") +
theme(axis.text.x = element_text(size =10),
axis.text.y = element_text(size =10),
axis.title = element_text(face="bold", size =10),
legend.text = element_text( size =10),
plot.title = element_text( size =10)) + 
scale_color_manual(values = c("astro" = "#33A02C",
"endo"="#FF7F00",
"olig"="#A6CEE3",
"mg"="#E31A1C",
"da"="#1F78B4",
"nonda"="#FB9A99",
"opc"="#FDBF6F"))+
xlab("UMAP1") + ylab("UMAP2")

## plot both and save
p1 + p2 
ggsave(file = paste(output_dir,'/manuscript_figures/reference.pdf', sep=''), height = 4, width = 8)


##########################################################################################
## Figure 4E
##########################################################################################
## load Seurat object
seu_int<-readRDS('/path/to/step7/objs7/seu_step7.rds')

## plot
p2 <- DimPlot(seu_int, reduction = "umap", group.by = "clustering_3", label = FALSE, label.size = 3, repel = TRUE) + NoLegend() + ggtitle("") +
theme(axis.text.x = element_text(size =10),
axis.text.y = element_text(size =10),
axis.title = element_text(face="bold", size =10),
legend.text = element_text( size =10),
plot.title = element_text( size =10)) + 
scale_colour_manual(values = c("astro" = "#33A02C",
"endo"="#FF7F00",
"Olig"="#A6CEE3",
"mg"="#E31A1C",
"opc"="#FDBF6F",
"inhib"= "#B2DF8A",
"excit" ="#FB9A99",  
"peri" = "#CAB2D6",
"ependymal"= "#6A3D9A",
"GABA" = "#FFD92F",
"DaN" = "#1F78B4",
"CADPS2" = "#B15928")) +
xlab("UMAP1") + ylab("UMAP2")
ggsave(file = paste(output_dir,'/manuscript_figures/Annotate.pdf', sep=''), height = 4, width =4)




