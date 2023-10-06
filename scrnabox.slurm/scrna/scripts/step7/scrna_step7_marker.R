#!/usr/bin/env Rscript

###############################################################################
# step7 -- cluster annotation
###############################################################################

## set sample ID metadata column -- this is standard and does not require parameter modification
par_level_genotype <- "Sample_ID"

## load parameters
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
pipeline_home=args[3]


## load library
.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr', 'xlsx', 'Matrix')
lapply(packages, library, character.only = TRUE)

## load parameters
source(paste(output_dir,'/job_info/parameters/step7_par.txt',sep=""))

###############################################################################
# step7 -- find marker
###############################################################################
if (tolower(par_run_find_marker)=='yes') {

## load existing Seurat objects
sample_name<-list.files(path = paste(output_dir, "/step6/objs6",sep=""),pattern = "*.rds")
seu_int<-readRDS(paste(output_dir,'/step6/objs6/',sample_name, sep=''))

## set cell identity to the clustering resolution defined by the user
Idents(seu_int) <- par_level_cluster

## create directories for marker annotation method
## figures 
OUT_DIR_figs <- paste(output_dir,"/step7/figs7",sep='') 
OUT_dir_figs_marker <- paste(OUT_DIR_figs,"/marker/",sep='') 
dir.create(OUT_dir_figs_marker)
## info
OUT_DIR_info <- paste(output_dir,"/step7/info7",sep='') 
OUT_dir_info_marker <- paste(OUT_DIR_info,"/marker/",sep='') 
dir.create(OUT_dir_info_marker)

## set output directory
PWD=OUT_dir_info_marker
setwd(PWD)

## identify cluster marker
ClusterMarkers <- FindAllMarkers(seu_int, only.pos = TRUE)

## identify top n marker genes for each cluster -- the number of top marker genes for each cluster is defined by the user
top5 <- ClusterMarkers %>% group_by(cluster) %>% top_n(n=par_top_sel, wt = avg_log2FC)
write.csv(top5,"top_sel.csv")

## print cluster information
if(file.exists(paste(output_dir,'/step7/info7/marker/', "cluster_whole.xlsx",sep=''))){
    file.remove(paste(output_dir,'/step7/info7/marker/', "cluster_whole.xlsx",sep=''))
}
for (i in  sort(unlist(unique(seu_int[[par_level_cluster]])))) {
  N1.c0 <- ClusterMarkers %>% filter(cluster == i & avg_log2FC > 0)
  genes <- N1.c0$gene
  write.xlsx(genes, file="cluster_just_genes.xlsx",sheetName=paste0('cluster',i), row.names=FALSE,append=TRUE)
  write.xlsx(N1.c0, file="cluster_whole.xlsx",sheetName=paste0('cluster',i), row.names=FALSE,append=TRUE)
}

## print heatmap of top marker genes
heat_map<-DoHeatmap(seu_int, features = top5$gene, size=3, angle =90, group.bar.height = 0.02, group.by = par_level_cluster)
ggsave(file = paste(OUT_dir_figs_marker,'heatmap.pdf', sep=''), dpi = 300, height = 11, width = 12, unit = 'in' )

## save cluster marker information
saveRDS(ClusterMarkers, paste(OUT_dir_info_marker,"ClusterMarkers.rds", sep=""))
write.csv(ClusterMarkers,paste(OUT_dir_info_marker,"ClusterMarkers.csv",sep=''))

## save session info
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step7/info7/sessionInfo_find_marker.txt', sep=""))
if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}
}

###############################################################################
# step7 -- enrichR marker GSEA
###############################################################################

if (tolower(par_run_enrichR)=='yes') {

## load libraries
packages<-c('Seurat','ggplot2', 'dplyr','xlsx','enrichR')
lapply(packages, library, character.only = TRUE)

## set working directory to step7 directory
PWD=paste(output_dir,'/step7/', sep='')
setwd(PWD)

## read in the cluster marker csv file csv file
cluster_marker <- read.delim(paste(output_dir,'step7/info7/marker/ClusterMarkers.csv', sep = ""), header = T, sep = ",") 
  
## create directories based on the clusters present in the csv file
  dir.create("enrichR") 
  for (i in unique(cluster_marker$cluster)) {
    dir.create(paste0(PWD,"/enrichR/",par_level_cluster,"/clust",i))
  }

  for (i in unique(cluster_marker$cluster)) {
    setwd(PWD)
    setwd(paste0('./enrichR/',par_level_cluster,'clust',i))
    N1.c0 <- cluster_marker %>% filter(cluster == i & avg_log2FC > 0)
    genes <- N1.c0$gene
    N1.c0.Er <- enrichr(genes, databases = par_db)
    if(is.null(N1.c0.Er)) next
    plotEnrich(N1.c0.Er[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
    ggsave(file = "plotenrich1.pdf")
    plotEnrich(N1.c0.Er[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
    ggsave(file = "plotenrich2.pdf")
    plotEnrich(N1.c0.Er[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
    ggsave(file = "plotenrich3.pdf")
    ### You don't need to make the tables or else this should all be saved
    N1.Er.genes.1 <- N1.c0.Er[[1]] %>% select(Term, Genes, Combined.Score)
    write.csv(N1.Er.genes.1,"Er.genes.1.csv")
    N1.Er.genes.2 <- N1.c0.Er[[2]] %>% select(Term, Genes, Combined.Score)
    write.csv(N1.Er.genes.2,"Er.genes.2.csv")
    N1.Er.genes.3 <- N1.c0.Er[[3]] %>% select(Term, Genes, Combined.Score)
    write.csv(N1.Er.genes.3,"Er.genes.3.csv")
  }

## write session info
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step7/info7/sessionInfo_enrich.txt', sep=""))

if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}
}


###############################################################################
# step7 -- module score
###############################################################################

if (tolower(par_run_module_score)=='yes') {

## load existing Seurat objects
sample_name<-list.files(path = paste(output_dir, "/step6/objs6",sep=""),pattern = "*.rds")
seu_int<-readRDS(paste(output_dir,'/step6/objs6/',sample_name, sep=''))

## set cell identity to the clustering resolution defined by the user
Idents(seu_int) <- par_level_cluster

## create directories for module score
## figures 
OUT_DIR_figs <- paste(output_dir,"/step7/figs7",sep='') 
OUT_dir_figs_module_score <- paste(OUT_DIR_figs,"/module_score/",sep='') 
dir.create(OUT_dir_figs_module_score)
## info
OUT_DIR_info <- paste(output_dir,"/step7/info7",sep='') 
OUT_dir_info_module_score <- paste(OUT_DIR_info,"/module_score/",sep='') 
dir.create(OUT_dir_info_module_score)

## set output directory
PWD=OUT_dir_info_module_score
setwd(PWD)

## compute module score
source(paste(pipeline_home,'/general_codes/module_score.R',sep=''))

## save session info
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step7/info7/sessionInfo_module_score.txt', sep=""))
if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}
}

###############################################################################
# step7 -- reference based annotation
###############################################################################

if (tolower(par_run_reference)=='yes') {

################### import the right Seurat object ###################
## load name of existing Seurat objects
sample_name<-list.files(path = paste(output_dir, "/step6/objs6",sep=""),pattern = "*.rds")

if(file.exists(paste(output_dir,'/step7/objs7/','seu_step7.rds', sep = ""))){
    seu_int<-readRDS(paste(output_dir,'/step7/objs7/','seu_step7.rds', sep=''))
}else{
    seu_int<-readRDS(paste(output_dir,'/step6/objs6/',sample_name, sep=''))
}
################### ############################## ###################

## set output directory
PWD=paste(output_dir,'/step7/objs7/', sep='')
setwd(PWD)

## create reference-based directory
#figures
OUT_DIR_figs <- paste(output_dir,"/step7/figs7",sep='') 
OUT_dir_figs_reference <- paste(OUT_DIR_figs,"/reference_based_annotation/",sep='') 
dir.create(OUT_dir_figs_reference)
#info
OUT_DIR_info <- paste(output_dir,"/step7/info7",sep='') 
OUT_dir_info_reference <- paste(OUT_DIR_info,"/reference_based_annotation/",sep='') 
dir.create(OUT_dir_info_reference)

## set user defined clustering resolution
seu_int <- SetIdent(seu_int, value = par_level_cluster)

## load reference Seurat object
reference0 <-readRDS(par_reference)
DefaultAssay(reference0) <- "RNA" ## new code

## load parallelization parameters
options(future.globals.maxSize = par_futureglobalsmaxSize)

# perform standard preprocessing on reference object
reference0<- NormalizeData(reference0)
reference0 <- FindVariableFeatures(reference0)
reference0<- ScaleData(reference0)
reference0 <- RunPCA(object = reference0, assay = "RNA", npcs = par_FindTransferAnchors_dim)

## find transfer anchors between reference and query Seurat objects
transfer.anchors <- FindTransferAnchors(reference = reference0, query = seu_int, dims = 1:par_FindTransferAnchors_dim, reference.reduction = "pca")

## add reference-based annotations to the qeury object
eval(parse(text = paste('predictions <- TransferData(anchorset = transfer.anchors, refdata = reference0$',par_level_celltype ,',dims = 1:',par_FindTransferAnchors_dim,')', sep='')))
seu_int <- AddMetaData(object = seu_int, metadata = predictions)

# Add metadata column for reference object
seu_int$temp_temp_2 <- seu_int@meta.data$predicted.id
name_meta <- names(seu_int@meta.data) 
length <- length(name_meta)
name_meta[length] <- paste(par_reference_name, "_predictions", sep = "")
names(seu_int@meta.data) <- name_meta

## save query Seurat object with reference annotation predicitions
saveRDS(seu_int,paste(output_dir,'/step7/objs7','/seu_step7.rds', sep=''))

## save metadata information
write.csv(colnames(seu_int[[]]), file= paste(output_dir,'/step7/info7/meta_info_seu_step7',".txt", sep=""))

## Print a umap projection showing the predicted cell types on the query object 
reference0 <- RunUMAP(reference0, dims = 1:par_FindTransferAnchors_dim, reduction = "pca", return.model = TRUE)
seu_int <- MapQuery(anchorset = transfer.anchors, reference = reference0, query = seu_int,
    refdata = list(celltype = par_level_celltype), reference.reduction = "pca", reduction.model = "umap")
p1 <- DimPlot(reference0, reduction = "umap", group.by = par_level_celltype, label = TRUE, label.size = 3,repel = TRUE, raster=FALSE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(seu_int, reduction = "umap", group.by = "predicted.id", label = TRUE, label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2
ggsave(file = paste(OUT_dir_figs_reference,par_reference_name,'_UMAP_transferred_labels.pdf', sep=''), dpi = 300, height = 7, width = 14, unit = 'in')

## print summary table
df <- data.frame(seu_int@meta.data)
cluster_list <- list(par_level_cluster)
keep_list <- append(cluster_list, "predicted.id")
df <- df[, (colnames(df) %in% keep_list)]
df_summary <- table(df[,1], df[,2])
df_summary <- data.frame(df_summary)
colnames(df_summary) <- c("cluster", "cell_type", "number_of_cells")
    
if(file.exists(paste(output_dir,'/step7/info7/reference_based_annotation/',par_reference_name, '_prediction_summary.xlsx', sep=""))){
file.remove(paste(output_dir,'/step7/info7/reference_based_annotation/',par_reference_name, '_prediction_summary.xlsx', sep=""))
}
    
write.xlsx(df_summary, file=paste(output_dir,'/step7/info7/reference_based_annotation/',par_reference_name, '_prediction_summary.xlsx', sep=""),row.names=FALSE)

## save RNA expression matrix
if (tolower(par_save_RNA)=='yes') {
    mat <- GetAssayData(object = seu_int, assay = "RNA", slot = "data")
    writeMM(mat,file= paste(output_dir,'/step7/info7/seu',"_RNA.txt", sep=""))
}

## save metadata 
if (tolower(par_save_metadata)=='yes') {
    write.csv(seu_int[[]], file = paste(output_dir,'/step7/info7/seu_MetaData.txt', sep=""), quote = TRUE, sep = ",")
}

## save session info
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step7/info7/sessionInfo_reference.txt', sep=""))
if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}
}

###############################################################################
# step7 -- visualize marker
###############################################################################

if (tolower(par_run_visualize_markers)=='yes') {

## load existing Seurat objects
sample_name<-list.files(path = paste(output_dir, "/step6/objs6",sep=""),pattern = "*.rds")
seu_int<-readRDS(paste(output_dir,'/step6/objs6/',sample_name, sep=''))

## set cell identity to the clustering resolution defined by the user
Idents(seu_int) <- par_level_cluster

## create directories for visualize features annotation method
## figures 
OUT_DIR_figs <- paste(output_dir,"/step7/figs7",sep='') 
OUT_dir_figs_visualize_features <- paste(OUT_DIR_figs,"/visualize_features/",sep='') 
dir.create(OUT_dir_figs_visualize_features)

## set output directory
PWD=OUT_dir_figs_visualize_features
setwd(PWD)

###################
## visulize select features with a list in the parameters
###################
if (exists("par_select_features_list")) {
#define the features
select_features <- par_select_features_list

#violin plot
vln_plt <- VlnPlot(seu_int, features = select_features, group.by = par_level_cluster, pt.size = 0)
ggsave(file = paste(OUT_dir_figs_visualize_features,'list_violin_plot.pdf', sep=''), dpi = 300, height = 15, width = 15, unit = 'in')

#feature plot
Feature_plt <- FeaturePlot(seu_int, features = select_features, raster = FALSE)
ggsave(file = paste(OUT_dir_figs_visualize_features ,'list_feature_plot.pdf', sep=''), dpi = 300, height = 15, width = 20, unit = 'in')

#dotplot
dot_plt <- DotPlot(seu_int, features = select_features, group.by = par_level_cluster) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(file = paste(OUT_dir_figs_visualize_features ,'list_dot_plot.pdf', sep=''))
}

###################
## visulize select features with a csv
###################

if (exists("par_select_features_csv")) {
gene_sets <- read.delim(par_select_features_csv, header = T, sep = ",", na.strings=c("","NA")) 

#convert user inputed dataframe to named list
gene_lists <- list()                   # Create empty list
for(i in 1:ncol(gene_sets)) {             # Using for-loop to add columns to list
  gene_lists[[i]] <- gene_sets[ , i]
}
names(gene_lists) <- colnames(gene_sets)#set names of list

## remove NA from list
gene_lists <- lapply(gene_lists, function(x) x[!is.na(x)])

for (i in 1:length(gene_lists)) {
    # Add module scores for each gene list
    #try(
            #violin plot
            vln_plt <- VlnPlot(seu_int, features = gene_lists[[i]], group.by = par_level_cluster, pt.size = 0)
            ggsave(file = paste(OUT_dir_figs_visualize_features, names(gene_lists[i]),'_violin_plot.pdf', sep=''), dpi = 300, height = 15, width = 15, unit = 'in')

            #feature plot
            Feature_plt <- FeaturePlot(seu_int, features = gene_lists[[i]], raster = FALSE)
            ggsave(file = paste(OUT_dir_figs_visualize_features, names(gene_lists[i]) ,'_feature_plot.pdf', sep=''), dpi = 300, height = 15, width = 20, unit = 'in')

            #dotplot
            dot_plt <- DotPlot(seu_int, features = gene_lists[[i]], group.by = par_level_cluster) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
            ggsave(file = paste(OUT_dir_figs_visualize_features, names(gene_lists[i]),'_dot_plot.pdf', sep=''))
       #)
    }
}
 

## save session info
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step7/info7/sessionInfo_viualize_features.txt', sep=""))
if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}
}


###############################################################################
# step7 -- annotate
###############################################################################

if (tolower(par_run_annotate)=='yes') {

################### import the right Seurat object ###################
## load name of existing Seurat objects
sample_name<-list.files(path = paste(output_dir, "/step6/objs6",sep=""),pattern = "*.rds")

if(file.exists(paste(output_dir,'/step7/objs7/','seu_step7.rds', sep = ""))){
    seu_int<-readRDS(paste(output_dir,'/step7/objs7/','seu_step7.rds', sep=''))
}else{
    seu_int<-readRDS(paste(output_dir,'/step6/objs6/',sample_name, sep=''))
}
################### ############################## ###################

## set cell identity to the clustering resolution defined by the user
Idents(seu_int) <- par_annotate_resolution

## create directories for annotation cluster 
## figures 
OUT_DIR_figs <- paste(output_dir,"/step7/figs7",sep='') 
OUT_dir_figs_annotate <- paste(OUT_DIR_figs,"/annotate/",sep='') 
dir.create(OUT_dir_figs_annotate)

## set output directory
#PWD=OUT_dir_figs_annotate
#setwd(PWD)

## add cluster annotation
cluster.ids<-par_annotate_labels

## set cluster resolution and rename cluster identities
seu_int <- SetIdent(seu_int, value = par_annotate_resolution)
names(cluster.ids) <- levels(seu_int)    
seu_int <- RenameIdents(seu_int, cluster.ids) 
seu_int$temp_temp_1 <- Idents(seu_int)

name_meta <- names(seu_int@meta.data) 
length <- length(name_meta)
name_meta[length] <- par_name_metadata
names(seu_int@meta.data) <- name_meta

##save RDS object
saveRDS(seu_int,paste(output_dir,'/step7/objs7','/seu_step7.rds', sep=''))

## print UMAP with cluster annotation
DimPlot(seu_int, reduction = "umap", label = TRUE, pt.size = 0.5, raster = FALSE) + NoLegend()
ggsave(file = paste(OUT_dir_figs_annotate,par_name_metadata,'_cluster_annotation.pdf', sep=''))

## print UMAP splitted
DimPlot(seu_int, reduction = "umap", split.by = "orig.ident", label = TRUE, pt.size = 0.5, raster = FALSE) + NoLegend()
ggsave(file = paste(OUT_dir_figs_annotate,par_name_metadata,'_split_cluster_annotation.pdf', sep=''), dpi = 300, height = 5, width = 10, unit = 'in')

## save metadata information
write.csv(colnames(seu_int[[]]), file= paste(output_dir,'/step7/info7/meta_info_seu_step7',".txt", sep=""))

## save rna expression matrix
if (tolower(par_save_RNA)=='yes') {
    mat <- GetAssayData(object = seu_int, assay = "RNA", slot = "data")
    writeMM(mat,file= paste(output_dir,'/step7/info7/seu',"_RNA.txt", sep=""))
}

## save metadata dataframe
if (tolower(par_save_metadata)=='yes') {
    write.csv(seu_int[[]], file = paste(output_dir,'/step7/info7/seu_MetaData.txt', sep=""), quote = TRUE, sep = ",")
}

## write session info
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step7/info7/sessionInfo_annotate.txt', sep=""))

if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}
}
