#!/usr/bin/env Rscript

####################
# step7 -- marker
####################
# Parameters: par_top_sel, par_level_cluster, par_level_genotype

## load parameters
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
pipeline_home=args[3]

## load library
.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr', 'xlsx')
lapply(packages, library, character.only = TRUE)

## load list of existing Seurat objects
sample_name<-list.files(path = paste(output_dir, "/step6/objs6",sep=""),pattern = "*.rds")
seu_int<-readRDS(paste(output_dir,'/step6/objs6/',sample_name, sep=''))

## load parameters
source(paste(output_dir,'/job_info/parameters/step7_par.txt',sep=""))

## set cell identity to the clustering resolution defined by the user
Idents(seu_int) <- par_level_cluster

## create directories for each cluster annotation method
## figures 
OUT_DIR_figs <- paste(output_dir,"/step7/figs7",sep='') 
OUT_dir_figs_marker <- paste(OUT_DIR_figs,"/marker/",sep='') 
dir.create(OUT_dir_figs_marker)
OUT_dir_figs_module_score <- paste(OUT_DIR_figs,"/module_score/",sep='') 
dir.create(OUT_dir_figs_module_score)
OUT_dir_figs_visualize_features <- paste(OUT_DIR_figs,"/visualize_select_features/",sep='') 
dir.create(OUT_dir_figs_visualize_features)
## info
OUT_DIR_info <- paste(output_dir,"/step7/info7",sep='') 
OUT_dir_info_marker <- paste(OUT_DIR_info,"/marker/",sep='') 
dir.create(OUT_dir_info_marker)
OUT_dir_info_module_score <- paste(OUT_DIR_info,"/module_score/",sep='') 
dir.create(OUT_dir_info_module_score)

## set output directory
PWD=OUT_dir_info_marker
setwd(PWD)

## identify cluster marker
ClusterMarkers <- FindAllMarkers(seu_int, only.pos = TRUE)

## identify top n marker genes for each cluster -- the number of top marker genes for each cluster is defined by the user
top5 <- ClusterMarkers %>% group_by(cluster) %>% top_n(n=par_top_sel, wt = avg_log2FC)
write.csv(top5,"top_sel.csv")

## print cluster information
if(file.exists(paste(output_dir,'/step7/info7/', "cluster_whole.xlsx",sep=''))){
    file.remove(paste(output_dir,'/step7/info7/', "cluster_whole.xlsx",sep=''))
}

if(file.exists(paste(output_dir,'/step7/info7/', "cluster_just_genes.xlsx",sep=''))){
    file.remove(paste(output_dir,'/step7/info7/', "cluster_just_genes.xlsx",sep=''))
}
for (i in  sort(unlist(unique(seu_int[[par_level_cluster]])))) {
  N1.c0 <- ClusterMarkers %>% filter(cluster == i & avg_log2FC > 0)
  genes <- N1.c0$gene
  write.xlsx(genes, file="cluster_just_genes.xlsx",sheetName=paste0('cluster',i), row.names=FALSE,append=TRUE)
  write.xlsx(N1.c0, file="cluster_whole.xlsx",sheetName=paste0('cluster',i), row.names=FALSE,append=TRUE)
}

## print heatmap for top markers genes
heat_map<-DoHeatmap(seu_int, features = top5$gene, size=3, angle =90, group.bar.height = 0.02, group.by = par_level_cluster)
ggsave(file = paste(OUT_dir_figs_marker,'heatmap.pdf', sep=''))

## print UMAP split by cluster
dim_plot<-DimPlot(seu_int,group.by = par_level_cluster)
ggsave(file = paste(output_dir,'/step7/figs7','/umap.pdf', sep=''))

## print UMAP split by sample identity
dim_plot<-DimPlot(seu_int,split.by = par_level_genotype, ncol = 3)
ggsave(file = paste(output_dir,'/step7/figs7','/umap_splitted.pdf', sep=''))

## save cluster marker information
saveRDS(ClusterMarkers, paste(OUT_dir_info_marker,"ClusterMarkers.rds", sep=""))
write.csv(ClusterMarkers,paste(OUT_dir_info_marker,"ClusterMarkers.csv",sep=''))

## compute module score
if (tolower(par_compute_module_score)=='yes') {
source(paste(pipeline_home,'/general_codes/module_score.R',sep=''))
}

## visulize select features
if (tolower(par_visualize_select_features)=='yes') {
select_features <- par_select_features
## violin plot
vln_plt <- VlnPlot(seu_int, features = select_features, group.by = par_level_cluster )
ggsave(file = paste(OUT_dir_figs_visualize_features,'select_feature_violin_plot.pdf', sep=''))
## feature plot
Feature_plt <- FeaturePlot(seu_int, features = select_features)
ggsave(file = paste(OUT_dir_figs_visualize_features ,'select_feature_feature_plot.pdf', sep=''))
#dotplot
dot_plt <- DotPlot(seu_int, features = select_features, group.by = par_level_cluster)
ggsave(file = paste(OUT_dir_figs_visualize_features ,'select_feature_dot_plot.pdf', sep=''))
}

## save session info
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step7/info7/sessionInfo_marker.txt', sep=""))

if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}


