#!/usr/bin/env Rscript

####################
# step7 -- find marker
####################

## set sample ID metadata column -- this is standard and does not require parameter modification
par_level_genotype <- "Sample_ID"

## load parameters
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
pipeline_home=args[3]

## load library
.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr', 'xlsx')
lapply(packages, library, character.only = TRUE)

## load existing Seurat objects
sample_name<-list.files(path = paste(output_dir, "/step6/objs6",sep=""),pattern = "*.rds")
seu_int<-readRDS(paste(output_dir,'/step6/objs6/',sample_name, sep=''))

## load parameters
source(paste(output_dir,'/job_info/parameters/step7_par.txt',sep=""))

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
