#!/usr/bin/env Rscript

####################
# step7amarker
## Marker

args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr', 'xlsx')
lapply(packages, library, character.only = TRUE)

sample_name<-list.files(path = paste(output_dir, "/step6/objs6",sep=""),pattern = "*.rds")
seu_int<-readRDS(paste(output_dir,'/step6/objs6/',sample_name, sep=''))

source(paste(output_dir,'/job_info/parameters/step7_par.txt',sep=""))

PWD=paste(output_dir,'/step7/info7/', sep='')
setwd(PWD)

Idents(seu_int) <- level_cluster
ClusterMarkers <- FindAllMarkers(seu_int, only.pos = TRUE)
write.csv(ClusterMarkers,paste(output_dir,'/step7/objs7/', "ClusterMarkers.csv",sep=''))
top5 <- ClusterMarkers %>% group_by(cluster) %>% top_n(n=top_sel, wt = avg_log2FC)
write.csv(top5,"top_sel.csv")
for (i in  sort(unlist(unique(seu_int[[level_cluster]])))) {
  N1.c0 <- ClusterMarkers %>% filter(cluster == i & avg_log2FC > 0)
  genes <- N1.c0$gene
  write.xlsx(genes, file="cluster_just_genes.xlsx",sheetName=paste0('cluster',i), row.names=FALSE,append=TRUE)
  write.xlsx(N1.c0, file="cluster_whole.xlsx",sheetName=paste0('cluster',i), row.names=FALSE,append=TRUE)
}


heat_map<-DoHeatmap(seu_int, features = top5$gene, size=3, angle =90, group.bar.height = 0.02, group.by = level_cluster)
ggsave(file = paste(output_dir,'/step7/figs7','/heatmap.pdf', sep=''))

dim_plot<-DimPlot(seu_int,group.by = level_cluster)
ggsave(file = paste(output_dir,'/step7/figs7','/umap.pdf', sep=''))

dim_plot<-DimPlot(seu_int,split.by = level_genotype)
ggsave(file = paste(output_dir,'/step7/figs7','/umap_splitted.pdf', sep=''))

saveRDS(ClusterMarkers, paste(output_dir,'/step7/info7',"/ClusterMarkers.rds", sep=""))
write.xlsx(ClusterMarkers, file=paste(output_dir,'/step7/info7',"/cluster_just_genes.xlsx", sep=""))

writeLines(capture.output(sessionInfo()), paste(output_dir,'/step7/info7/sessionInfo_marker.txt', sep=""))
