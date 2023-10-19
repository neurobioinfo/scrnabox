#!/usr/bin/env Rscript

###############################################################################
# step7 -- cluster annotation
###############################################################################

## set sample ID metadata column -- this is standard and does not require parameter modification
#par_level_genotype <- "Sample_ID"

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
