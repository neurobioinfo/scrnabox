#!/usr/bin/env Rscript

####################
# step7 -- cluster marker GSEA
####################

## load parameters
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
.libPaths(r_lib_path)

## load libraries
packages<-c('Seurat','ggplot2', 'dplyr','xlsx','enrichR')
lapply(packages, library, character.only = TRUE)

## load in parameters files
source(paste(output_dir,'/job_info/parameters/step7_par.txt',sep=""))

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

