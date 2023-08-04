#!/usr/bin/env Rscript

####################
# step7
# paramters: par_futureglobalsmaxSize, level_cluster

args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr','xlsx','enrichR')
lapply(packages, library, character.only = TRUE)

sample_name<-list.files(path = paste(output_dir, "/step6/objs6",sep=""),pattern = "*.rds")

source(paste(output_dir,'/job_info/parameters/step7_par.txt',sep=""))
PWD=paste(output_dir,'/step7/info7/', sep='')
setwd(PWD)

options(future.globals.maxSize = par_futureglobalsmaxSize)

PSUE=paste(output_dir,'/step6/objs6/',sample_name, sep='')
ClusterMarkers<-read.csv(file = paste(output_dir,'/step7/objs7/', "ClusterMarkers.csv",sep=''),row.names = 1)
# scrnaboxR::annotation(level_cluster,ClusterMarkers,PWD,PSUE,top_sel,db)



setwd(PWD)
seu.q6 <- readRDS(PSUE)

if(dir.exists(paste(output_dir,'/step7/info7/annot_enrich',sep=''))){
    unlink(paste(output_dir,'/step7/info7/annot_enrich',sep=''), recursive=TRUE)
}

dir.create("annot_enrich")
for (i in sort(unlist(unique(seu.q6[[par_level_cluster]])))) {
  dir.create(paste0(PWD,"/annot_enrich", "/clust",i))
}
Idents(seu.q6) <- par_level_cluster
setEnrichrSite("Enrichr") # Human genes
for (i in sort(unlist(unique(seu.q6[[par_level_cluster]])))) {
  setwd(PWD)
  setwd(paste0('./annot_enrich/clust',i))
  N1.c0 <- ClusterMarkers %>% filter(cluster == i & avg_log2FC > 0)
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


writeLines(capture.output(sessionInfo()), paste(output_dir,'/step7/info7/sessionInfo_enrich.txt', sep=""))
if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}
