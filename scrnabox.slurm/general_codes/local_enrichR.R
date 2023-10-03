#set up the environment
library(Seurat)
library(dplyr)
library(ggplot2)
library(enrichR)

############################################
#set up the parameters
############################################
PWD <-'/Users/mfiorini/Desktop/scRNA_pipeline/enrichr1.38_practuce'
cluster_marker <- '/Users/mfiorini/Desktop/scRNA_pipeline/ClusterMarkers.csv'
db <- c('Descartes_Cell_Types_and_Tissue_2021','CellMarker_Augmented_2021','Azimuth_Cell_Types_2021')

############################################
#prepare the function
############################################
#set up the function
annotation<-function(PWD,cluster_marker,db) {
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(enrichR)
  setwd(PWD)

  ## need to read in the cluster marker csv file csv file
  cluster_marker <- read.delim(cluster_marker, header = T, sep = ",") 
  
  
  #create directories based on the clusters present in the csv file
  dir.create("annot_enrich") 
  for (i in unique(cluster_marker$cluster)) {
    dir.create(paste0(PWD,"/annot_enrich", "/clust",i))
  }

  for (i in unique(cluster_marker$cluster)) {
    setwd(PWD)
    setwd(paste0('./annot_enrich/clust',i))
    N1.c0 <- cluster_marker %>% filter(cluster == i & avg_log2FC > 0)
    genes <- N1.c0$gene
    N1.c0.Er <- enrichr(genes, databases = db)
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
}


############################################
#Annotate
############################################
annotation(PWD,cluster_marker,db)

