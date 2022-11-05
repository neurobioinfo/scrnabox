# Generate the annotation
annotation<-function(level_cluster,PWD,PSUE,top_sel,db) {
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(xlsx)
  library(enrichR)
  setwd(PWD)
  seu.q6 <- readRDS(PSUE)
  # mkdir -p annot_enrich
  dir.create("annot_enrich")
  for (i in sort(unlist(unique(seu.q6[[level_cluster]])))) {
    dir.create(paste0(PWD,"/annot_enrich", "/clust",i))
  }
  Idents(seu.q6) <- level_cluster
  ClusterMarkers <- FindAllMarkers(seu.q6, only.pos = TRUE)
  top5 <- ClusterMarkers %>% group_by(cluster) %>% top_n(n=top_sel, wt = avg_log2FC)
  write.csv(top5,"top_sel.csv")
  for (i in  sort(unlist(unique(seu.q6[[level_cluster]])))) {
    N1.c0 <- ClusterMarkers %>% filter(cluster == i & avg_log2FC > 0)
    genes <- N1.c0$gene
    write.xlsx(genes, file="cluster_just_genes.xlsx",sheetName=paste0('cluster',i), row.names=FALSE,append=TRUE)
    write.xlsx(N1.c0, file="cluster_whole.xlsx",sheetName=paste0('cluster',i), row.names=FALSE,append=TRUE)
  }
  heat_map<-DoHeatmap(seu.q6, features = top5$gene, size=3, angle =90, group.bar.height = 0.02, group.by = level_cluster)
  ggsave(file = "heatmap.pdf")
  dim_plot<-DimPlot(seu.q6,group.by = level_cluster)
  ggsave(file = "umap.pdf")
  setEnrichrSite("Enrichr") # Human genes
  for (i in sort(unlist(unique(seu.q6[[level_cluster]])))) {
    setwd(PWD)
    setwd(paste0('./annot_enrich/clust',i))
    N1.c0 <- ClusterMarkers %>% filter(cluster == i & avg_log2FC > 0)
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
