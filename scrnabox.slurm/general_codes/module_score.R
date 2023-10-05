####################
# Module score
####################

## load library
library(Seurat)
library(ggplot2)
library(vctrs)

## set cluster level
level_cluster=par_level_cluster #users must define their desired resolution in the parameters

## create output dir
OUT_DIR <- paste(output_dir,"/step7",sep='') 
dir.create(OUT_DIR)

## set default assay
DefaultAssay(seu_int) <- "RNA"

## set idents
Idents(seu_int) <- level_cluster

## read csv file which is manually prepared by the user. I have attached an example of this file in my message
gene_sets <- read.delim(par_module_score, header = T, sep = ",", na.strings=c("","NA")) 

#convert user inputed dataframe to named list
gene_lists <- list()                   # Create empty list
for(i in 1:ncol(gene_sets)) {             # Using for-loop to add columns to list
  gene_lists[[i]] <- gene_sets[ , i]
}
names(gene_lists) <- colnames(gene_sets)#set names of list

## remove NA from list
gene_lists <- lapply(gene_lists, function(x) x[!is.na(x)])

## calculate module score function
calculate_module_scores <- function(seu_int, gene_lists, OUT_DIR) {
  result <- list()
  for (i in names(gene_lists)) {
    # Add module scores for each gene list
    try( 
      seu_int <- AddModuleScore(
        seu_int,
        features = gene_lists[i],
        pool = NULL,
        nbin = 24,
        ctrl = 4,
        k = FALSE,
        name=paste0("moduleset_",i, "_")
      ), silent = TRUE)
    # set metadata into a dataframe
    df <- data.frame(seu_int@meta.data)
    names(df)
    df1 <- df[ , paste0("moduleset_",i, "_1") ]
    df1
    df1 <- data.frame(df1)
    rownames(df1) <- colnames(x = seu_int)
    # add to Seurat metadata
    seu_int <- AddMetaData(seu_int, df1, col.name = paste0("module_scores_",i))
    
    # print feature plot
    try( 
      FeaturePlot(seu_int,
                  features = paste0("module_scores_",i), label = TRUE, repel = TRUE, cols = c("white", "blue")), silent = TRUE)
    ggsave(file = paste(OUT_DIR, "/figs7/module_score/module_score_",i, ".pdf", sep='')) #change this step 7 output directory
    
    # add cluster annotation to cell-wise expression
    df1$cluster <- seu_int[[par_level_cluster]]
    colnames(df1)
    
    # calculate cluster-wise average module score 
    mean_scores_cluster <- df1 %>%
      group_by(cluster) %>%
      dplyr::summarize(Mean = mean(df1, na.rm=TRUE))
    
    result[[i]] <- mean_scores_cluster$Mean
    names(result[[i]]) <- mean_scores_cluster$cluster
  }
  # convert list object to dataframe
  df <- data.frame(result)
  df$clusters <- as.numeric(rownames(df)) -1
  # write to csv file
  write.csv(df, paste(OUT_DIR,"/info7/module_score/geneset_by_cluster.csv", sep=""))
  
}

## perform module score computation
calculate_module_scores(seu_int, gene_lists, OUT_DIR)
