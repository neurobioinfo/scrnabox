#add module score
#load packages
library(Seurat)
library(ggplot2)
library(vctrs)

#read in seurat object
# seu_int <- seu_int #we will take step 6 seurat object
# par_level_cluster='integrated_snn_res.1'
level_cluster=par_level_cluster #users must define their desired resolution in the parameters
OUT_DIR <- paste(output_dir,"/step6/module_score",sep='') # this will be step 7 output dir in the pipeline
dir.create(OUT_DIR)
Idents(seu_int) <- level_cluster

#read csv file which is manually prepared by the user. I have attached an example of this file in my message
gene_sets <- read.delim(par_module_score, header = T, sep = ",", na.strings=c("","NA")) #users must provide this, we can let them define the location in the parameters

#convert user inputed dataframe to named list
gene_lists <- list()                   # Create empty list
for(i in 1:ncol(gene_sets)) {             # Using for-loop to add columns to list
  gene_lists[[i]] <- gene_sets[ , i]
}
names(gene_lists) <- colnames(gene_sets)#set names of list

#remove NA from list
gene_lists <- lapply(gene_lists, function(x) x[!is.na(x)])

calculate_module_scores <- function(seu_int, gene_lists, OUT_DIR) {
  result <- list()
  for (i in names(gene_lists)) {
    # Add module scores for each gene list
    try( 
    seu_int <- AddModuleScore(
      seu_int,
      features = gene_lists[[i]],
      pool = NULL,
      nbin = 24,
      ctrl = 4,
      k = FALSE,
      name=paste0("moduleset_",i, "_")
    ), silent = TRUE)
    #set metadata into a dataframe
    df <- data.frame(seu_int@meta.data)
    #select the columns of interest
    df1 <- df[ , grepl("^moduleset_", names( df ), perl = TRUE ) ]
    
    #calculate mean expression across the gene set
    df1$mean <- rowMeans(df1)
    mean_scores<- df1$mean
    names(mean_scores) <- colnames(x = seu_int)
    mean_scores <- data.frame(mean_scores)
    rownames(mean_scores) <- colnames(x = seu_int)
    #add to metadat
    seu_int <- AddMetaData(seu_int, mean_scores, col.name = paste0("mean_scores_",i))
    
    #print feature plot
    try( 
    FeaturePlot(seu_int,
                features = paste0("mean_scores_",i), label = TRUE, repel = TRUE), silent = TRUE)
    ggsave(file = paste(OUT_DIR, "/module_score",i, ".png", sep='')) #change this step 7 output directory
    
    #add cluster annotation to cell-wise expression
    mean_scores$cluster <- seu_int[[par_level_cluster]]

    #calculate cluster-wise average expression 
    mean_scores_cluster <- mean_scores %>%
      group_by(cluster) %>%
      dplyr::summarize(Mean = mean(mean_scores, na.rm=TRUE))
    
    result[[i]] <- mean_scores_cluster$Mean
    names(result[[i]]) <- mean_scores_cluster$cluster
  }
  #convert list object to dataframe
  # result
  df <- data.frame(result)
  df$clusters <- rownames(df)
  # df
  #write to csv file
  write.csv(df, paste(OUT_DIR,"/geneset_by_cluster.csv", sep=""))
  # seu_int
  
}

calculate_module_scores(seu_int, gene_lists, OUT_DIR)


#parameters
#level_cluster
#gene set --> location of csv file (users can create this in excel and then save it as a csv file)

