#!/usr/bin/env Rscript

####################
# step8 -- DGE
####################

## load parameters
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]

## load library
.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr','stringi','tidyverse','Matrix', 'ggrepel','DESeq2','EnhancedVolcano')
lapply(packages, library, character.only = TRUE)

## load parameters
source(paste(output_dir,'/job_info/parameters/step8_par.txt',sep=""))

############################################################################
## Wilcoxon all cells
############################################################################

if (tolower(par_run_wilcoxon_all_cells)=='yes') {

################### import the right Seurat object ###################
## load name of existing Seurat objects
sample_name<-list.files(path = paste(output_dir, "/step7/objs7",sep=""),pattern = "*.rds")

if(file.exists(paste(output_dir,'/step8/objs8/','seu_step8.rds', sep = ""))){
    seu_int<-readRDS(paste(output_dir,'/step8/objs8/','seu_step8.rds', sep=''))
}else{
    seu_int<-readRDS(paste(output_dir,'/step7/objs7/',sample_name, sep=''))
}
################### ############################## ###################

## set default assay to RNA
DefaultAssay(seu_int) <- "RNA"

### sample-sample dge wilcoxon analysis ####
dd<-read.csv(paste(output_dir,'/job_info/parameters/step8_contrast_wilcoxon_all_cells.txt',sep='/'), sep="")


### create directories
## base directory
OUT_dir_sample <- paste(output_dir,"/step8/wilcoxon_all_cells",sep='') 
dir.create(OUT_dir_sample)

## contrast-specific
for(i in 1:nrow(dd)){  
dir.create(paste(OUT_dir_sample,"/",dd[i,1],sep=''))
}

for(i in 1:nrow(dd)){  
Idents(seu_int) <- dd[i,2]    
DGE <- FindMarkers(seu_int, ident.1 = dd[i,3], ident.2 = dd[i,4],  logfc.threshold = 0)
#write dge
write.csv(DGE, file = paste(OUT_dir_sample,"/",dd[i,1],"/", dd[i,1],'_DEG.csv', sep=""), quote = FALSE, sep = ",")

# volcano plot
# volcano plot
## print volcano plot
    EnhancedVolcano(DGE,
    lab = rownames(DGE),
    x = 'avg_log2FC',
    y = 'p_val',
    #ylim = c(0,20),
    #xlim = c(-2,2),
    pCutoff = 0.001,
    FCcutoff = 1,
    pointSize = 2.0,
    labSize = 5, 
    #legendLabSize = 20,
    #subtitleLabSize = 20,
    #legendIconSize = 10,
    )
ggsave(file = paste(OUT_dir_sample,"/",dd[i,1],"/", dd[i,1],'_volcano_plot.pdf', sep=""))

}
###################

## print RDS object
saveRDS(seu_int, paste(output_dir,'/step8/objs8',"/seu_step8.rds", sep=""))

## save rna expression matrix
if (tolower(par_save_RNA)=='yes') {
    mat <- GetAssayData(object = seu_int, assay = "RNA", slot = "data")
    #  write.csv(mat, paste(output_dir,'/step2/info2/',sample_name[i],"_RNA.csv", sep=""))
    writeMM(mat,file= paste(output_dir,'/step8/info8/seu',"_RNA.txt", sep=""))
}

## save metadata dataframe
if (tolower(par_save_metadata)=='yes') {
    write.csv(seu_int[[]], file = paste(output_dir,'/step8/info8/seu_MetaData.txt', sep=""), quote = TRUE, sep = ",")
}

## write session info
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step8/info8/sessionInfo_wilcoxon_all_cells.txt', sep=""))
if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}
}

############################################################################
## Wilcoxon cell type groups
############################################################################

if (tolower(par_run_wilcoxon_celltype_groups)=='yes') {

################### import the right Seurat object ###################
## load name of existing Seurat objects
sample_name<-list.files(path = paste(output_dir, "/step7/objs7",sep=""),pattern = "*.rds")

if(file.exists(paste(output_dir,'/step8/objs8/','seu_step8.rds', sep = ""))){
    seu_int<-readRDS(paste(output_dir,'/step8/objs8/','seu_step8.rds', sep=''))
}else{
    seu_int<-readRDS(paste(output_dir,'/step7/objs7/',sample_name, sep=''))
}
################### ############################## ###################


## set default assay to RNA
DefaultAssay(seu_int) <- "RNA"

### sample-cell dge analysis ###
dd<-read.csv(paste(output_dir,'/job_info/parameters/step8_contrast_wilcoxon_celltype_groups.txt',sep='/'), sep="")

### create directories
## base directory
OUT_dir_sample <- paste(output_dir,"/step8/wilcoxon_celltype_groups",sep='') 
dir.create(OUT_dir_sample)

## contrast-specific
for(i in 1:nrow(dd)){  
dir.create(paste(OUT_dir_sample,"/",dd[i,1],sep=''))
}

for(i in 1:nrow(dd)){
Idents(seu_int) <- dd[i,2]    
celltype.sub.seu <- subset(seu_int, idents = dd[i,3])
Idents(celltype.sub.seu) <- dd[i,4]    
DGE <- FindMarkers(celltype.sub.seu, ident.1 = dd[i,5], ident.2 = dd[i,6],  logfc.threshold = 0)
#write dge
write.csv(DGE, file = paste(OUT_dir_sample,"/",dd[i,1],"/", dd[i,1],'_DEG.csv', sep=""), quote = FALSE, sep = ",")

# volcano plot
## print volcano plot
    EnhancedVolcano(DGE,
    lab = rownames(DGE),
    x = 'avg_log2FC',
    y = 'p_val',
    #ylim = c(0,20),
    #xlim = c(-2,2),
    pCutoff = 0.001,
    FCcutoff = 1,
    pointSize = 2.0,
    labSize = 5, 
    #legendLabSize = 20,
    #subtitleLabSize = 20,
    #legendIconSize = 10,
    )
ggsave(file = paste(OUT_dir_sample,"/",dd[i,1],"/", dd[i,1],'_volcano_plot.pdf', sep=""))


}
################### ########################## ###################

## print RDS object
saveRDS(seu_int, paste(output_dir,'/step8/objs8',"/seu_step8.rds", sep=""))

## save rna expression matrix
if (tolower(par_save_RNA)=='yes') {
    mat <- GetAssayData(object = seu_int, assay = "RNA", slot = "data")
    #  write.csv(mat, paste(output_dir,'/step2/info2/',sample_name[i],"_RNA.csv", sep=""))
    writeMM(mat,file= paste(output_dir,'/step8/info8/seu',"_RNA.txt", sep=""))
}

## save metadata dataframe
if (tolower(par_save_metadata)=='yes') {
    write.csv(seu_int[[]], file = paste(output_dir,'/step8/info8/seu_MetaData.txt', sep=""), quote = TRUE, sep = ",")
}

## write session info
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step8/info8/sessionInfo_wilcoxon_celltype_groups.txt', sep=""))
if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}
}



############################################################################
## pseudo bulk cell type groups
############################################################################

if (tolower(par_run_pseudo_bulk_celltype_groups)=='yes') {

################### import the right Seurat object ###################
## load name of existing Seurat objects
sample_name<-list.files(path = paste(output_dir, "/step7/objs7",sep=""),pattern = "*.rds")

if(file.exists(paste(output_dir,'/step8/objs8/','seu_step8.rds', sep = ""))){
    seu<-readRDS(paste(output_dir,'/step8/objs8/','seu_step8.rds', sep=''))
}else{
    seu<-readRDS(paste(output_dir,'/step7/objs7/',sample_name, sep=''))
}
################### ############################## ###################

## read in the pseudobulk contrast matrix
dd<-read.csv(paste(output_dir,'/job_info/parameters/step8_contrast_pseudo_bulk_celltype_groups.txt',sep='/'), sep="")
## store the contrast name
cont_name <- dd[1,1]
## create a psudobulk folder
sudo_dir <- paste(output_dir,'/step8/pseudo_bulk_celltype_groups', sep = "")
dir.create(sudo_dir)
outputpath <- sudo_dir
dir.create(paste(sudo_dir,"/",cont_name, sep = ""))
dir.create(paste(sudo_dir,"/",cont_name,"/figs", sep = ""))
dir.create(paste(sudo_dir,"/",cont_name,"/info", sep = ""))


## see how many variables the user inputted, not including contrast name
n_vars <- ncol(dd)
vars <- dd[1,c(2:n_vars)]

#Get the aggregated expression and calculate differential gene expression
Idents(seu) <- dd$CellType #cell types annotation
sum_counts <- AggregateExpression(seu, assay = "RNA", group.by = dput(as.character(vars)))
# this creates a list
sum_counts_df <- as.data.frame(sum_counts$RNA)
# print Agregated expression summary
write.csv(sum_counts_df, paste(sudo_dir,"/",cont_name,"/Aggregated_expression_summary.csv", sep= ""))

# Prep contrast matrix
df <- sum_counts_df
sample <- as.data.frame(colnames(df))
colnames(sample) <- "sample"
sample_names <- sample$sample
sample_info <- strsplit(sample_names, "_")
sample_info <- do.call(rbind, sample_info)

## create an empty list
list.results <- vector("list", length(vars))
vars <- colnames(dd[,c(2:n_vars)])
names(list.results) <- dput(as.character(vars)) ###

# Create columns for all variables
for (i in 1:length(list.results)){
list.results[[i]] <- sample_info[, i]
}


test <- data.frame(list.results)
test2 <- cbind(sample, test)
sample <- test2
meta.df <- sample
df.meta <- meta.df

Celltypes <- unique(df.meta$CellType)
df.trans <- as.data.frame(t(df))
# Create an empty list to store the results for each cell type
list.results <- vector("list", length(Celltypes))

# Define all_contrasts based on unique levels in the MainContrast column
all_contrasts <- unique(df.meta$MainContrast)                    

# Loop through each cell type and perform DESeq analysis
for (i in seq_along(Celltypes)) {
  # Subset the expression dataframe by the current cell type
  print(Celltypes[i])
  df_sub <- df.trans[grepl(paste0("^", Celltypes[i], "_"), rownames(df.trans)), ]
  print(dim(df_sub))
  # test one cell group
  #i = "NPC-div"
  df.meta_sub <- df.meta[df.meta$CellType == Celltypes[i], ]
  print(dim(df.meta_sub))
  # Prepare the DESeq object
  dft <- as.data.frame(t(df_sub)) # Transpose the subset dataframe to get genes as rows and samples as columns
  dfi <- lapply(dft, as.integer)
  dfi <- as.data.frame(dfi)
  rownames(dfi) <- rownames(dft)
  
  # Check if the MainContrast has at least two unique values
  if (length(unique(df.meta_sub$MainContrast)) >= 2) {
    # Create the DESeqDataSet object using the subset dataframe and metadata
    dds <- DESeqDataSetFromMatrix(countData = dfi, colData = df.meta_sub, design = ~MainContrast)
    # Perform DESeq analysis
    dds <- DESeq(dds)
    # Store the DESeq results in the list with the cell type as the list index
    res <- results(dds)
    list2 <- list()
    list2[["dds"]] <- dds
    list2[["results"]] <- res
    
    # Initialize an empty list to store the results for each contrast
    if (length(all_contrasts) > 1) {
      all_results <- list()
      # Loop through each contrast and calculate the results
      for (j in 1:(length(all_contrasts) - 1)) {
        for (k in (j + 1):length(all_contrasts)) {
          contrast_level1 <- all_contrasts[j]
          contrast_level2 <- all_contrasts[k]

          # Check if both levels have at least one sample
          if (any(dds$MainContrast %in% c(contrast_level1, contrast_level2))) {
            # Check if contrast levels are different
            if (contrast_level1 != contrast_level2) {
              contrast_name <- paste("MainContrast", contrast_level1, "vs", contrast_level2)
              # Filter rows with complete cases for the current contrast levels
              complete_cases <- complete.cases(dds$MainContrast, dds$MainContrast %in% c(contrast_level1, contrast_level2))
              # Subset the DESeq object and calculate the contrast results
              dds_sub <- dds[complete_cases, ]

              # Check if both levels still exist in the subset after removing missing values
              if (contrast_level1 %in% unique(dds_sub$MainContrast) && contrast_level2 %in% unique(dds_sub$MainContrast)) {
                contrast_result <- results(dds_sub, contrast = c("MainContrast", contrast_level1, contrast_level2), name = contrast_name)
                all_results[[contrast_name]] <- contrast_result
              } else {
                message(paste("Skipping", contrast_name, "due to missing contrast levels in MainContrast."))
              }
            } else {
              message(paste("Skipping", contrast_level1, "vs", contrast_level2, "since they are the same level."))
            }
          }
        }
      }
      
      # Add the contrast results to the list
      list2[["contrast_results"]] <- all_results
    }
  
    # Store the list for the current cell type in the appropriate slot
    list.results[[i]] <- list2
  } else {
    message(paste("Skipping", Celltypes[i], "due to insufficient unique values in MainContrast."))
  }
}


names(list.results) <- Celltypes
length(list.results) # there are ten features, one for each cell type

# make contrast name
dir.create(paste(sudo_dir,"/",cont_name, sep = ""))

# make a saving loop 
celltypes.res <- names(list.results)
for (i in celltypes.res){
  for (j in length(names(list.results[[i]][["contrast_results"]]))){
    result <- as.data.frame(list.results[[i]][["contrast_results"]][[j]])
    celltype <- i
    contrast <- names(list.results[[i]][["contrast_results"]])
    write.csv(result, paste(sudo_dir,"/",cont_name,"/info","/DGE_",celltype,contrast,".csv", sep = ""))
    ## print volcano plot
    EnhancedVolcano(result,
    lab = rownames(result),
    x = 'log2FoldChange',
    y = 'pvalue',
    #ylim = c(0,20),
    #xlim = c(-2,2),
    pCutoff = 0.001,
    FCcutoff = 1,
    pointSize = 2.0,
    labSize = 5, 
    #legendLabSize = 20,
    #subtitleLabSize = 20,
    #legendIconSize = 10,
    ) 
    ggsave(file = paste(sudo_dir,"/",cont_name,"/figs","/DGE_",celltype,contrast,".pdf", sep = ""))
  }
}


#A summarize results function
summarize_contrast <- function(result, 
                               adjp_threshold = 0.05,
                               logfoldchange_thesh = 0) {
  num_de_genes <- sum(result$padj <= adjp_threshold, na.rm = TRUE)
  num_downregulated <- sum(result$padj <= adjp_threshold & result$log2FoldChange <= -logfoldchange_thesh, na.rm = TRUE)
  num_upregulated <- sum(result$padj <= adjp_threshold & result$log2FoldChange >= logfoldchange_thesh, na.rm = TRUE)
  
  return(data.frame(NumGenes = num_de_genes, NumDownregulated = num_downregulated, NumUpregulated = num_upregulated))
}


# Initialize an empty list to store the summary results for each cell type and contrast
summary_results_list <- vector("list", length(Celltypes))
names(summary_results_list) <- Celltypes

# Loop through each cell type and summarize the results for each contrast
for (i in seq_along(Celltypes)) {
  # Get the results for the current cell type
  all_results_celltype <- list.results[[Celltypes[i]]][["contrast_results"]]
  
  # Initialize lists to store the summary results for each contrast
  contrasts <- names(all_results_celltype)
  num_contrasts <- length(contrasts)
  celltype_summary <- data.frame(
    Celltype = rep(Celltypes[i], num_contrasts),
    Contrast = contrasts,
    DGE_total = numeric(num_contrasts),
    DGE_up = numeric(num_contrasts),
    DGE_down = numeric(num_contrasts)
  )
  
  # Loop through each contrast and summarize the results
  for (j in seq_along(contrasts)) {
    contrast_name <- contrasts[j]
    contrast_result <- all_results_celltype[[contrast_name]]
    summary_result <- summarize_contrast(contrast_result,
                                         adjp_threshold = 0.1,
                                         logfoldchange_thesh = 0)
    
    # Store the summary results for the current contrast
    celltype_summary[j, "DGE_total"] <- summary_result$NumGenes
    celltype_summary[j, "DGE_up"] <- summary_result$NumUpregulated
    celltype_summary[j, "DGE_down"] <- summary_result$NumDownregulated
  }
  
  # Store the summary results for the current cell type
  summary_results_list[[Celltypes[i]]] <- celltype_summary
}

# Combine all the summary results into a single data frame
summary_table <- do.call(rbind, summary_results_list)
write.csv(summary_table,paste(sudo_dir,"/",cont_name,"/PseudoBulk_DGEsummarytable.csv", sep = ""))
# Print the summary table
print(summary_table)

## write session info
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step8/info8/sessionInfo_pseudo_bulk_celltype_groups.txt', sep=""))
if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}

}


############################################################################
## pseudo bulk all cells
############################################################################

if (tolower(par_run_pseudo_bulk_all_cells)=='yes') {

################### import the right Seurat object ###################
## load name of existing Seurat objects
sample_name<-list.files(path = paste(output_dir, "/step7/objs7",sep=""),pattern = "*.rds")

if(file.exists(paste(output_dir,'/step8/objs8/','seu_step8.rds', sep = ""))){
    seu<-readRDS(paste(output_dir,'/step8/objs8/','seu_step8.rds', sep=''))
}else{
    seu<-readRDS(paste(output_dir,'/step7/objs7/',sample_name, sep=''))
}
################### ############################## ###################

## read in the pseudobulk contrast matrix
dd<-read.csv(paste(output_dir,'/job_info/parameters/step8_contrast_pseudo_bulk_all_cells.txt',sep='/'), sep="")
## store the contrast name
cont_name <- dd[1,1]
## create a psudobulk folder
sudo_dir <- paste(output_dir,'/step8/pseudo_bulk_all_cells', sep = "")
dir.create(sudo_dir)
outputpath <- sudo_dir
dir.create(paste(sudo_dir,"/",cont_name, sep = ""))
#dir.create(paste(sudo_dir,"/",cont_name,"/figs", sep = ""))
#dir.create(paste(sudo_dir,"/",cont_name,"/info", sep = ""))


## see how many variables the user inputted, not including contrast name
n_vars <- ncol(dd)
vars <- dd[1,c(2:n_vars)]

#Get the aggregated expression and calculate differential gene expression
#Idents(seu) <- dd$CellType #cell types annotation
sum_counts <- AggregateExpression(seu, assay = "RNA", group.by = dput(as.character(vars)))
# this creates a list
sum_counts_df <- as.data.frame(sum_counts$RNA)
# print Agregated expression summary
write.csv(sum_counts_df, paste(sudo_dir,"/",cont_name,"/Aggregated_expression_summary.csv", sep= ""))

# Prep contrast matrix
df <- sum_counts_df
sample <- as.data.frame(colnames(df))
colnames(sample) <- "sample"
sample_names <- sample$sample
sample_info <- strsplit(sample_names, "_")
sample_info <- do.call(rbind, sample_info)

## create an empty list
list.results <- vector("list", length(vars))
vars <- colnames(dd[,c(2:n_vars)])
names(list.results) <- dput(as.character(vars)) ###

# Create columns for all variables
for (i in 1:length(list.results)){
list.results[[i]] <- sample_info[, i]
}


test <- data.frame(list.results)
test2 <- cbind(sample, test)
sample <- test2
meta.df <- sample
df.meta <- meta.df

df.meta$CellType <- "AllCells"

Celltypes <- unique(df.meta$CellType)
df.trans <- as.data.frame(t(df))
# Create an empty list to store the results for each cell type
list.results <- vector("list", length(Celltypes))

# Define all_contrasts based on unique levels in the MainContrast column
all_contrasts <- unique(df.meta$MainContrast)                    


# Loop through each cell type and perform DESeq analysis
for (i in seq_along(Celltypes)) {
  # Subset the expression dataframe by the current cell type
  print(Celltypes[i])
  df_sub <- df.trans[rownames(df.trans), ]
  print(dim(df_sub))
  # test one cell group
  #i = "NPC-div"
  df.meta_sub <- df.meta[df.meta$CellType == Celltypes[i], ]
  print(dim(df.meta_sub))
  # Prepare the DESeq object
  dft <- as.data.frame(t(df_sub)) # Transpose the subset dataframe to get genes as rows and samples as columns
  dfi <- lapply(dft, as.integer)
  dfi <- as.data.frame(dfi)
  rownames(dfi) <- rownames(dft)
  
  # Check if the MainContrast has at least two unique values
  if (length(unique(df.meta_sub$MainContrast)) >= 2) {
    # Create the DESeqDataSet object using the subset dataframe and metadata
    dds <- DESeqDataSetFromMatrix(countData = dfi, colData = df.meta_sub, design = ~MainContrast)
    # Perform DESeq analysis
    dds <- DESeq(dds)
    # Store the DESeq results in the list with the cell type as the list index
    res <- results(dds)
    list2 <- list()
    list2[["dds"]] <- dds
    list2[["results"]] <- res
    
    # Initialize an empty list to store the results for each contrast
    if (length(all_contrasts) > 1) {
      all_results <- list()
      # Loop through each contrast and calculate the results
      for (j in 1:(length(all_contrasts) - 1)) {
        for (k in (j + 1):length(all_contrasts)) {
          contrast_level1 <- all_contrasts[j]
          contrast_level2 <- all_contrasts[k]

          # Check if both levels have at least one sample
          if (any(dds$MainContrast %in% c(contrast_level1, contrast_level2))) {
            # Check if contrast levels are different
            if (contrast_level1 != contrast_level2) {
              contrast_name <- paste("MainContrast", contrast_level1, "vs", contrast_level2)
              # Filter rows with complete cases for the current contrast levels
              complete_cases <- complete.cases(dds$MainContrast, dds$MainContrast %in% c(contrast_level1, contrast_level2))
              # Subset the DESeq object and calculate the contrast results
              dds_sub <- dds[complete_cases, ]

              # Check if both levels still exist in the subset after removing missing values
              if (contrast_level1 %in% unique(dds_sub$MainContrast) && contrast_level2 %in% unique(dds_sub$MainContrast)) {
                contrast_result <- results(dds_sub, contrast = c("MainContrast", contrast_level1, contrast_level2), name = contrast_name)
                all_results[[contrast_name]] <- contrast_result
              } else {
                message(paste("Skipping", contrast_name, "due to missing contrast levels in MainContrast."))
              }
            } else {
              message(paste("Skipping", contrast_level1, "vs", contrast_level2, "since they are the same level."))
            }
          }
        }
      }
      
      # Add the contrast results to the list
      list2[["contrast_results"]] <- all_results
    }
  
    # Store the list for the current cell type in the appropriate slot
    list.results[[i]] <- list2
  } else {
    message(paste("Skipping", Celltypes[i], "due to insufficient unique values in MainContrast."))
  }
}


names(list.results) <- Celltypes
length(list.results) # there are ten features, one for each cell type

# make a saving loop 
celltypes.res <- names(list.results)
for (i in celltypes.res){
  for (j in length(names(list.results[[i]][["contrast_results"]]))){
    result <- as.data.frame(list.results[[i]][["contrast_results"]][[j]])
    celltype <- i
    contrast <- names(list.results[[i]][["contrast_results"]])
    write.csv(result, paste(sudo_dir,"/",cont_name,"/DGE_",celltype,contrast,".csv", sep = ""))
    ## print volcano plot
    EnhancedVolcano(result,
    lab = rownames(result),
    x = 'log2FoldChange',
    y = 'pvalue',
    #ylim = c(0,20),
    #xlim = c(-2,2),
    pCutoff = 0.001,
    FCcutoff = 1,
    pointSize = 2.0,
    labSize = 5, 
    #legendLabSize = 20,
    #subtitleLabSize = 20,
    #legendIconSize = 10,
    ) 
    ggsave(file = paste(sudo_dir,"/",cont_name,"/DGE_",celltype,contrast,".pdf", sep = ""))
  }
}


#A summarize results function
summarize_contrast <- function(result, 
                               adjp_threshold = 0.05,
                               logfoldchange_thesh = 0) {
  num_de_genes <- sum(result$padj <= adjp_threshold, na.rm = TRUE)
  num_downregulated <- sum(result$padj <= adjp_threshold & result$log2FoldChange <= -logfoldchange_thesh, na.rm = TRUE)
  num_upregulated <- sum(result$padj <= adjp_threshold & result$log2FoldChange >= logfoldchange_thesh, na.rm = TRUE)
  
  return(data.frame(NumGenes = num_de_genes, NumDownregulated = num_downregulated, NumUpregulated = num_upregulated))
}


# Initialize an empty list to store the summary results for each cell type and contrast
summary_results_list <- vector("list", length(Celltypes))
names(summary_results_list) <- Celltypes

# Loop through each cell type and summarize the results for each contrast
for (i in seq_along(Celltypes)) {
  # Get the results for the current cell type
  all_results_celltype <- list.results[[Celltypes[i]]][["contrast_results"]]
  
  # Initialize lists to store the summary results for each contrast
  contrasts <- names(all_results_celltype)
  num_contrasts <- length(contrasts)
  celltype_summary <- data.frame(
    Celltype = rep(Celltypes[i], num_contrasts),
    Contrast = contrasts,
    DGE_total = numeric(num_contrasts),
    DGE_up = numeric(num_contrasts),
    DGE_down = numeric(num_contrasts)
  )
  
  # Loop through each contrast and summarize the results
  for (j in seq_along(contrasts)) {
    contrast_name <- contrasts[j]
    contrast_result <- all_results_celltype[[contrast_name]]
    summary_result <- summarize_contrast(contrast_result,
                                         adjp_threshold = 0.1,
                                         logfoldchange_thesh = 0)
    
    # Store the summary results for the current contrast
    celltype_summary[j, "DGE_total"] <- summary_result$NumGenes
    celltype_summary[j, "DGE_up"] <- summary_result$NumUpregulated
    celltype_summary[j, "DGE_down"] <- summary_result$NumDownregulated
  }
  
  # Store the summary results for the current cell type
  summary_results_list[[Celltypes[i]]] <- celltype_summary
}

# Combine all the summary results into a single data frame
summary_table <- do.call(rbind, summary_results_list)
write.csv(summary_table,paste(sudo_dir,"/",cont_name,"/PseudoBulk_DGEsummarytable.csv", sep = ""))
# Print the summary table
print(summary_table)

## write session info
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step8/info8/sessionInfo_pseudo_bulk_celltype_groups.txt', sep=""))
if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}

}