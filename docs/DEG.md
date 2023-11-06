## DEG analysis
This tutorial provides the code used for producing Figure 5 of our [pre-print manuscript](); downstream analyses of the differential gene expression (DGE) results for the midbrain dataset ([Smajic et al. 2022](https://academic.oup.com/brain/article/145/3/964/6469020)). The tutorial can be broken down into three sections: <br />

 - [DEG summary](#1-deg-summary)<br />
 - [Cell-based DGE enrichment analysis](#2-cell-based-dge-enrichment-analysis)<br />
 - [Sample-based DGE enrichment analysis](#3-sample-based-dge-enrichment-analysis).<br />

The outputs from Step 8 of the [scRNAbox analuysis of the midbrain dataset](Dataset1.md) tutorial are used as input.
- - - -
## Set up
**Set seed**

    set.seed(1234)

**Load libraries**

    library(ggrepel)
    library(ggplot2)
    library(stringr)
    BiocManager::install("clusterProfiler")
    BiocManager::install("pathview")
    BiocManager::install("enrichplot")
    library(clusterProfiler)
    library(enrichplot)
    organism = "org.Hs.eg.db"
    BiocManager::install(organism, character.only = TRUE)
    library(organism, character.only = TRUE)
    require(DOSE)
    library(stringr)
    library(dplyr)

- - - -

## 1. DEG summary

### 1.1. load cell-based all cells DEG file

    ## load DEG csv file
    all <- read.delim("/Users/mfiorini/Desktop/scRNA_pipeline/Manuscript/Manuscript_Figures/Smajic3/Step5/wilcoxon_all_cells/HCvPD/HCvPD_DEG.csv", header = T, sep = ",")

### 1.2. load cell-based cell type groups DEG files

    ## list of cell types used for DEG contrasts
    cell_based_cell_types <- c('Astrocytes', 'DaN', 'Endothelial', 'Ependymal', 'Excitatory', 'GABA', 'Inhibitory', 'Microglia', 'Oligocendrocytes', 'OPC', 'Pericytes')

    ## create empty list object
    list <- list()

    ## load DEG csv files
    for(i in unique(cell_based_cell_types)){
      list[i] <- list(read.delim(paste0("/Users/mfiorini/Desktop/scRNA_pipeline/Manuscript/Manuscript_Figures/Smajic3/Step5/wilcoxon_celltype_groups/", i,"PDvHC/",i,"PDvHC_DEG.csv"), header = T, sep = ","))
    }

    ## create data frame for each cell type
    list2env(list,envir=.GlobalEnv)

### 1.3. load sample-based all cells DEG file

    ## load DEG csv file
    all_sample <- read.delim("/Users/mfiorini/Desktop/scRNA_pipeline/Manuscript/Manuscript_Figures/Smajic3/Step5/pseudo_bulk_all_cells/PDvControl/DGE_AllCellsMainContrast HC vs PD.csv", header = T, sep = ",")

### 1.4. load sample-based cell type groups DEG files

    ## list of cell types used for DEG contrasts
    sample_based_cell_types <- c('astro', 'DaN', 'endo', 'ependymal', 'excit', 'GABA', 'inhib', 'mg', 'Olig', 'opc', 'peri')

    ## create empty list object
    list <- list()

    ## load DEG csv files
    for(i in unique(sample_based_cell_types)){
      list[paste0(i,"_sample")] <- list(read.delim(paste0("/Users/mfiorini/Desktop/scRNA_pipeline/Manuscript/Manuscript_Figures/Smajic3/Step5/pseudo_bulk_celltype_groups/PDvControlbulk/info/DGE_",i,"MainContrast HC vs PD.csv"), header = T, sep = ","))
    }

    ## create data frame for each cell type
    list2env(list,envir=.GlobalEnv)

### 1.5. Clean up cell-based DEG files

    ## make list of cell type dataframes
    CellType_df <- list(Astrocytes, DaN, Endothelial, Ependymal, Excitatory, GABA, Inhibitory, Microglia, Oligocendrocytes, OPC, Pericytes, all)

    ## set names of list as cell types
    names(CellType_df) <- c('Astrocytes', 'DaN', 'Endothelial', 'Ependymal', 'Excitatory', 'GABA', 'Inhibitory', 'Microglia', 'Oligocendrocytes', 'OPC', 'Pericytes', "all")

    ## list of columns that we want to keep
    cols_keep <- c("X","p_val","avg_log2FC","p_val_adj")

    ## process each dataframe to only keep DEGs with a p-value < 0.05 and log2FC > 1 or < -1
    list <- list()

    for (i in unique(names(CellType_df))){
      j <- data.frame(CellType_df[i])
      colnames(j) <- c("X", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")
      j <- j[,colnames(j) %in% cols_keep]
      j <- subset(j, p_val <= 0.05)
      j <- subset(j, avg_log2FC >= 1 | avg_log2FC <= -1)
      if (nrow(j > 0)){
      j$cell_type <- paste0(i)
      j$method <- "cell_based"
      list[i] <- list(j)
      }else{
        print(paste0("did not find any DEGs for ", i))
      }
      list[i]
    }

    ## create data frame for each cell type
    list2env(list,envir=.GlobalEnv)

### 1.6. Clean up sample-based DEG files

    ## make list of cell type dataframes
    CellType_df <- list(astro_sample, DaN_sample, endo_sample, ependymal_sample, excit_sample, GABA_sample, inhib_sample, mg_sample, Olig_sample, opc_sample, peri_sample, all_sample)

    ## set names of list as cell types
    names(CellType_df) <- c('astro_sample', 'DaN_sample', 'endo_sample', 'ependymal_sample', 'excit_sample', 'GABA_sample', 'inhib_sample', 'mg_sample', 'Olig_sample', 'opc_sample', 'peri_sample', "all_sample")

    ## list of columns that we want to keep
    cols_keep <- c("X","pvalue","log2FoldChange","padj")

    ## process each dataframe to only keep DEGs with a p-value and log2FC > 1 or < -1
    list <- list()

    for (i in unique(names(CellType_df))){
      j <- data.frame(CellType_df[i])
      colnames(j) <- c('X', 'baseMean', 'log2FoldChange', 'lfcSE' ,'stat', 'pvalue', 'padj')
      j <- j[,colnames(j) %in% cols_keep]
      j <- subset(j, pvalue <= 0.05)
      j <- subset(j, log2FoldChange >= 1 | log2FoldChange <= -1)
      if (nrow(j > 0)){
      j$cell_type <- paste0(i)
      j$method <- "sample_based"
      list[i] <- list(j)
      }else{
        print(paste0("did not find any DEGs for ", i))
      }
      list[i]
    }

    ## create data frame for each cell type
    list2env(list,envir=.GlobalEnv)

### 1.7. Compute number of DEGs identified for each type using cell-based DGE analysis

    ## bind dataframes
    # do not include any cell types that identified 0 DEGS (e.g. OPC)
    bind_cell_based <- rbind(Astrocytes, DaN, Endothelial, Ependymal, Excitatory, GABA, Inhibitory, Microglia, Oligocendrocytes,Pericytes, all)

    ## define colours
    bind_cell_based$col[bind_cell_based$avg_log2FC > 0] <- "indianred3"
    bind_cell_based$col[bind_cell_based$avg_log2FC < 0] <- "dodgerblue2"
    bind_cell_based$col[bind_cell_based$avg_log2FC > 0 & bind_cell_based$p_val_adj < 0.05] <- "red4"
    bind_cell_based$col[bind_cell_based$avg_log2FC < 0 & bind_cell_based$p_val_adj < 0.05] <- "navy"
    bind_cell_based$direction[bind_cell_based$avg_log2FC > 0] <- "up"
    bind_cell_based$direction[bind_cell_based$avg_log2FC < 0] <- "down"

    ## plot
    cell_based_DEG_counts <- ggplot(bind_cell_based, aes(x = direction, fill = col)) +
      geom_bar() +
      theme_classic() + 
      theme(axis.text.x = element_blank(),
            axis.text.y = element_text(size = 12),
            axis.ticks.x = element_blank(),
            strip.background = element_rect(colour=NA, fill=NA),
            strip.text = element_text(size = 12),
            legend.position = "none",
            axis.title.y = element_text(face="bold", size = 12),
            axis.title.x = element_text(face="bold", size = 12),
            plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) +
      scale_fill_manual(values = c("dodgerblue2", "indianred3", "navy","red4")) +
      facet_wrap(~cell_type,  ncol=12, strip.position = "bottom") +
      ylab("Number of DEGs") +
      xlab("Cell type") + 
      scale_y_continuous(limits = c(0, 60)) +
      ggtitle ("Cell-based: MAST")
    cell_based_DEG_counts

<p align="center">
<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/0fc6bd1e-89de-400b-8970-ff2e404d43ac"> <br /> 
</p>

**Figure 1. Number of DEGs identified by cell-based DGE analysis for each cell type.** DEGs were calculated using cells as replicates and the MAST framework. Bar chart showing the number of DEGs identified with a log 2 fold-change < -1 (blue) and > 1 (red) and p-values < 0.05. Bonferroni adjusted p-values < 0.05 are indicated by the darker shade.



### 1.8. Compute number of DEGs identified for each cell type using sample-based DGE analysis

    ## bind dataframes
    bind_sample_based <- rbind(astro_sample, DaN_sample, endo_sample, ependymal_sample, excit_sample, GABA_sample, inhib_sample, mg_sample, Olig_sample, opc_sample, peri_sample, all_sample)

    ## define colours
    bind_sample_based$col[bind_sample_based$log2FoldChange > 0] <- "indianred3"
    bind_sample_based$col[bind_sample_based$log2FoldChange < 0] <- "dodgerblue2"
    bind_sample_based$col[bind_sample_based$log2FoldChange > 0 & bind_sample_based$padj < 0.05] <- "red4"
    bind_sample_based$col[bind_sample_based$log2FoldChange < 0 & bind_sample_based$padj < 0.05] <- "navy"
    bind_sample_based$direction[bind_sample_based$log2FoldChange > 0] <- "up"
    bind_sample_based$direction[bind_sample_based$log2FoldChange < 0] <- "down"

    ## rename cell types
    bind_sample_based$cell_type[bind_sample_based$cell_type == "astro_sample"] <- "Astrocyte"
    bind_sample_based$cell_type[bind_sample_based$cell_type == "DaN_sample"] <- "DaN"
    bind_sample_based$cell_type[bind_sample_based$cell_type == "endo_sample"] <- "Endothelial"
    bind_sample_based$cell_type[bind_sample_based$cell_type == "ependymal_sample"] <- "Ependymal"
    bind_sample_based$cell_type[bind_sample_based$cell_type == "excit_sample"] <- "Excitatory"
    bind_sample_based$cell_type[bind_sample_based$cell_type == "GABA_sample"] <- "GABA"
    bind_sample_based$cell_type[bind_sample_based$cell_type == "inhib_sample"] <- "Inhibitory"
    bind_sample_based$cell_type[bind_sample_based$cell_type == "mg_sample"] <- "Microglia"
    bind_sample_based$cell_type[bind_sample_based$cell_type == "Olig_sample"] <- "Oligoden."
    bind_sample_based$cell_type[bind_sample_based$cell_type == "opc_sample"] <- "OPC"
    bind_sample_based$cell_type[bind_sample_based$cell_type == "peri_sample"] <- "Pericyte"
    bind_sample_based$cell_type[bind_sample_based$cell_type == "all_sample"] <- "All Cells"


    ## plot
    sample_based_DEG_counts <- ggplot(bind_sample_based, aes(x = direction, fill = col)) +
      geom_bar() +
      theme_classic() + 
      theme(axis.text.x = element_blank(),
            axis.text.y = element_text(size = 12),
            axis.ticks.x = element_blank(),
            strip.background = element_rect(colour=NA, fill=NA),
            strip.text = element_text(size = 12),
            legend.position = "none",
            axis.title.y = element_text(face="bold", size = 12),
            axis.title.x = element_text(face="bold", size = 12),
            plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) +
      scale_fill_manual(values = c("dodgerblue2", "indianred3", "navy","red4")) +
      facet_wrap(~cell_type,  ncol=12, strip.position = "bottom") +
      ylab("Number of DEGs") +
      xlab("Cell type") + 
      ggtitle ("Sample-based: DESeq2")
    sample_based_DEG_counts

<p align="center">
<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/86e4f2dc-6cd7-4b96-bbab-5a60193ad744"> <br /> 
</p>

**Figure 2. Number of DEGs identified by sample-based DGE analysis for each cell type.**. DEGs were calculated using samples as replicates and the DESeq2 framework. Bar chart showing the number of DEGs identified with a log 2 fold-change < -1 (blue) and > 1 (red) and p-values < 0.05. Bonferroni adjusted p-values < 0.05 are indicated by the darker shade.

### 1.9. DEG overlap analysis

    ## list of cell-based dataframes
    cell_based_dfs <- list(Astrocytes, DaN, Endothelial, Ependymal, Excitatory, GABA, Inhibitory, Microglia, Oligocendrocytes, OPC, Pericytes, all)

    ## set names of named list
    names(cell_based_dfs) <- c('Astrocytes_CellBased', 'DaN_CellBased', 'Endothelial_CellBased', 'Ependymal_CellBased', 'Excitatory_CellBased', 'GABA_CellBased', 'Inhibitory_CellBased', 'Microglia_CellBased', 'Oligocendrocytes_CellBased', 'OPC_CellBased', 'Pericytes_CellBased', 'all_CellBased')

    ## list of sample-based dataframes
    sample_based_dfs <- list(astro_sample, DaN_sample, endo_sample, ependymal_sample, excit_sample, GABA_sample, inhib_sample, mg_sample, Olig_sample, opc_sample, peri_sample, all_sample)

    ## set names of named list
    names(sample_based_dfs) <- c('Astrocytes_SampleBased', 'DaN_SampleBased', 'Endothelial_SampleBased', 'Ependymal_SampleBased', 'Excitatory_SampleBased', 'GABA_SampleBased', 'Inhibitory_SampleBased', 'Microglia_SampleBased', 'Oligocendrocytes_SampleBased', 'OPC_SampleBased', 'Pericytes_SampleBased', 'all_SampleBased')

    ## combine the list
    combined_list <- append(cell_based_dfs,sample_based_dfs )

    ## list of unique cell types
    cell_types <- unique(str_extract(names(combined_list), "[^_]+"))

    ## compute overlap 
    list <- list()

    for (i in unique(cell_types)){
      df <- data.frame(names(combined_list))
      df$num <- rownames(df)
      df_lim <- df[grep(i, df$names.combined_list.), ]
      
      cell_num <- as.numeric(df_lim$num[grep("CellBased", df_lim$names.combined_list.)])
      sample_num <- as.numeric(df_lim$num[grep("SampleBased", df_lim$names.combined_list.)])

      cell <- data.frame(combined_list[as.numeric(cell_num)])
      colnames(cell) <- c('X', 'p_val', 'avg_log2FC', 'p_val_adj', 'cell_type', 'method')
      sample <- data.frame(combined_list[as.numeric(sample_num)])
      colnames(sample) <- c("X", "log2FoldChange", "pvalue", "padj", "cell_type", "method" )

      # both padjusted
      cell_temp <- cell$X[cell$p_val_adj <= 0.05]
      cell_temp <-  cell_temp[!is.na(cell_temp)]
      sample_temp  <- sample$X[sample$padj <= 0.05]
      sample_temp <-  sample_temp[!is.na(sample_temp)]
      length_adj <- length(intersect(sample_temp, cell_temp))
      
      # both pvalue
      cell_temp <- cell$X[cell$p_val <= 0.05]
      cell_temp <-  cell_temp[!is.na(cell_temp)]
      sample_temp  <- sample$X[sample$pvalue <= 0.05]
      sample_temp <-  sample_temp[!is.na(sample_temp)]
      length_p <- length(intersect(sample_temp, cell_temp))

      ## cell-based pvalue
      cell_temp <- cell$X
      cell_temp <-  cell_temp[!is.na(cell_temp)]
      length_cell_base <- length(cell_temp) 
      
      ## sample-based pvalue
      sample_temp <- sample$X
      sample_temp <-  sample_temp[!is.na(sample_temp)]
      length_sample_base <- length(sample_temp)

      ## cell-based padjust
      cell_temp <- cell$X[cell$p_val_adj <= 0.05]
      cell_temp <-  cell_temp[!is.na(cell_temp)]
      length_cell_adj <- length(cell_temp)
      
      ## sample-based padjust
      sample_temp <- sample$X[sample$padj <= 0.05]
      sample_temp <-  sample_temp[!is.na(sample_temp)]
      length_sample_adj <- length(sample_temp)

      ## dataframe
      class <- c('AdjustedP_Overlap', 'P_Overlap', 'P_CellBased', 'P_SampleBased', 'AdjustedP_CellBased',    'AdjustedP_SampleBased')
      val <- c(length_adj,length_p,length_cell_base,length_sample_base, length_cell_adj,length_sample_adj)
      temp_df <- data.frame(class, val)
      temp_df$celltype <- i
      list[i] <- list(temp_df)
    } 

    ## set list elements to dataframe
    list2env(list,envir=.GlobalEnv)

    ## bind dataframe
    df <- bind_rows(mget(cell_types))
      
    ## set fator level
    df$class <- factor(df$class, levels = c("P_CellBased", "P_SampleBased", "P_Overlap", "AdjustedP_CellBased", "AdjustedP_SampleBased", "AdjustedP_Overlap"))

    ## plot
    ggplot(df, aes(x = celltype, y = class, fill = val)) +geom_tile(col = "black", fill = "white") + theme_bw()+ geom_text(aes(label = val))+
      theme( axis.text = element_text(size = 12),
             axis.title = element_text(size = 12, face = "bold"),
             legend.title = element_text(face = "bold"),
             axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
             axis.text.y = element_text(colour = "black"))+
      scale_colour_manual(values = c("black", "black"), legend) +
      scale_y_discrete(expand = c(0,0)) + #labels = c("Pseudo-bulk: p-value < 0.05",
      scale_x_discrete(expand = c(0,0)) +
      xlab("Cell Type") +
      ylab("")  + guides(colour = "none")

<p align="center">
<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/ec4bc949-776c-4dff-9d80-99ffa8c08c7c"> <br /> 
</p>

**Figure 3. Number of DEG identified by cell-based DGE-MAST, sample-based DGE-DESeq2, or both frameworks across all cell types**. Only DEGs with olg 2 fold-change < -1 and > 1 are included.

- - - -
## 2. Cell-based DGE enrichment analysis

### 2.1. load cell-based all cells DEG file

    ## load DEG csv file
    all <- read.delim("/Users/mfiorini/Desktop/scRNA_pipeline/Manuscript/Manuscript_Figures/Smajic3/Step5/wilcoxon_all_cells/HCvPD/HCvPD_DEG.csv", header = T, sep = ",")

### 2.2. load cell-based cell type groups DEG files

    ## list of cell types used for DEG contrasts
    cell_based_cell_types <- c('Astrocytes', 'DaN', 'Endothelial', 'Ependymal', 'Excitatory', 'GABA', 'Inhibitory', 'Microglia', 'Oligocendrocytes', 'OPC', 'Pericytes')

    ## create empty list object
    list <- list()

    ## load DEG csv files
    for(i in unique(cell_based_cell_types)){
      list[i] <- list(read.delim(paste0("/Users/mfiorini/Desktop/scRNA_pipeline/Manuscript/Manuscript_Figures/Smajic3/Step5/wilcoxon_celltype_groups/", i,"PDvHC/",i,"PDvHC_DEG.csv"), header = T, sep = ","))
    }

    ## create data frame for each cell type
    list2env(list,envir=.GlobalEnv)

### 2.3. Clean up cell-based DEG files

    ## make list of cell type dataframe
    CellType_df <- list(Astrocytes, DaN, Endothelial, Ependymal, Excitatory, GABA, Inhibitory, Microglia, Oligocendrocytes, OPC, Pericytes, all)

    ## set names of list as cepp types
    names(CellType_df) <- c('Astrocytes', 'DaN', 'Endothelial', 'Ependymal', 'Excitatory', 'GABA', 'Inhibitory', 'Microglia', 'Oligocendrocytes', 'OPC', 'Pericytes', "all")

    ## list of columns that we want to keep
    cols_keep <- c("X","p_val","avg_log2FC","p_val_adj")

    ## process each dataframe to only keep DEGs with a p-value < 0.05 and log2FC > 1
    list_CellBased <- list()

    for (i in unique(names(CellType_df))){
      j <- data.frame(CellType_df[i])
      colnames(j) <- c("X", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")
      j <- j[,colnames(j) %in% cols_keep]
      j <- subset(j, p_val <= 0.05)
      j <- subset(j, avg_log2FC >= 1 | avg_log2FC <= -1)
      if (nrow(j > 0)){
      j$cell_type <- paste0(i)
      j$method <- "cell_based"
      list_CellBased[i] <- list(j)
      }else{
        print(paste0("did not find any DEGs for ", i))
      }
      list_CellBased[i]
    }

    ## create data frame for each cell type
    list2env(list_CellBased,envir=.GlobalEnv)

### 2.4. EnrichR: Cell type groups

    ## do not run enrichr for all cell types
    test <- data.frame(names(list_CellBased))
    test <- test %>% 
      filter(!grepl('all', names.list_CellBased.))
    cell_type_num <- as.numeric(rownames(test))

    ## create empty list
    CellBased_gse_list <- list()

    ## GSE loop
    for (i in cell_type_num){
      ## create cell type specific dataframe
      j <- data.frame(list_CellBased[i])
      colnames(j) <- c('X', 'p_val', 'avg_log2FC', 'p_val_adj', 'cell_type', 'method')
      ## we want the log2 fold change 
      gene_list <- j$avg_log2FC
      ## assign gene to log2fc
      names(gene_list) <- j$X
      ## omit any NA values 
      gene_list<-na.omit(gene_list)
      ## sort the list in decreasing order
      gene_list = sort(gene_list, decreasing = TRUE)
      ## perform gseGO
      gse_bulk<- gseGO(geneList=gene_list, 
                                ont ="BP", 
                                keyType = "SYMBOL", 
                                nPerm = 1000, 
                                minGSSize = 3, 
                                maxGSSize = 800, 
                                pvalueCutoff = 1.00, 
                                verbose = TRUE, 
                                OrgDb = organism, 
                                pAdjustMethod = "none")
      ## organize results
      if (nrow(gse_bulk@result) > 0){
      ## count the gene number
      gene_count<- gse_bulk@result %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
      ## merge with the original dataframe
      dot_df<- left_join(gse_bulk@result, gene_count, by = "ID") %>% mutate(GeneRatio = count/setSize)
      dot_df <- dot_df[order(dot_df$pvalue),]
      dot_df$Celltype <- names(list_CellBased[i])
      ## store in list 
      CellBased_gse_list[i] <- list(dot_df)
      }else{
        print(paste0('no term enriched under specific pvalueCutoff for ',names(list_CellBased[i])))
      }
    }

    ## remove NULL (cell types with no enrichment results)
    CellBased_gse_list <- CellBased_gse_list[!sapply(CellBased_gse_list,is.null)]

    ## set names of list element to cell types
    names <- list()
    for (i in 1:length(CellBased_gse_list)){
    j <- data.frame(CellBased_gse_list[i])
    names[i] <- unique(j$Celltype)
    }
    names(CellBased_gse_list) <- names

    ## set to dataframe
    list2env(CellBased_gse_list,envir=.GlobalEnv)


### 2.5. EnrichR: All cells

    ## retrieve logFC
    original_gene_list <- all$avg_log2FC
    ## assign logFC to gene
    names(original_gene_list) <- all$X
    ## omit any NA values 
    gene_list<-na.omit(original_gene_list)
    ## sort the list in decreasing order 
    gene_list = sort(gene_list, decreasing = TRUE)
    ## perform gseGO
    gse_all<- gseGO(geneList=gene_list, 
                     ont ="BP", 
                     keyType = "SYMBOL", 
                     nPerm = 1000, 
                     minGSSize = 3, 
                     maxGSSize = 800, 
                     pvalueCutoff = 1.00, 
                     verbose = TRUE, 
                     OrgDb = organism, 
                     pAdjustMethod = "none")


    ## count the gene number
    gene_count<- gse_all@result %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)

    ## merge with the original dataframe
    dot_df<- left_join(gse_all@result, gene_count, by = "ID") %>% mutate(GeneRatio = count/setSize)

    ## Take 5 most significantly enriched terms
    dot_df <- dot_df[order(dot_df$pvalue),]
    all_terms <- dot_df$Description[c(1:5)]
    dot_df_all <- dot_df
    dot_df_all$Celltype <- "all"

### 2.6. Plot GSE result for cell-based DEG

    ## bind cell type dataframes
    namesX <- unlist(names)
    df <- bind_rows(mget(namesX))

    ## bind celltype dataframe with all cells
    df_total <- rbind(df, dot_df_all)

    ## subset to only inlclude top GSE terms from all cells
    df_total <- subset(df_total, Description %in% all_terms)
    length(unique(df_total$Description))

    ## [1] 5

    ## rename cell type names
    df_total$Celltype[df_total$Celltype == "all"] <- "All cells"
    df_total$Celltype[df_total$Celltype == "Oligodendrocytes"] <- "Oligoden."

    ## plot
    gse_CellBased<- ggplot(df_total, aes(x = GeneRatio, y=Description)) + 
      geom_point(aes(size = count, color = pvalue)) +
      theme_bw() +
      theme(strip.background = element_rect(colour=NA, fill=NA),
            strip.text = element_text(size = 12),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.title.x = element_text(size = 12, face = "bold"),
            axis.text.y = element_text(size = 10),
            axis.text.x = element_text(size = 12),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12)) +
      scale_colour_gradient( low="red", limits=c(0.0,1)) +
      facet_wrap(~Celltype,  ncol=12) +
      ylab(NULL) +
      ggtitle("Cell-based: MAST")+
      scale_x_continuous(limits=c(0.1,1), labels = c("0.0","", "0.5","", "1.0"))
    gse_CellBased

<p align="center">
<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/9ffb6fb9-6cec-40ba-9ded-740dbba20c13"> <br /> 
</p>

**Figure 4. Cell-based DGE: Enrichment of the top 5 GO terms for GO-Biological Processes calculated for all cell types together across cell types.**  DEGs with p-values < 0.05 and log 2 fold-change < -1 and > 1 were used as the input for gene set enrichment analysis. The gene ratio, gene count, and p-value of the five terms in each cell type are shown.

- - - -

## 3. Sample-based DGE enrichment analysis

### 3.1. load sample-based all cells DEG file

    ## load DEG csv file
    all_sample <- read.delim("/Users/mfiorini/Desktop/scRNA_pipeline/Manuscript/Manuscript_Figures/Smajic3/Step5/pseudo_bulk_all_cells/PDvControl/DGE_AllCellsMainContrast HC vs PD.csv", header = T, sep = ",")

### 3.2. load sample-based cell type groups DEG files

    ## list of cell types used for DEG contrasts
    sample_based_cell_types <- c('astro', 'DaN', 'endo', 'ependymal', 'excit', 'GABA', 'inhib', 'mg', 'Olig', 'opc', 'peri')

    ## create empty list object
    list <- list()

    ## load DEG csv files
    for(i in unique(sample_based_cell_types)){
      list[paste0(i,"_sample")] <- list(read.delim(paste0("/Users/mfiorini/Desktop/scRNA_pipeline/Manuscript/Manuscript_Figures/Smajic3/Step5/pseudo_bulk_celltype_groups/PDvControlbulk/info/DGE_",i,"MainContrast HC vs PD.csv"), header = T, sep = ","))
    }

    ## create data frame for each cell type
    list2env(list,envir=.GlobalEnv)


### 3.3. Clean up sample-based DEG files

    ## make list of cell type dataframe
    CellType_df <- list(astro_sample, DaN_sample, endo_sample, ependymal_sample, excit_sample, GABA_sample, inhib_sample, mg_sample, Olig_sample, opc_sample, peri_sample, all_sample)

    ## set names of list as cepp types
    names(CellType_df) <- c('astro_sample', 'DaN_sample', 'endo_sample', 'ependymal_sample', 'excit_sample', 'GABA_sample', 'inhib_sample', 'mg_sample', 'Olig_sample', 'opc_sample', 'peri_sample', "all_sample")

    ## list of columns that we want to keep
    cols_keep <- c("X","pvalue","log2FoldChange","padj")

    ## process each dataframe to only keep DEGs with a p-value < 0.05 log2FC > 1
    list_SampleBased <- list()

    for (i in unique(names(CellType_df))){
      j <- data.frame(CellType_df[i])
      colnames(j) <- c('X', 'baseMean', 'log2FoldChange', 'lfcSE' ,'stat', 'pvalue', 'padj')
      j <- j[,colnames(j) %in% cols_keep]
      j <- subset(j, pvalue <= 0.05)
      j <- subset(j, log2FoldChange >= 1 | log2FoldChange <= -1)
      if (nrow(j > 0)){
      j$cell_type <- paste0(i)
      j$method <- "sample_based"
      list_SampleBased[i] <- list(j)
      }else{
        print(paste0("did not find any DEGs for ", i))
      }
      list_SampleBased[i]
    }

    ## create data frame for each cell type
    list2env(list_SampleBased,envir=.GlobalEnv)

    ## <environment: R_GlobalEnv>

### 3.4. EnrichR: Cell type groups

    ## do not run enrichr for all cell types
    test <- data.frame(names(list_SampleBased))
    test <- test %>% 
      filter(!grepl('all', list_SampleBased))
    cell_type_num <- as.numeric(rownames(test))

    ## create empty list
    SampleBased_gse_list <- list()

    ## GSE loop
    for (i in cell_type_num){
      ## create cell type specific dataframe
      j <- data.frame(list_SampleBased[i])
      colnames(j) <- c('X', 'log2FoldChange', 'pvalue', 'padj', 'cell_type', 'method')
      ## we want the log2 fold change 
      gene_list <- j$log2FoldChange
      ## assign gene to log2fc
      names(gene_list) <- j$X
      ## omit any NA values 
      gene_list<-na.omit(gene_list)
      ## sort the list in decreasing order
      gene_list = sort(gene_list, decreasing = TRUE)
      ## perform gseGO
      gse_bulk<- gseGO(geneList=gene_list, 
                                ont ="BP", 
                                keyType = "SYMBOL", 
                                nPerm = 1000, 
                                minGSSize = 3, 
                                maxGSSize = 800, 
                                pvalueCutoff = 1.00, 
                                verbose = TRUE, 
                                OrgDb = organism, 
                                pAdjustMethod = "none")
      ## organize results
      if (nrow(gse_bulk@result) > 0){
      ## count the gene number
      gene_count<- gse_bulk@result %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
      ## merge with the original dataframe
      dot_df<- left_join(gse_bulk@result, gene_count, by = "ID") %>% mutate(GeneRatio = count/setSize)
      dot_df <- dot_df[order(dot_df$pvalue),]
      dot_df$Celltype <- names(list_SampleBased[i])
      ## store in list 
      SampleBased_gse_list[i] <- list(dot_df)
      }else{
        print(paste0('no term enriched under specific pvalueCutoff for ',names(list_SampleBased[i])))
      }
    }



    ## remove NULL (cell types with no enrichment results)
    SampleBased_gse_list <- SampleBased_gse_list[!sapply(SampleBased_gse_list,is.null)]

    ## set names of list element to cell types
    names <- list()
    for (i in 1:length(SampleBased_gse_list)){
    j <- data.frame(SampleBased_gse_list[i])
    names[i] <- unique(j$Celltype)
    }
    names(SampleBased_gse_list) <- names

    ## set to dataframe
    list2env(SampleBased_gse_list,envir=.GlobalEnv)

    ## <environment: R_GlobalEnv>

### 3.5. EnrichR: All cells

    ## retrieve logFC
    original_gene_list <- all_sample$log2FoldChange
    ## assign logFC to gene
    names(original_gene_list) <- all_sample$X
    ## omit any NA values 
    gene_list<-na.omit(original_gene_list)
    ## sort the list in decreasing order 
    gene_list = sort(gene_list, decreasing = TRUE)
    ## perform gseGO
    gse_all<- gseGO(geneList=gene_list, 
                     ont ="BP", 
                     keyType = "SYMBOL", 
                     nPerm = 1000, 
                     minGSSize = 3, 
                     maxGSSize = 800, 
                     pvalueCutoff = 1.00, 
                     verbose = TRUE, 
                     OrgDb = organism, 
                     pAdjustMethod = "none")


    ## count the gene number
    gene_count<- gse_all@result %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)

    ## merge with the original dataframe
    dot_df<- left_join(gse_all@result, gene_count, by = "ID") %>% mutate(GeneRatio = count/setSize)

    ## Take 5 most significantly enriched terms
    dot_df <- dot_df[order(dot_df$pvalue),]
    all_terms <- dot_df$Description[c(1:5)]
    dot_df_all <- dot_df
    dot_df_all$Celltype <- "all"

### 3.6. Plot GSE result for cell-based DEG

    ## bind cell type dataframes
    namesX <- unlist(names)
    df <- bind_rows(mget(namesX))

    ## bind celltype dataframe with all cells
    df_total <- rbind(df, dot_df_all)

    ## subset to only inlclude top GSE terms from all cells
    df_total <- subset(df_total, Description %in% all_terms)
    length(unique(df_total$Description))

    ## rename cell type names
    df_total$Celltype[df_total$Celltype == "all"] <- "All cells"
    df_total$Celltype[df_total$Celltype == "DaN_sample"] <- "DaN"
    df_total$Celltype[df_total$Celltype == "endo_sample"] <- "Endothelial"
    df_total$Celltype[df_total$Celltype == "ependymal_sample"] <- "Ependymal"
    df_total$Celltype[df_total$Celltype == "excit_sample"] <- "Excitatory"
    df_total$Celltype[df_total$Celltype == "GABA_sample"] <- "GABA"
    df_total$Celltype[df_total$Celltype == "inhib_sample"] <- "Inhibitory"
    df_total$Celltype[df_total$Celltype == "mg_sample"] <- "Microglia"
    df_total$Celltype[df_total$Celltype == "Olig_sample"] <- "Oligoden."
    df_total$Celltype[df_total$Celltype == "opc_sample"] <- "OPC"
    df_total$Celltype[df_total$Celltype == "peri_sample"] <- "Pericytes"
    df_total$Celltype[df_total$Celltype == "astro_sample"] <- "Astrocyte"

    ## plot
    gse_CellBased<- ggplot(df_total, aes(x = GeneRatio, y=Description)) + 
      geom_point(aes(size = count, color = pvalue)) +
      theme_bw() +
      theme(strip.background = element_rect(colour=NA, fill=NA),
            strip.text = element_text(size = 12),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.title.x = element_text(size = 12, face = "bold"),
            axis.text.y = element_text(size = 10),
            axis.text.x = element_text(size = 12),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12)) +
      scale_colour_gradient( low="red", limits=c(0.0,1)) +
      facet_wrap(~Celltype,  ncol=12) +
      ylab(NULL) +
      ggtitle("Sample-based: DESeq2")+
      scale_x_continuous(limits=c(0.1,1), labels = c("0.0","", "0.5","", "1.0"))
    gse_CellBased

<p align="center">
<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/1e9126f3-cb17-4c80-8850-c5b457936370"> <br /> 
</p>

**Figure 5. Sample-based DGE: Enrichment of the top 5 GO terms for GO-Biological Processes calculated for all cell types together across cell types.**  DEGs with p-values < 0.05 and log 2 fold-change < -1 and > 1 were used as the input for gene set enrichment analysis. The gene ratio, gene count, and p-value of the five terms in each cell type are shown.

