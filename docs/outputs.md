## Outputs of each Analytical Step in the scRNAbox pipeline

- [Standard scRNAseq Analysis Track](#standard-scrnaseq-analysis-track)
- [Cell Hashtag scRNAseq Analysis Track](#cell-hashtag-scrnaseq-analysis-track)

## Standard scRNAseq Analysis Track
#### Step 1: FASTQ pre-processing
All of the outputs of the CellRanger _counts_ pipeline are produced. For more information on the outputs, please visit the CellRanger [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count).

#### Step 2: Ambient RNA removal and create Seurat object
|Output type |Name|Description|
|:--|:--|:--|
|Figure|vioplot_sample_name.png | Sample-specific violin plot showing the distribution of cells according to QC metrics|
|Info|sample_name_RNA.txt|Sample-specific sparse matrix of RNA assay|
|Info|MetaData_sample_name.txt|Sample-specific dataframe showing the Seurat object metadata|
|Info|meta_info_sample_name.txt|Sample-specific text file showing the column names of the Seurat object metadata|
|Info|summary_sample_name.txt|Sample-specific text file showing the summary of QC metrics (Minimum, 1st Quartile, Median, Mean, 3rd Quartile, Maximum) | 
|Info|sessionInfo.txt|Session information for the R session|
|Data object|sample_name.rds|Sample-specific intermediate Seurat RDS object|


#### Step 3: Quality control and filtering
|Output type |Name|Description|
|:--|:--|:--|
|Figure|cellcycle_sample_name.png | Sample-specific violin plot showing the distribution of cells according to cell cycle S and G2M scores|
|Figure|dimplot_pca_sample_name.png | Sample-specific PCA showing the first two PCs|
|Figure|dimplot_umap_sample_name.png | Sample-specific UMAP|
|Figure|QC_vioplot_sample_name.png | Sample-specific violin plot showing the distribution of cells according to QC metrics after filtering|
|Figure|VariableFeaturePlot_sample_name.png | Sample-specific figure showing the most variably expressed genes|
|Info|sample_name_RNA.txt|Sample-specific sparse matrix of RNA assay|
|Info|MetaData_sample_name.txt|Sample-specific dataframe showing the Seurat object metadata|
|Info|meta_info_sample_name.txt|Sample-specific text file showing the column names of the Seurat object metadata|
|Info|most_variable_genes_sample_name.txt|Sample-specific text file showing the column names of the Seurat object metadata|
|Info|summary_sample_name.txt|Sample-specific text file showing the summary of QC metrics (Minimum, 1st Quartile, Median, Mean, 3rd Quartile, Maximum) | 
|Info|sessionInfo.txt|Session information for the R session|
|Data object|sample_name.rds|Sample-specific intermediate Seurat RDS object|

#### Step 4: Doublet removal
|Output type |Name|Description|
|:--|:--|:--|
|Figure|sample_nameDF.classifications.png | Sample-specific UMAP plot showing droplet classifications (singlet or doublet)|
|Info|sample_name_RNA.txt|Sample-specific sparse matrix of RNA assay|
|Info|MetaData_sample_name.txt|Sample-specific dataframe showing the Seurat object metadata|
|Info|meta_info_sample_name.txt|Sample-specific text file showing the column names of the Seurat object metadata|
|Info|sessionInfo.txt|Session information for the R session|
|Data object|sample_name.rds|Sample-specific intermediate Seurat RDS object|


#### Step 5: Integration and linear dimensional reduction
|Output type |Name|Description|
|:--|:--|:--|
|Figure|DimPlot_pca.png | PCA showing the first two PCs, colour-coded by sample|
|Figure|DimPlot_umap.png | UMAP, colour-coded by sample|
|Figure|elbow.png | Elbow plot to visualize the percentage of variance explained by each PC |
|Figure|Jackstraw_plot.png| Jackstraw plot to visualize the distribution of p-values for each PC|
|Info|seu_int_RNA.txt|Sparse matrix of integrated assay|
|Info|seu_int_MetaData.txt|Dataframe showing the Seurat object metadata|
|Info|meta_info_seu_step5.csv|Text file showing the column names of the Seurat object metadata|
|Info|sessionInfo.txt|Session information for the R session|
|Data object|seu_step5.rds|Integrated intermediate Seurat RDS object|

#### Step 6: Clustering
|Output type |Name|Description|
|:--|:--|:--|
|Figure|clustree_int.png  | Clustree plot showing the stability across the user-defied clustering resolutions|
|Figure|integrated_snn_res.png| UMAP at the user defined clustering-resolution|
|Figure|ARI.png | Mean and standard deviation of the Adjusted Rand Index (ARI) between clustering pairs at a user-defined resolution |
|Info|**test.xlsx**| Excel file showing the mean and standard deviation of the ARI between clustering pairs at a user-defined resolution |
|Info|seu_RNA.txt|Sparse matrix of integrated assay|
|Info|seu_MetaData.txt|Dataframe showing the Seurat object metadata|
|Info|meta_info.csv |Text file showing the column names of the Seurat object metadata|
|Info|sessionInfo.txt|Session information for the R session|
|Data object|seu_step6.rds|Intermediate Seurat RDS object|


#### Step 7: Cluster annotation
Cluster annotation method|Output type |Name|Description|
|:--|:--|:--|:--|
|General|Figure|umap.pdf | UMAP plot of integrated assay at the user-defined clustering resolution used for cluster annotation|
|General|Figure|umap_splitted.pdf | UMAP plot of integrated assay at the user-defined clustering resolution used for cluster annotation, split by sample|
|Cluster marker GSEA (Method 1)|Figure|heatmap.pdf | Heatmap showing the expression of the top marker genes across cells, stratified by cluster|
|Cluster marker GSEA (Method 1)|Figure|plotenrich.pdf | Barplot showing the 20 most enriched terms for a particular cluster and cell type library|
|Module score (Method 2)|Figure|module_score_gene_set.png | UMAP plot showing the module score across cells for user-defined gene sets|
|Reference-based annotations (Method 3)|Figure|UMAP_transferred_labels.pdf |UMAP plots showing the cluster annotations from the reference Seurat object projected onto the query Seurat object|
|Visualize select features|Figure|select_feature_dot_plot.pdf|Dotplot showing the expression of user-defined features at the cluster level|
|Visualize select features|Figure|select_feature_violin_plot.pdf |Violin plot showing the expression of user-defined features at the cluster level|
|Visualize select features|Figure|select_feature_feature_plot.pdf  |UMAP plots showing the expression of user-defined features at the cell level|
|General|Info|meta_info_seu_step7.txt |Text file showing the column names of the Seurat object metadata|
|General|Info|sessionInfo_marker.txt|Session information for the R session|
|Cluster marker GSEA (Method 1)|Info|cluster_just_genes.xlsx  |Excel file showing the marker genes for each cluster|
|Cluster marker GSEA (Method 1)|Info|cluster_whole.xlsx  |Excel file showing the marker genes and corresponding summary statistics for each cluster|
|Cluster marker GSEA (Method 1)|Info|ClusterMarkers.csv  |csv file showing the marker genes and corresponding summary statistics for each cluster|
|Cluster marker GSEA (Method 1)|Info| top_sel.csv |csv file showing the top n marker genes for each cluster. The user defined n in the execution parameters|
|Cluster marker GSEA (Method 1)|Info| Er.genes.csv|Enrichment terms and the corresponding statistics for a particular cluster and cell type library|
|Module score (Method 2)|Info| geneset_by_cluster.csv|Mean module score across clusters for each user-defined gene set|
|General|Data object|seu_step7.rds|Intermediate Seurat RDS object|
|Cluster marker GSEA (Method 1)|Data object| ClusterMarkers.rds |RDS object containing the marker genes for each cluster|

#### Step 8: Differential gene expression contrasts


## Cell Hashtag scRNAseq Analysis Track
#### Step 1: FASTQ pre-processing
All of the outputs of the CellRanger _counts_ pipeline are produced. For more information on the outputs, please visit the CellRanger [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count).

#### Step 2: Ambient RNA removal and create Seurat object
|Output type |Name|Description|
|:--|:--|:--|
|Figure|vioplot_sample_name.png | Sample-specific violin plot showing the distribution of cells according to QC metrics|
|Info|sample_name_RNA.txt|Sample-specific sparse matrix of RNA assay|
|Info|MetaData_sample_name.txt|Sample-specific dataframe showing the Seurat object metadata|
|Info|meta_info_sample_name.txt|Sample-specific text file showing the column names of the Seurat object metadata|
|Info|summary_sample_name.txt|Sample-specific text file showing the summary of QC metrics (Minimum, 1st Quartile, Median, Mean, 3rd Quartile, Maximum) | 
|Info|sessionInfo.txt|Session information for the R session|
|Data object|sample_name.rds|Sample-specific intermediate Seurat RDS object|


#### Step 3: Quality control and filtering
|Output type |Name|Description|
|:--|:--|:--|
|Figure|cellcycle_sample_name.png | Sample-specific violin plot showing the distribution of cells according to cell cycle S and G2M scores|
|Figure|dimplot_pca_sample_name.png | Sample-specific PCA showing the first two PCs|
|Figure|dimplot_umap_sample_name.png | Sample-specific UMAP|
|Figure|QC_vioplot_sample_name.png | Sample-specific violin plot showing the distribution of cells according to QC metrics after filtering|
|Figure|VariableFeaturePlot_sample_name.png | Sample-specific figure showing the most variably expressed genes|
|Info|sample_name_RNA.txt|Sample-specific sparse matrix of RNA assay|
|Info|MetaData_sample_name.txt|Sample-specific dataframe showing the Seurat object metadata|
|Info|meta_info_sample_name.txt|Sample-specific text file showing the column names of the Seurat object metadata|
|Info|most_variable_genes_sample_name.txt|Sample-specific text file showing the column names of the Seurat object metadata|
|Info|summary_sample_name.txt|Sample-specific text file showing the summary of QC metrics (Minimum, 1st Quartile, Median, Mean, 3rd Quartile, Maximum) | 
|Info|sessionInfo.txt|Session information for the R session|
|Data object|sample_name.rds|Sample-specific intermediate Seurat RDS object|

#### Step 4: Demultiplexing and doublet removal



#### Step 5: Integration and linear dimensional reduction
|Output type |Name|Description|
|:--|:--|:--|
|Figure|DimPlot_pca.png | PCA showing the first two PCs, colour-coded by sample|
|Figure|DimPlot_umap.png | UMAP, colour-coded by sample|
|Figure|elbow.png | Elbow plot to visualize the percentage of variance explained by each PC |
|Figure|Jackstraw_plot.png| Jackstraw plot to visualize the distribution of p-values for each PC|
|Info|seu_int_RNA.txt|Sparse matrix of integrated assay|
|Info|seu_int_MetaData.txt|Dataframe showing the Seurat object metadata|
|Info|meta_info_seu_step5.csv|Text file showing the column names of the Seurat object metadata|
|Info|sessionInfo.txt|Session information for the R session|
|Data object|seu_step5.rds|Integrated intermediate Seurat RDS object|

#### Step 6: Clustering
|Output type |Name|Description|
|:--|:--|:--|
|Figure|clustree_int.png  | Clustree plot showing the stability across the user-defied clustering resolutions|
|Figure|integrated_snn_res.png| UMAP at the user defined clustering-resolution|
|Figure|ARI.png | Mean and standard deviation of the Adjusted Rand Index (ARI) between clustering pairs at a user-defined resolution |
|Info|**test.xlsx**| Excel file showing the mean and standard deviation of the ARI between clustering pairs at a user-defined resolution |
|Info|seu_RNA.txt|Sparse matrix of integrated assay|
|Info|seu_MetaData.txt|Dataframe showing the Seurat object metadata|
|Info|meta_info.csv |Text file showing the column names of the Seurat object metadata|
|Info|sessionInfo.txt|Session information for the R session|
|Data object|seu_step6.rds|Intermediate Seurat RDS object|


#### Step 7: Cluster annotation
Cluster annotation method|Output type |Name|Description|
|:--|:--|:--|:--|
|General|Figure|umap.pdf | UMAP plot of integrated assay at the user-defined clustering resolution used for cluster annotation|
|General|Figure|umap_splitted.pdf | UMAP plot of integrated assay at the user-defined clustering resolution used for cluster annotation, split by sample|
|Cluster marker GSEA (Method 1)|Figure|heatmap.pdf | Heatmap showing the expression of the top marker genes across cells, stratified by cluster|
|Cluster marker GSEA (Method 1)|Figure|plotenrich.pdf | Barplot showing the 20 most enriched terms for a particular cluster and cell type library|
|Module score (Method 2)|Figure|module_score_gene_set.png | UMAP plot showing the module score across cells for user-defined gene sets|
|Reference-based annotations (Method 3)|Figure|UMAP_transferred_labels.pdf |UMAP plots showing the cluster annotations from the reference Seurat object projected onto the query Seurat object|
|Visualize select features|Figure|select_feature_dot_plot.pdf|Dotplot showing the expression of user-defined features at the cluster level|
|Visualize select features|Figure|select_feature_violin_plot.pdf |Violin plot showing the expression of user-defined features at the cluster level|
|Visualize select features|Figure|select_feature_feature_plot.pdf  |UMAP plots showing the expression of user-defined features at the cell level|
|General|Info|meta_info_seu_step7.txt |Text file showing the column names of the Seurat object metadata|
|General|Info|sessionInfo_marker.txt|Session information for the R session|
|Cluster marker GSEA (Method 1)|Info|cluster_just_genes.xlsx  |Excel file showing the marker genes for each cluster|
|Cluster marker GSEA (Method 1)|Info|cluster_whole.xlsx  |Excel file showing the marker genes and corresponding summary statistics for each cluster|
|Cluster marker GSEA (Method 1)|Info|ClusterMarkers.csv  |csv file showing the marker genes and corresponding summary statistics for each cluster|
|Cluster marker GSEA (Method 1)|Info| top_sel.csv |csv file showing the top n marker genes for each cluster. The user defined n in the execution parameters|
|Cluster marker GSEA (Method 1)|Info| Er.genes.csv|Enrichment terms and the corresponding statistics for a particular cluster and cell type library|
|Module score (Method 2)|Info| geneset_by_cluster.csv|Mean module score across clusters for each user-defined gene set|
|General|Data object|seu_step7.rds|Intermediate Seurat RDS object|
|Cluster marker GSEA (Method 1)|Data object| ClusterMarkers.rds |RDS object containing the marker genes for each cluster|

#### Step 8: Differential gene expression contrasts



