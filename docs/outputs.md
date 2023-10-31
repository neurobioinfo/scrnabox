## Outputs of each step in the scRNAbox pipeline
- [Introduction](#introduction)
- [Outputs](#outputs)
    - [Step 1: FASTQ to expression matrix](#step-1-fastq-to-gene-expression-matrix)
    - [Step 2: Create Seurat object and remove ambient RNA](#step-2-create-seurat-object-and-remove-ambient-rna)
    - [Step 3: Quality control and generation of filtered data objects](#step-3-quality-control-and-generation-of-filtered-data-objects)
    - [Step 4: Doublet removal (standard track)](#step-4-doublet-removal-standard-track)
    - [Step 4: Demultiplexing and doublet detection (HTO track)](#step-4-demultiplexing-and-doublet-detection-hto-track)
    - [Step 5: Creation of a single Seurat object from all samples](#step-5-creation-of-a-single-seurat-object-from-all-samples)
    - [Step 6: Clustering](#step-6-clustering)
    - [Step 7: Cluster annotation](#step-7-cluster-annotation)
    - [Step 8: Differential gene expression](#step-8-differential-gene-expression)

- - - -

## Introduction
The outputs of each Step of the scRNAbox pipeline are deposited into a step-specific folder in the working directory which contains three sub folders:
```
working_directory
└──step1
    ├── figs1
    ├── info1
    └── objs1
```
 - The `figs/` folder contains figures;
 - The `info/` folder contains text files and tables;
 - The `objs/` folder contains intermediate Seurat RDS objects.

 **Note:** If users re-run an Analytical Step, the outputs from the previous run will automatically be overwritten. If you do not want to lose the outputs from a previous run, it is important to copy the materials to a separate directory. One exception to this is when annotating data in Step 7; users can re-run the **annotate** step as many times as they wish and each interation will add a new metadata column to the already existing Seurat object.
- - - -
## Outputs
#### Step 1: FASTQ to gene expression matrix 
All of the outputs of the CellRanger _counts_ pipeline are produced. For more information on the outputs, please visit the CellRanger [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count).
- - - -

#### Step 2: Create Seurat object and remove ambient RNA
|Output type |Name|Description|
|:--|:--|:--|
|Figure|ambient_RNA_estimation_sample_name.pdf| Sample-specific probability density plot showing the ambient RNA estimation. For more information see [here](https://cran.r-project.org/web/packages/SoupX/vignettes/pbmcTutorial.html)|
|Figure|ambient_RNA_markers_sample_name.pdf| Sample-specific figure showing the marker genes used for ambient RNA estimation. For more information see [here]https://cran.r-project.org/web/packages/SoupX/vignettes/pbmcTutorial.html)|
|Figure|vioplot_sample_name.pdf | Sample-specific violin plot showing the distribution of cells according to QC metrics|
|Figure|zoomed_in_vioplot_sample_name.pdf | Sample-specific violin plot showing the distribution of cells according to QC metrics. The minimum value to the mean is shown.|
|Figure|cell_cycle_dim_plot_sample_name.pdf | Sample-specific principal component analysis of cell-cycle genes, colour-coded by the cell cycle score of each cell.|
|Info|sample_name_ambient_rna_summary.rds| Sample-specific summary of ambient RNA estimation by SoupX|
|Info|sample_name_RNA.txt|Sample-specific sparse matrix of RNA assay|
|Info|estimated_ambient_RNA_sample_name.txt| Sample-specific ambient RNA estimation.|
|Info|MetaData_sample_name.txt|Sample-specific dataframe showing the Seurat object metadata|
|Info|meta_info_sample_name.txt|Sample-specific text file showing the column names of the Seurat object metadata|
|Info|summary_sample_name.txt|Sample-specific text file showing the summary of QC metrics (Minimum, 1st Quartile, Median, Mean, 3rd Quartile, Maximum) | 
|Info|sessionInfo.txt|Session information for the R session|
|Data object|sample_name.rds|Sample-specific intermediate Seurat RDS object|

- - - -

#### Step 3: Quality control and generation of filtered data objects
|Output type |Name|Description|
|:--|:--|:--|
|Figure|dimplot_pca_sample_name.pdf | Sample-specific PCA showing the first two PCs|
|Figure|elbow_sample_name.pdf | Elbow plot to visualize the percentage of variance explained by each PC |
|Figure|filtered_QC_vioplot_sample_name.pdf | Sample-specific violin plot showing the distribution of cells according to QC metrics after filtering|
|Figure|VariableFeaturePlot_sample_name.pdf | Sample-specific figure showing the most variably expressed genes|
|Info|sample_name_RNA.txt|Sample-specific sparse matrix of RNA assay|
|Info|MetaData_sample_name.txt|Sample-specific dataframe showing the Seurat object metadata|
|Info|meta_info_sample_name.txt|Sample-specific text file showing the column names of the Seurat object metadata|
|Info|most_variable_genes_sample_name.txt|Sample-specific text file showing the column names of the Seurat object metadata|
|Info|summary_sample_name.txt|Sample-specific text file showing the summary of QC metrics (Minimum, 1st Quartile, Median, Mean, 3rd Quartile, Maximum) | 
|Info|sessionInfo.txt|Session information for the R session|
|Data object|sample_name.rds|Sample-specific intermediate Seurat RDS object|

- - - -

#### Step 4: Doublet removal (standard track)
|Output type |Name|Description|
|:--|:--|:--|
|Figure|sample_nameDF.classifications.pdf | Sample-specific UMAP plot showing droplet classifications (singlet or doublet)|
|Figure|sample_doublet_summary.pdf | Sample-specific violin plot showing pANN value across singlet and doublet assignments; sample-specific bar plot showing the number of singlets and doublets.|
|Info|n_predicted_doublets_sample_name.txt|Sample-specific text file showing the number of identified doublets.|
|Info|sample_name_RNA.txt|Sample-specific sparse matrix of RNA assay|
|Info|MetaData_sample_name.txt|Sample-specific dataframe showing the Seurat object metadata|
|Info|meta_info_sample_name.txt|Sample-specific text file showing the column names of the Seurat object metadata|
|Info|sessionInfo.txt|Session information for the R session|
|Data object|sample_name.rds|Sample-specific intermediate Seurat RDS object|

- - - -

#### Step 4: Demultiplexing and doublet detection (HTO track)
|Output type |Name|Description|
|:--|:--|:--|
|Figure|run_name_DotPlot_HTO_MSD.pdf|Run-specific dot plot showing the enrichment of barcode-labels across cell assignments |
|Figure|run_name_Heatmap_HTO_MSD.pdf |Run-specific heatmap showing the enrichment of barcode-labels across cell assignments  |
|Figure|run_name_Ridgeplot_HTO_MSD.pdf |Run-specific ridge plot showing the enrichment of barcode-labels across cell assignments |
|Figure|run_name_HTO_dimplot_pca_.pdf | Run-specific PCA of antibody assay |
|Figure|run_name_HTO_dimplot_umap_.pdf | Run-specific UMAP of antibody assay|
|Figure|run_name_nCounts_RNA_MSD.pdf | Run-specific violin plot showing the number of unque transcripts across cell assignments|
|Info| run_name.rds_old_antibody_label_MULTIseqDemuxHTOcounts.csv| Run-specific list of sample-specific barcode labels used in the experiment|
|Info| run_name_MULTIseqDemuxHTOcounts.csv| Run-specific number of cells assigned to each sample|
|Info| run_namefiltered_MULTIseqDemuxHTOcounts.csv| Run-specific number of cells assigned to each sample after removal of doublet and negative droplets|
|Info| run_name_meta_info_.txt|Run-specific text file showing the column names of the Seurat object metadata|
|Info| run_name_MetaData.txt|Run-specific dataframe showing the Seurat object metadata|
|Info| run_name_RNA.txt|Run-specific sparse matrix of RNA assay|
|Info|sessionInfo.txt|Session information for the R session|
|Data object|run_name.rds|Run-specific intermediate Seurat RDS object|

- - - -

#### Step 5: Creation of a single Seurat object from all samples 
|Output type |Name|Description|
|:--|:--|:--|
|Figure|intergrated_DimPlot_pca.pdf | PCA showing the first two PCs of integrated assay, colour-coded by sample|
|Figure|integrated_DimPlot_umap.pdf | UMAP of integrated assay, colour-coded by sample|
|Figure|integrated_elbow.pdf | Elbow plot to visualize the percentage of variance explained by each PC for the integrated assay |
|Figure|integrated_Jackstraw.pdf |Jackstraw plot to visualize the distribution of p-values for each PC for the integrated assay |
|Figure|merge_DimPlot_pca.pdf | PCA showing the first two PCs of merged object, colour-coded by sample|
|Figure|merge_DimPlot_umap.pdf | UMAP of merged object, colour-coded by sample|
|Figure|merge_elbow.pdf | Elbow plot to visualize the percentage of variance explained by each PC for the merged object |
|Figure|merge_Jackstraw.pdf |Jackstraw plot to visualize the distribution of p-values for each PC for the merged object |
|Info|seu_int_RNA.txt|Sparse matrix of integrated assay|
|Info|seu_int_MetaData.txt|Dataframe showing the integrated object metadata |
|Info|integrated_meta_info_seu_step5.csv|Text file showing the column names of the integrated object metadata|
|Info|seu_merge_RNA.txt|Sparse matrix of merged data object|
|Info|seu_merge_MetaData.txt|Dataframe showing the merged object metadata |
|Info|merge_meta_info_seu_step5.csv|Text file showing the column names of the merged object metadata|
|Info|sessionInfo.txt|Session information for the R session|
|Data object|seu_step5.rds|Integrated intermediate Seurat RDS object|

- - - -

#### Step 6: Clustering
|Output type |Name|Description|
|:--|:--|:--|
|Figure|clustree_int.pdf | Clustree plot showing the stability across the user-defied clustering resolutions|
|Figure|integrated_snn_res.pdf| UMAP at the user defined clustering-resolution|
|Figure|ARI.pdf | Mean and standard deviation of the Adjusted Rand Index (ARI) between clustering pairs at a user-defined resolution |
|Info|clustering_ARI.xlsx| Excel file showing the mean and standard deviation of the ARI between clustering pairs at a user-defined resolution |
|Info|seu_RNA.txt|Sparse matrix of integrated assay|
|Info|seu_MetaData.txt|Dataframe showing the Seurat object metadata|
|Info|meta_info.csv |Text file showing the column names of the Seurat object metadata|
|Info|sessionInfo.txt|Session information for the R session|
|Data object|seu_step6.rds|Intermediate Seurat RDS object|

- - - -

#### Step 7: Cluster annotation
|Cluster annotation method|Output type |Name|Description|
|:--|:--|:--|:--|
|**Tool 1:** Cluster marker GSEA |Figure|heatmap.pdf | Heatmap showing the expression of the top marker genes across cells, stratified by cluster|
|**Tool 1:** Cluster marker GSEA |Figure|plotenrich.pdf | Barplot showing the 20 most enriched terms for a particular cluster and cell type library|
|**Tool 2:** Profile known markers|Figure|module_score_gene_set.pdf | UMAP plot showing the module score across cells for user-defined gene sets|
|**Tool 2:** Profile known markers|Figure|select_feature_dot_plot.pdf|Dotplot showing the expression of user-defined features at the cluster level|
|**Tool 2:** Profile known markers|Figure|select_feature_violin_plot.pdf |Violin plot showing the expression of user-defined features at the cluster level|
|**Tool 2:** Profile known markers|Figure|select_feature_feature_plot.pdf  |UMAP plots showing the expression of user-defined features at the cell level|
|**Tool 3:** Reference-based annotations|Figure|UMAP_transferred_labels.pdf |UMAP plots showing the cluster annotations from the reference Seurat object projected onto the query Seurat object|
|**Annotate**|Figure|clustering_name_cluster_annotation.pdf  |UMAP plot of the integrated assay showing the cluster annotation |
|**Annotate**|Figure|clustering_name_split_cluster_annotation.pdf |UMAP plot of the integrated assay showing the cluster annotation, split by sample|
|**General**|Info|meta_info_seu_step7.txt |Text file showing the column names of the Seurat object metadata|
|**General**|Info|sessionInfo_marker.txt|Session information for the R session|
|**Tool 1:** Cluster marker GSEA |Info|cluster_just_genes.xlsx  |Excel file showing the marker genes for each cluster|
|**Tool 1:** Cluster marker GSEA |Info|cluster_whole.xlsx  |Excel file showing the marker genes and corresponding summary statistics for each cluster|
|**Tool 1:** Cluster marker GSEA |Info|ClusterMarkers.csv  |csv file showing the marker genes and corresponding summary statistics for each cluster|
|**Tool 1:** Cluster marker GSEA |Info| top_sel.csv |csv file showing the top n marker genes for each cluster. The user defined n in the execution parameters|
|**Tool 1:** Cluster marker GSEA |Info| Er.genes.csv|Enrichment terms and the corresponding statistics for a particular cluster and cell type library|
|**Tool 1:** Cluster marker GSEA |Data object| ClusterMarkers.rds |RDS object containing the marker genes for each cluster|
|**Tool 2:** Module score |Info| geneset_by_cluster.csv|Mean module score across clusters for each user-defined gene set|
|**Tool 3:** Reference-based annotations|Info| reference_predictions_summary.xlsx |Number of cells from each cluster assigned a particular annotation based of the reference|
|**General**|Data object|seu_step7.rds|Intermediate Seurat RDS object|

- - - -

#### Step 8: Differential gene expression
|DGE method|Output type |Name|Description|
|:--|:--|:--|:--|
|**general**|Figure|contrast_name.pdf | Volcano plot of showing differentially expressed genes|
|**Cell-based DGE**|Info|contrast_name_DEG.csv| Differentially exppresed genes identified for the user-defined contrast|
|**Sample-based DGE**|Info|Aggregated_expression_summary.csv | Aggregated counts across user-defined sample groups|
|**Sample-based DGE**|Info|SampleBased_DGEsummarytable.csv | Number of differentially expressed genes in the positive and negative direction for each user-defined contrast|
|**Sample-based DGE**|Info|DGE_contrast_name.csv| Differentially exppresed genes identified for the user-defined contrast|
|**General**|Info|seu_RNA.txt|Sparse matrix of integrated assay|
|**General**|Info|seu_MetaData.txt|Dataframe showing the Seurat object metadata|
|**General**|Info|meta_info.csv |Text file showing the column names of the Seurat object metadata|
|**General**|Info|sessionInfo.txt|Session information for the R session|
|**General**|Data object|seu_step8.rds|Intermediate Seurat RDS object|

