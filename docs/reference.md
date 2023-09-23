## Adjustable execution parameters for the scRNAbox pipeline

- [Standard scRNAseq Analysis Track](#standard-scrnaseq-analysis-track)
- [Cell Hashtag scRNAseq Analysis Track](#cell-hashtag-scrnaseq-analysis-track)

## Standard scRNAseq Analysis Track
#### Step 1: FASTQ pre-processing
|Parameter|Default|Description|
|:--|:--|:--|
|REF_DIR_GRCH|NULL|Path to reference genome for FASTQ alignment. 10X Genomics reference genomes are available for download. For more information see their [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ct).|
|R1LENGTH|NULL|Minimum number of bases to retain for R1 sequence of gene expression|
|MEMPERCORE|30|For clusters whose job managers do not support memory requests, it is possible to request memory in the form of cores. This option will scale up the number of threads requested via the __MRO_THREADS__ variable according to how much memory a stage requires when given to the ratio of memory on your nodes.|

#### Step 2: Ambient RNA removal and create Seurat object
|Parameter|Default|Description|
|:--|:--|:--|
|par_save_RNA| No| Whether or not to export an RNA expression matrix|
|par_save_metadata| No|Whether or not to export a metadata dataframe|
|par_ambient_RNA| Yes|Whether or not to correct the feature-barcode expression matrices for ambient RNA contamination|
|par_count_matrices| NULL|If users skipped Step 1, the may provide the path to a directory that contains existing feature-barcode expression matrices to initiate the pipeline at Step 2 |
|par_min.cells_L| 0|Only retain genes expressed in a minimum number of cells|
|par_min.features_L| 0|Only retain cells expressing a minimum number of genes|

#### Step 3: Quality control and filtering
|Parameter|Default|Description|
|:--|:--|:--|
|par_save_RNA| No| Whether or not to export an RNA expression matrix|
|par_save_metadata| No|Whether or not to export a metadata dataframe|
|par_seurat_object| NULL |If users skipped Steps 1 and 2, the may provide the path to a directory that contains existing Seurat objects to initiate the pipeline at Step 3|
|par_nFeature_RNA_L|NULL |Only retain cells expressing a minimum number of genes|
|par_nFeature_RNA_U|NULL |Only retain cells expressing a maximum number of genes|
|par_nCount_RNA_L|NULL |Only retain cells with a minimum number of unique transcripts|
|par_nCount_RNA_U|NULL |Only retain cells with a maximum number of unique transcripts|
|par_mitochondria_percent_L|NULL | Only retain cells with a minimum percentage of mitochondrial genes|
|par_mitochondria_percent_U|NULL |Only retain cells with a maximum percentage of mitochondrial genes|
|par_ribosomal_percent_L|NULL |Only retain cells with a minimum percentage of ribosomal genes|
|par_ribosomal_percent_U|NULL |Only retain cells with a maximum percentage of ribosomal genes|
|par_log10GenesPerUMI_L|NULL | Only retain cells with a minimum number of genes per unique molecular identifier|
|par_log10GenesPerUMI_U|NULL | Only retain cells with a maximum number of genes per unique molecular identifier|
|par_remove_mitochondrial_genes|Yes| Whether or not to remove mitochondrial genes|
|par_remove_ribosomal_genes|Yes| Whether or not to remove ribosomal genes|
|par_remove_genes|NULL|If users want to remove specific genes from their data, they may define a list of gene identifiers|
|par_normalization.method|LogNormalize|Method to use for normalization|
|par_scale.factor|10000|Scale factor for scaling the data|
|par_selection.method|vst|Method for choosing the top variable features|
|par_nfeatures|2500|Number of features to select as top variable features|
|par_top|10|Number of most variable features to be reported in the csv file|
|par_npcs_pca|30|Total Number of principal components to compute and store for principal component analysis (PCA)|
|par_cells|500|Number of cells to include in Seurat's _dimheatmap_ function|
|par_dims|12|Number of dimensions to include in Seurat's _dimheatmap_ function|
|par_dims_umap|10|Number of dimensions to use as input features for uniform manifold approximation and projection (UMAP)|
|par_n.neighbors|65|Number of neighboring points to use in local approximations of manifold structure|

#### Step 4: Doublet removal
|Parameter|Default|Description|
|:--|:--|:--|
|par_save_RNA| No| Whether or not to export an RNA expression matrix|
|par_save_metadata| No|Whether or not to export a metadata dataframe|
|par_dropDN| Yes| Whether or not to remove predicted doublets from downstream analyses|
|par_PCs|20| The number of statistically significant principal components. Can be informed by elbow plot produced in Step 3|
|par_pN|0.25| The number of artificial doublets to generate. DoubletFinderr is largely invariant to this parameter. We suggest keeping 0.25|
|par_sct|FALSE|Logical representing whether SCTransform was used during original Seurat object pre-processing|
|par_sample_names|NULL| A list of sample names for each sample in the experiement, corresponding to the expected doublet rates listed in the parameter below. Sample names should be the same as those used to produce the `samples_info` folder during the setup procedures.|
|par_expected_doublet_rate|NULL| A list of expected doublet rates for each sample, corresponding to the sample names listed in the above parameter|

#### Step 5: Integration and linear dimensional reduction
|Parameter|Default|Description|
|:--|:--|:--|
|par_save_RNA| No| Whether or not to export an RNA expression matrix|
|par_save_metadata| No|Whether or not to export a metadata dataframe|
|par_skip_integration| No| Whether or not to skip integration. This is applicable for experiments that comprises of only one sequencing run.|
|par_FindIntegrationAnchors_dim|25|Which dimensions to use from the canonical correlation analysis (CCA) to specify the neighbor search space|
|par_DefaultAssay|RNA|The assay to perform normalization, scaling, and linear dimensiona reduction on. For most use cases this will be RNA.|
|par_normalization.method|LogNormalize|Method to use for normalization|
|par_selection.method|vst|Method for detecting top variable features|
|par_nfeatures|2500|Number of features to select as top variable features|
|par_RunUMAP_n.neighbors|65|Number of neighboring points used in local approximations of manifold structure|
|par_RunPCA_npcs|30| Total Number of principal components to compute and store for principal component analysis (PCA)
|par_RunUMAP_dims|10| Number of dimensions to use as input features for uniform manifold approximation and projection (UMAP)|
|par_compute_jackstraw |No|Whether or not to perform JackStraw computation. This computation takes a long time.|

#### Step 6: Clustering
|Parameter|Default|Description|
|:--|:--|:--|
|par_save_RNA| No| Whether or not to export an RNA expression matrix|
|par_save_metadata| No|Whether or not to export a metadata dataframe|
|par_skip_step5|No|Whether or not the user skipped integration in Step 5| 
|par_FindNeighbors_dims|30| Number of dimensions from linear dimensional reduction used as input to identify neighbours. Can be informed by the elbow and Jackstraw plots produced in Step 5|
|par_FindNeighbors_k.param|60|Defines k for the k-nearest neighbor algorithm|
|par_FindNeighbors_prune.SNN|1/15|Sets the cutoff for acceptable Jaccard index when computing the neighborhood overlap for the shared nearest-neighbour (SNN) construction
|par_FindClusters_resolution|0.1 to 0.9, in intervals of 0.1|Value of the clustering resolution parameter. You may provide multiple resolution values|
|par_compute_ARI|Yes| Whether or not you want to compute the Adjusted Rand Index (ARI) between clusters at a given clustering resolution|
|par_RI_reps|100|Number of iterations for clustering the data at a given resolution in order to calculate the ARI|

#### Step 7: Cluster annotation
|Parameter|Default|Description (the cluster annotation method associated with the parameter is shown)|
|:--|:--|:--|
|par_save_RNA| No| Whether or not to export an RNA expression matrix|
|par_save_metadata| No|Whether or not to export a metadata dataframe|
|par_level_cluster| integrated_snn_res.0.7| The cluster resolution that you want to use for downstream analyses. If you skipped integration in Step 5, use par_level_cluster='RNA_snn_res.0.7', if you want to proceed with a clustering resolution of 0.7|
|par_top_sel|5|**Method 1:** Number of top markers to identify based on avg_log2FC|
|par_db|Descartes_Cell_Types_and_Tissue_2021,<br /> CellMarker_Augmented_2021,<br />Azimuth_Cell_Types_2021|**Method 1:** Character vector of EnrichR databases that define cell types. The top marker genes for each cluster will be tested for enrichment across these databases.|
|par_compute_module_score|Yes|**Method 2:** Whether or not to perform the module score computation|
|par_compute_module_score|Yes|**Method 2:** Whether or not to perform the module score computation|
|par_module_score|NULL|**Method 2:** Path defining the location of the directory that contains the csv file of the gene sets used to compute the module score|
|par_reference|NULL|**Method 3:** Path defining the location of the reference Seurat object|
|par_level_celltype|NULL|**Method 3:** The name of the metadata column in the reference Seurat object that defines cell types|
|par_FindTransferAnchors_dim|10|**Method 3:** Number of dimensions from linear dimensional reduction used to find transfer anchors between the reference and query Seurat objects|
|par_futureglobalsmaxSize|50000 * 1024^2|**Method 3:** This will increase your RAM usage so set this number mindfully|
|par_visualize_select_features|No|**Visualize select features:** Whether or not to visualize select features|
|par_select_features|NULL|**Visualize select features:** list of gene identifiers to visualize the expression of select features|

#### Step 8: Differential gene expression contrasts
|Parameter|Default|Description (the cluster annotation method associated with the parameter is shown)|
|:--|:--|:--|
|par_save_RNA| No| Whether or not to export an RNA expression matrix|
|par_save_metadata| No|Whether or not to export a metadata dataframe|
|par_level_cluster|integrated_snn_res.0.7| The cluster resolution that you used in the cluster annotation module. If you skipped integration in Step 5, use par_level_cluster='RNA_snn_res.0.7', if you want to proceed with a clustering resolution of 0.7|
|par_step8_clus_label|NULL|List of user-currated cluster labels obtained from the annotation module. Make sure to have the same number of labels as clusters at the desired clustering resolution.|
|par_new_genotype|yes|Whether or not you want to add new sample labels to simplify the contrasts. For example, you may wish to set both control1 and control2 as control.|
|par_old_sample_label|NULL|list of old sample labels (i.e. those used to create the samples_info folder in the setup procedures)|
|par_new_sample_label|NULL|list of new sample labels corresponding to the old sample labels defined in the parameter above|

## Cell Hashtag scRNAseq Analysis Track
#### Step 1: FASTQ pre-processing
|Parameter|Default|Description|
|:--|:--|:--|
|REF_DIR_GRCH|NULL|Path to reference genome for FASTQ alignment. 10X Genomics reference genomes are available for download. For more information see their [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ct).|
|R1LENGTH|NULL|Minimum number of bases to retain for R1 sequence of gene expression|
|MEMPERCORE|30|For clusters whose job managers do not support memory requests, it is possible to request memory in the form of cores. This option will scale up the number of threads requested via the __MRO_THREADS__ variable according to how much memory a stage requires when given to the ratio of memory on your nodes.|

#### Step 2: Ambient RNA removal and create Seurat object
|Parameter|Default|Description|
|:--|:--|:--|
|Save_RNA| No| Whether or not to export an RNA expression matrix|
|Save_metadata| No|Whether or not to export a metadata dataframe|
|count_matrices| NULL|If users skipped Step 1, the may provide the path to a directory that contains existing feature-barcode expression matrices to initiate the pipeline at Step 2 |
|min.cells_L| 0|Only retain genes expressed in a minimum number of cells|
|min.features_L| 0|Only retain cells expressing a minimum number of genes|

#### Step 3: Quality control and filtering
|Parameter|Default|Description|
|:--|:--|:--|
|par_save_RNA| No| Whether or not to export an RNA expression matrix|
|par_save_metadata| No|Whether or not to export a metadata dataframe|
|par_seurat_object| NULL |If users skipped Steps 1 and 2, the may provide the path to a directory that contains existing Seurat objects to initiate the pipeline at Step 3|
|par_nFeature_RNA_L|NULL |Only retain cells expressing a minimum number of genes|
|par_nFeature_RNA_U|NULL |Only retain cells expressing a maximum number of genes|
|par_nCount_RNA_L|NULL |Only retain cells with a minimum number of unique transcripts|
|par_nCount_RNA_U|NULL |Only retain cells with a maximum number of unique transcripts|
|par_mitochondria_percent_L|NULL | Only retain cells with a minimum percentage of mitochondrial genes|
|par_mitochondria_percent_U|NULL |Only retain cells with a maximum percentage of mitochondrial genes|
|par_ribosomal_percent_L|NULL |Only retain cells with a minimum percentage of ribosomal genes|
|par_ribosomal_percent_U|NULL |Only retain cells with a maximum percentage of ribosomal genes|
|par_log10GenesPerUMI_L|NULL | Only retain cells with a minimum number of genes per unique molecular identifier|
|par_log10GenesPerUMI_U|NULL | Only retain cells with a maximum number of genes per unique molecular identifier|
|par_remove_mitochondrial_genes|Yes| Whether or not to remove mitochondrial genes|
|par_remove_ribosomal_genes|Yes| Whether or not to remove ribosomal genes|
|par_remove_genes|NULL|If users want to remove specific genes from their data, they may define a list of gene identifiers|
|par_normalization.method|LogNormalize|Method to use for normalization|
|par_scale.factor|10000|Scale factor for scaling the data|
|par_selection.method|vst|Method for choosing the top variable features|
|par_nfeatures|2500|Number of features to select as top variable features|
|par_top|10|Number of most variable features to be reported in the csv file|
|par_npcs_pca|30|Total Number of principal components to compute and store for principal component analysis (PCA)|
|par_cells|500|Number of cells to include in Seurat's _dimheatmap_ function|
|par_dims|12|Number of dimensions to include in Seurat's _dimheatmap_ function|
|par_dims_umap|10|Number of dimensions to use as input features for uniform manifold approximation and projection (UMAP)|
|par_n.neighbors|65|Number of neighboring points to use in local approximations of manifold structure|

#### Step 4: Doublet removal
|Parameter|Default|Description|
|:--|:--|:--|
|par_save_RNA| No| Whether or not to export an RNA expression matrix|
|par_save_metadata| No|Whether or not to export a metadata dataframe|
|par_normalization.method|CLR|Method for normalizing the HTO assay|
|par_scale.factor| 1000|Scale factor for scaling the HTO assay|
|par_selection.method|vst|Method for selecting the most variable features in the HTO assay|
|par_nfeatures|5|Number of features to select as top variable features for the HTO assay. This value is dependent on the number of sample specific barcodes used in the experiment|
|par_dims_umap|5|Number of dimensions to use as input features for uniform manifold approximation and projection (UMAP) of HTO assay|
|par_n.neighbor|65|Number of neighboring points to use in local approximations of manifold structure|
|par_dimensionality_reduction|Yes|Whether or not to perform linear dimensionality reduction on the HTO assay|
|par_npcs_pca|30|Total Number of principal components to compute and store for principal component analysis (PCA) of HTO assay|
|par_dropDN|Yes|Whether or not to remove predicted doublets and negatives from downstream analyses|
|par_label_dropDN|Doublet, Negative| Labels used to identify doublet and negative droplets|
|par_quantile|0.9|The quantile to use for droplet classification using _MULTIseqDemux_|
|par_autoThresh|TRUE| Whether or not to perform automated threshold finding to define the best quantile for droplet classification using _MULTIseqDemux_|
|par_maxiter|5|Maximum number of iterations to use if autoThresh = TRUE|
|par_RidgePlot_ncol|3|Number of columns used to display RidgePlots, which visualizes the enrichment of barcode labels across samples|
|par_old_antibody_label|NULL| If you wish to rename the barcode labels, first list the existing barcode labels in this parameter. old antibody labels can be identified in the "_old_antibody_label_MULTIseqDemuxHTOcounts" file produced by running Step 4 msd|
|par_new_antibody_label|NULL|If you wish to rename the barcode labels, list the new labels corresponding to the old labels listed in the parameter above|

#### Step 5: Integration and linear dimensional reduction
|Parameter|Default|Description|
|:--|:--|:--|
|par_save_RNA| No| Whether or not to export an RNA expression matrix|
|par_save_metadata| No|Whether or not to export a metadata dataframe|
|par_skip_integration| No| Whether or not to skip integration. This is applicable for experiments that comprises of only one sequencing run.|
|par_FindIntegrationAnchors_dim|25|Which dimensions to use from the canonical correlation analysis (CCA) to specify the neighbor search space|
|par_DefaultAssay|RNA|The assay to perform normalization, scaling, and linear dimensiona reduction on. For most use cases this will be RNA.|
|par_normalization.method|LogNormalize|Method to use for normalization|
|par_scale.factor|1000|Scale factor for scaling the data|
|par_selection.method|vst|Method for detecting top variable features|
|par_nfeatures|2500|Number of features to select as top variable features|
|par_RunUMAP_n.neighbors|65|Number of neighboring points used in local approximations of manifold structure|
|par_RunPCA_npcs|30| Total Number of principal components to compute and store for principal component analysis (PCA)
|par_RunUMAP_dims|10| Number of dimensions to use as input features for uniform manifold approximation and projection (UMAP)|
|par_compute_jackstraw |No|Whether or not to perform JackStraw computation. This computation takes a long time.|

#### Step 6: Clustering
|Parameter|Default|Description|
|:--|:--|:--|
|par_save_RNA| No| Whether or not to export an RNA expression matrix|
|par_save_metadata| No|Whether or not to export a metadata dataframe|
|par_skip_step5|No|Whether or not the user skipped integration in Step 5| 
|par_FindNeighbors_dims|30| Number of dimensions from linear dimensional reduction used as input to identify neighbours. Can be informed by the elbow and Jackstraw plots produced in Step 5|
|par_FindNeighbors_k.param|60|Defines k for the k-nearest neighbor algorithm|
|par_FindNeighbors_prune.SNN|1/15|Sets the cutoff for acceptable Jaccard index when computing the neighborhood overlap for the shared nearest-neighbour (SNN) construction
|par_FindClusters_resolution|0.1 to 0.9, in intervals of 0.1|Value of the clustering resolution parameter. You may provide multiple resolution values|
|par_compute_ARI|Yes| Whether or not you want to compute the Adjusted Rand Index (ARI) between clusters at a given clustering resolution|
|par_RI_reps|100|Number of iterations for clustering the data at a given resolution in order to calculate the ARI|

#### Step 7: Cluster annotation
|Parameter|Default|Description (the cluster annotation method associated with the parameter is shown)|
|:--|:--|:--|
|par_save_RNA| No| Whether or not to export an RNA expression matrix|
|par_save_metadata| No|Whether or not to export a metadata dataframe|
|par_level_cluster| integrated_snn_res.0.7| The cluster resolution that you want to use for downstream analyses. If you skipped integration in Step 5, use par_level_cluster='RNA_snn_res.0.7', if you want to proceed with a clustering resolution of 0.7|
|par_level_genotype|MULTI_ID_Lables|Name of the metadata column in your Seurat object that describes the sample names|
|par_top_sel|5|**Method 1:** Number of top markers to identify based on avg_log2FC|
|par_db|Descartes_Cell_Types_and_Tissue_2021,<br /> CellMarker_Augmented_2021,<br />Azimuth_Cell_Types_2021|**Method 1:** Character vector of EnrichR databases that define cell types. The top marker genes for each cluster will be tested for enrichment across these databases.|
|par_compute_module_score|Yes|**Method 2:** Whether or not to perform the module score computation|
|par_compute_module_score|Yes|**Method 2:** Whether or not to perform the module score computation|
|par_module_score|NULL|**Method 2:** Path defining the location of the directory that contains the csv file of the gene sets used to compute the module score|
|par_reference|NULL|**Method 3:** Path defining the location of the reference Seurat object|
|par_level_celltype|NULL|**Method 3:** The name of the metadata column in the reference Seurat object that defines cell types|
|par_FindTransferAnchors_dim|10|**Method 3:** Number of dimensions from linear dimensional reduction used to find transfer anchors between the reference and query Seurat objects|
|par_futureglobalsmaxSize|50000 * 1024^2|**Method 3:** This will increase your RAM usage so set this number mindfully|
|par_visualize_select_features|No|**Visualize select features:** Whether or not to visualize select features|
|par_select_features|NULL|**Visualize select features:** list of gene identifiers to visualize the expression of select features|

#### Step 8: Differential gene expression contrasts
|Parameter|Default|Description (the cluster annotation method associated with the parameter is shown)|
|:--|:--|:--|
|par_save_RNA| No| Whether or not to export an RNA expression matrix|
|par_save_metadata| No|Whether or not to export a metadata dataframe|
|par_level_cluster|integrated_snn_res.0.7| The cluster resolution that you used in the cluster annotation module. If you skipped integration in Step 5, use par_level_cluster='RNA_snn_res.0.7', if you want to proceed with a clustering resolution of 0.7|
|par_step8_clus_label|NULL|List of user-currated cluster labels obtained from the annotation module. Make sure to have the same number of labels as clusters at the desired clustering resolution.|
|par_new_genotype|yes|Whether or not you want to add new sample labels to simplify the contrasts. For example, you may wish to set both control1 and control2 as control.|
|par_old_antibody_label|NULL|list of old sample labels (i.e. those used to create the samples_info folder in the setup procedures)|
|par_new_antibody_label|NULL|list of new sample labels corresponding to the old sample labels defined in the parameter above|
