---
layout: post
title:  A Guide to Analyzing HTO with the Scrnabox Pipeline
description: A short introduction to  Hashtag oligonucleotide analyzing using scrnabox pipeline
date: 2023-06-16
author: Saeid Amiri
published: true
tags: scRNA HTO
categories: 
comments: false
---
## ScRNAbox pipeline: Cell Hashtag scRNAseq Analysis Track

## Contents
- [Introduction](#introduction)
  - [Setup](#setup)
  - [Step 1: FASTQ to gene expression matrix](#step-1-fastq-to-gene-expression-matrix)
  - [Step 2: Create Seurat object and remove ambient RNA ](#step-2-create-seurat-object-and-remove-ambient-rna)   
  - [Step 3: Quality control and filtering](#step-3-quality-control-and-filtering)
  - [Step 4: Demultiplexing and doublet removal](#step-4-demultiplexing-and-doublet-removal)
  - [Step 5: Integration and linear dimensional reduction](#step-5-integration-and-linear-dimensional-reduction)
  - [Step 6: Clustering](#step-6-clustering)   
  - [step 7: Cluster annotation](#step-7-cluster-annotation)    
  - [step 8: Differential gene expression contrasts](#step-8-differential-gene-expression-contrasts)     


## Introduction 
This guide provides an overview of the Analytical Steps that comprise the Cell Hashtag Analysis Track of the scRNAbox pipeline. The Cell Hashtag Analysis Track is designed for multiplexed scRNAseq experiments, whereby samples are tagged with sample-specific barcodes, pooled, and sequenced together; thus, users should have FASTQ files that contain scRNAseq data from multiple samples. If instead samples were sequenced separately, users should leverage the [Standard scRNAseq](SCRNA.md) Analysis Track.<br /> 


The Analytical Steps involved in the Cell Hashtag Analysis Track of the scRNAbox pipeline are outlined in the figure below.<br />  

 <p align="center">
 <img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/4356d563-4cee-42fd-ac26-502f9c020209" width="550" height="100">
 </p>

**Note:** This tutorial assumes that `scrnabox.slurm`,`CellRanger`, `R`, and the required R packages have already been installed onto the HPC system. If this is not the case, please visit [Installation](installation.md) to do so before proceeding. If the required packages are installed, you can proceed to [Setup](#setup).

 - - - -

### Setup

Before running the pipeline, create a dedicated folder for the analysis (hereafter referred to as the working directory). Then, define the path of the working directory and the path to `scrnabox.slurm`:
```
mkdir ~/working_directory
cd ~/working_directory

export SCRNABOX_HOME=~/scrnabox.slurm
export SCRNABOX_PWD=~/working_directory
```
               
Next, run the pipeline initiation Step (`--steps 0`) and define the Standard scRNAseq Analysis Track (`--method HTO`) using the following command from the working directory:
```
cd ~/working_directory 

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 0 \
--method HTO
```

After running the pipeline initiation Step, the structure of the working directory should be as follows:
```
 working_directory
 └── job_info
     ├── configs
     ├── logs
     └── parameters

```
- The `configs/` directory contains the `scrnabox_config.ini` file which allows users to specify their job allocations (memory, threads, and walltime) for each Analytical Step using the Slurm Workload Manager; <br /> 
- The `logs/` directory records the events of each Analytical Step; <br />
- The `parameters/` directory contains adjustable, Step-specific text files which allow users to define the execution parameters for each Analytical Step. <br />

Next, navigate to the `scrnabox_config.ini` file in `~/working_directory/job_info/configs` to define the path to the R library, the version of R, and the path to CellRanger. For example: 

```
MODULECELLRANGER=mugqic/cellranger/5.0.1
R_VERSION=4.2.1
R_LIB_PATH=Path/to/R/library
```
**Note:** For more information, please see the [Job cofigurations](config.md) sections of the scRNAbox documentation.

Upon completing the setup procedures, users can run their analysis using the scRNAbox pipeline.

 - - - -

### scRNAbox Analytical Steps
Specific Analytical Steps are called using the `--steps` flag. There are three componets that correspond to each Analytical Step in the scRNAbox pipeline:

1) **Job configurations**; <br />
2) **Execution parameters**; <br />
3) **Outputs**.

Prior to running each Analytical Step, users should modify their **job configurations** using the `scrnabox_config.ini` located in `~/working_directory/job_info/configs`.  Similarly, users should modify the **execution parameters** prior to each Analytical Step using the parameters text files located in `~/working_directory/job_info/parameters`. The **outputs** of each Analytical Step are deposited into its respective folder within the working directory (e.g. `~/working_directory/step1`).

**Note:** For more information, please see the [Job cofigurations](config.md), [Execution parameters](reference.md) and [Outputs](outputs.md) sections of the scRNAbox documentation. For a detailed description of each Analytical Step please see our [pre-print manuscript](). 

- - - -

### Step 1: FASTQ to gene expression matrix
In this step, feature-barcode expression matrices are generated from FASTQ files using the CellRanger _counts_ pipeline. Prior to running CellRanger, `library.csv` and `feature_ref.csv` files must be prepared to define the FASTQ files and sample-specific barcodes for each sequencing run. ScRNAbox provides an option for automating this process or users may manually prepare the libraries. For more information, please see the the [CellRanger library preparation](library_prep.md) tutorial.<br />

The following parameters are adjustable for Step 1 (`~/working_directory/job_info/parameters/step1_par.txt`):

|Parameter|Default|Description|
|:--|:--|:--|
|par_automated_library_prep|Yes|Whether or not to perform automated library prep. Alternatively, you may set this parameter to "no" and manually prepare the libraries.|
|par_fastq_directory|NULL|Path to directory containing the FASTQ files. This directory should only contain FASTQ files for the experiment.|
|par_RNA_run_names|NULL|The names of the sequencing runs for the RNA assay|
|par_HTO_run_names|NULL|The names of the sequencing runs for the HTO assay|
|par_seq_run_names|NULL|The user-selected name for the sequencing run.  These names will be used to identify cells in the Seurat objects|
|par_paired_end_seq|Yes|Whether or not paired-end sequencing was performed|
|id|NULL|Barcode ID which will be used to track the feature counts|
|name|NULL|The user-selected name for the barcode identifier|
|read|R2|Which RNA sequencing read contains the barcode sequence. This value Will be either R1 or R2.|
|pattern|NULL|The pattern of the barcode identifiers|
|sequence|NULL|The nucleotide sequence associated with the barcode identifier|
|REF_DIR_GRCH|NULL|Path to reference genome for FASTQ alignment. 10X Genomics reference genomes are available for download. For more information see their [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ct).|
|R1LENGTH|NULL|Minimum number of bases to retain for R1 sequence of gene expression|
|MEMPERCORE|30|For clusters whose job managers do not support memory requests, it is possible to request memory in the form of cores. This option will scale up the number of threads requested via the __MRO_THREADS__ variable according to how much memory a stage requires when given to the ratio of memory on your nodes.|


**Note:** The execution parameters for each analystical step can be adjusted in the Step-specific text files located in `~/working_directory/job_info/parameters/`

Given that CellRanger runs a user interface and is not submitted as a Job, it is recommended to run Step 1 in a **'screen'** which will allow the the task to keep running if the connection is broken. To run Step 1, use the following command:
```
screen -S run_scrnabox
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 1
```

The resulting output files are deposited into `~/working_directory/step1`. The expression matrix, features, and barcode files outputed by CellRanger are located in `~/working_directory/step1/sample/ouput_folder/outs/raw_feature_bc_matrix`.

**Note:** If you do not have access to FASTQ files for your experiment, you may intiate the pipeline at which ever Analytical Step takes your data object as input. In the case where FASTQ files are not available, users do not have to create the `samples_info` folder. For more information see [Processed Data](PROC.md).  

- - - -


### Step 2: Create Seurat object and remove ambient RNA 
In this step, the CellRanger outputs generated in Step 1 (expression matrix, features, and barcodes) are used to create a Seurat object for each sequencing run. The ambient RNA quantity is estimated and there is an option to correct gene expression profiles for RNA contamination using SoupX ([Young et al. 2020](https://academic.oup.com/gigascience/article/9/12/giaa151/6049831)). Then, CellRanger (if not removing ambient RNA) or SoupX (if removing ambient RNA) feature-barcode expression matrices are transformed into Seurat objects.

The following parameters are adjustable for Step 2 (`~/working_directory/job_info/parameters/step2_par.txt`):

|Parameter|Default|Description|
|:--|:--|:--|
|Save_RNA| No| Whether or not to export an RNA expression matrix|
|Save_metadata| No|Whether or not to export a metadata dataframe|
|par_ambient_RNA| Yes|Whether or not to correct the feature-barcode expression matrices for ambient RNA contamination|
|par_count_matrices| NULL|If users skipped Step 1, the may provide the path to a directory that contains existing feature-barcode expression matrices to initiate the pipeline at Step 2 |
|par_min.cells_L| 0|Only retain genes expressed in a minimum number of cells|
|par_normalization.method|LogNormalize|Method to use for normalization|
|par_scale.factor|10000|Scale factor for scaling the data|
|par_selection.method|vst|Method for choosing the top variable features|
|par_nfeatures|2500|Number of features to select as top variable features|

To run Step 2, use the following command:
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 2
```

The resulting output files are deposited into `~/working_directory/step2`. For a description of the outputs see [here](outputs.md).

- - - -


### Step 3: Quality control and filtering
Low quality cells are filtered based on the user-defined thresholds for:

- the number of genes detected per cell;
- the number of unique transcripts detected per cell;
- the percentage of mitochondrial-encoded transcripts; 
- the percentage of ribosomal-encoded transcripts.

In addition, mitochondrial- and ribosomal-encoded genes can be filtered out, as well as a custom user-defined list of genes. Cell cycle genes can be regressed. Finally, normalization and scaling is performed on the individual Seurat objects prior to cell-cycle scoring. <br />
<br />
The following parameters are adjustable for Step 3 (`~/working_directory/job_info/parameters/step3_par.txt`):

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
|par_remove_mitochondrial_genes|Yes| Whether or not to remove mitochondrial genes|
|par_remove_ribosomal_genes|Yes| Whether or not to remove ribosomal genes|
|par_remove_genes|NULL|If users want to remove specific genes from their data, they may define a list of gene identifiers|
|par_regress_cell_cycle_genes|Yes|Whether or not to regress cell cycle genes|
|par_normalization.method|LogNormalize|Method to use for normalization|
|par_scale.factor|10000|Scale factor for scaling the data|
|par_selection.method|vst|Method for choosing the top variable features|
|par_nfeatures|2500|Number of features to select as top variable features|
|par_top|10|Number of most variable features to be reported in the csv file|
|par_npcs_pca|30|Total Number of principal components to compute and store for principal component analysis (PCA)|


To run Step 3, use the following command:
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 3 
```
The resulting output files are deposited into `~/working_directory/step3`. For a description of the outputs see [here](outputs.md).

- - - -

### Step 4: Demultiplexing and doublet removal
Seurat’s implementation (_MULTIseqDemux_) of the tag assignment algorithm outlined in Multi-seq is used to demultiplex pooled samples and identify doublets according to the expression matrices of the sample-specific barcodes ([McGinnis et al 2019](https://pubmed.ncbi.nlm.nih.gov/31209384/)).

The following parameters are adjustable for Step 4 (`~/working_directory/job_info/parameters/step4_par.txt`):

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

To perform demultiplexing and doublet detection, the first step is to obtain the barcode labels used in the analysis by running the following command:

```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 4 \
--msd T 
```
 The names of the existing barcode labels can be revised to be more descriptive in the execution parameters of this Step (par_old_antibody_label; par_new_antibody_label)

 Next, perform demultiplexing and doublet detection by running the following command:

```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 4 
```

The resulting output files are deposited into `~/working_directory/step4`. For a description of the outputs see [here](outputs.md).
- - - -

### Step 5: Integration and linear dimensional reduction
Individual Seurat objects are merged and integrated to enable the joint analysis across sample runs using Seurat's integration algorithm ([Stuart et al. 2019](https://pubmed.ncbi.nlm.nih.gov/31178118/)); if experiments are limited to a single sequencing run, the integration Step can be bypassed. However, Step 5 must be run because normalization, scaling, and linear dimensional reduction is then performed on the resulting Seurat object to inform the optimal parameters for clustering in Step 6.<br />
<br />
The following parameters are adjustable for Step 5:

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

To run Step 5, use the following command:
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 5 
```
The resulting output files are deposited into `~/working_directory/step5`. For a description of the outputs see [here](outputs.md).
- - - -
### Step 6: Clustering 
Clustering is performed to define groups of cells with similar expression profiles using the graph-based clustering approach implemented in the Seurat framework ([Butler et al. 2015](https://www.nature.com/articles/nbt.4096)).<br />
<br />
The following parameters are adjustable for Step 6:

|Parameter|Default|Description|
|:--|:--|:--|
|par_save_RNA| No| Whether or not to export an RNA expression matrix|
|par_save_metadata| No|Whether or not to export a metadata dataframe|
|par_skip_integration|No|Whether or not the user skipped integration in Step 5| 
|par_FindNeighbors_dims|30| Number of dimensions from linear dimensional reduction used as input to identify neighbours. Can be informed by the elbow and Jackstraw plots produced in Step 5|
|par_RunUMAP_dims|10| Number of dimensions to use as input features for uniform manifold approximation and projection (UMAP)|
|par_FindNeighbors_k.param|60|Defines k for the k-nearest neighbor algorithm|
|par_FindNeighbors_prune.SNN|1/15|Sets the cutoff for acceptable Jaccard index when computing the neighborhood overlap for the shared nearest-neighbour (SNN) construction
|par_FindClusters_resolution|0, 0.05, 0.2, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0|Value of the clustering resolution parameter. You may provide multiple resolution values|
|par_compute_ARI|Yes| Whether or not you want to compute the Adjusted Rand Index (ARI) between clusters at a given clustering resolution|
|par_RI_reps|100|Number of iterations for clustering the data at a given resolution in order to calculate the ARI|

To run Step 6, use the following command:
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 6 
```
The resulting output files are deposited into `~/working_directory/step6`. For a description of the outputs see [here](outputs.md).
- - - -
- - - -
### Step 7: Cluster annotation
ScRNAbox provides three Methods to assist in cluster annotation or identification of cell types:<br />
<br />
 **Method 1. Cluster marker gene set enrichment analysis (GSEA)**: Seurat's _FindAllMarkers_ function is used to identify differentially expressed marker genes (DEG) by the Wilcoxon rank-sum test ([Macosko et al. 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4481139/)). DEGs in the positive direction     (Log2 fold-change > 0.00) are then tested for enrichment across user-defined gene set libraries that define cell types using the EnrichR tool ([Chen et al. 2013](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-128)). To find marker genes use `par_run_find_marker = "yes"` and to perform GSEA use `par_run_enrichR = "yes"` in the execution parameters of Step 7. <br />
 <br />
 **Method 2. Module score**: Seurat’s implementation (_AddModuleScore_) of Tirosh et al.’s algorithm is used to comparatively quantify the expression of gene sets across clusters at the single-cell level ([Tirosh et al. 2016](https://pubmed.ncbi.nlm.nih.gov/27124452/). Users must          define their desired gene sets in the parameters file of Step 7. Gene sets should be defined in a csv file, where the column names correspond to the arbitrary name of the gene set and the corresponding rows define the genes in the gene set. For example: <br />
 

|gene_set_1|gene_set_2|gene_set_3|
|:--|:--|:--|
|CLDN5|IFITM3|TM4SF1|
|ANGPT2|EGFL7|IGFBP|
|FLT1|A2M|GIMAP7|
|DEPP1|SOD2|EMP2|
|TMSB10|PODXL|B2M|
|HLA-E|VWF|BST2|
|SLCO4A1|PECAM1|BSG|
|TGM2|SLC2A3|PARP14|
|IFI27|TSC22D1|NFKBIA|
|MT2A|HLA-B|ID3|
|EPAS1||
|IFITM2||
 
To compute module score use `par_run_module_score = "yes"` in the execution parameters of Step 7.


 **Method 3. Reference-based annotation**: Seurat's _FindTransferAnchors_ and _TransferData_ functions are used to leverage cell-type annotations from a reference Seurat object and generate annotation predictions for the query dataset ([Macosko et al. 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4481139/)). User's must define the location of their referene Seurat object in the parameters file of Step 7. To perform reference-based annotation use `par_run_reference = "yes"` in the execution parameters of Step 7.
<br />


 - - - -

In addition to the three Methods described above, users can **Visualize the expression of select features** and **Annotate** their Seurat object.

 **Visualize the expression of select features:** The expression of known markger genes can be visualized at the cluster or cell level via a user-provided list of gene identifiers in order to further inform the cell types that make up specific clusters. If users want to explore the expression of multiple gene lists they can do so using a csv file, where the column names correspond to the arbitrary name of the gene list and the corresponding rows define the genes in the list. For example:

|list_1|list_2|list_3|
|:--|:--|:--|
|CLDN5|IFITM3|TM4SF1|
|ANGPT2|EGFL7|IGFBP|
|FLT1|A2M|GIMAP7|
|DEPP1|SOD2|EMP2|
|TMSB10|PODXL|B2M|
|HLA-E|VWF|BST2|
|SLCO4A1|PECAM1|BSG|
|TGM2|SLC2A3|PARP14|
|IFI27|TSC22D1|NFKBIA|
|MT2A|HLA-B|ID3|
|EPAS1||
|IFITM2||

To visualize the expression of select features use `par_run_visualize_markers = "yes"` in the execution parameters of Step 7.

**Annotate:** Users can annotate their Seurat object and visualize their annotations by defining the clustering resolution that they wish to annotate and supplying a list of annotations corresponding to each cluster at the givenn resolution. The Annotation step can be performed multiple times as the previous annotations will not be removed from the Seurat object. To annotate you Seurat object use `par_run_annotate = "yes"` in the execution parameters of Step 7.

The following parameters are adjustable for Step 7:

|Annotation Method|Parameter|Default|Description|
|:--|:--|:--|:--|
|**Method 1**|par_run_find_marker|Yes|Whether or not to find marker genes for each cluster|
|**Method 1**|par_run_enrichR|Yes|Whether or not to run gene set enrichment analysis (GSEA) on the marker genes for each cluster using the EnrichR tools. Note that the HPC must have access to the internet to run GSEA.|
|**Method 2**|par_run_module_score|Yes|Whether or not to compute module score |
|**Method 3**|par_run_reference|Yes|Whether or not to perform reference-based annotation|
|**Visualize features**|par_run_visualize_markers|Yes| Whether or not to Visualize select features|
|**Annotate**|par_run_annotate|Yes|Whether or not to Annotate|
|**General**|par_save_RNA| No| Whether or not to export an RNA expression matrix|
|**General**|par_save_metadata| No|Whether or not to export a metadata dataframe|
|**General**|par_level_cluster| integrated_snn_res.0.7| The cluster resolution that you want to use for downstream analyses. If you skipped integration in Step 5, use par_level_cluster='RNA_snn_res.0.7', if you want to proceed with a clustering resolution of 0.7|
|**Method 1**|par_top_sel|5|Number of top markers to identify based on avg_log2FC|
|**Method 1**|par_db|Descartes_Cell_Types_and_Tissue_2021,<br /> CellMarker_Augmented_2021,<br />Azimuth_Cell_Types_2021|Character vector of EnrichR databases that define cell types. The top marker genes for each cluster will be tested for enrichment across these databases.|
|**Method 2**|par_module_score|NULL|Path defining the location of the directory that contains the csv file of the gene sets used to compute the module score|
|**Method 3**|par_reference|NULL| Path defining the location of the reference Seurat object|
|**Method 3**|par_reference_name|Reference| An arbitrary name for the reference object. This will be used to name the metadata slot.|
|**Method 3**|par_level_celltype|NULL|The name of the metadata column in the reference Seurat object that defines cell types|
|**Method 3**|par_FindTransferAnchors_dim|10| Number of dimensions from linear dimensional reduction used to find transfer anchors between the reference and query Seurat objects|
|**Method 3**|par_futureglobalsmaxSize|50000 * 1024^2|This will increase your RAM usage so set this number mindfully|
|**Visualize features**|par_select_features_list|NULL| A list of features to visualize|
|**Visualize features**|par_select_features_csv|NULL|If you want to define multiple lists of features to visualize, you can do so with a csv file. The header should contain the list names and all features belonging to the same list should be in the same column.|
|**Annotate**|par_annotate_resolution|NULL| Which clustering resolution you want to annotate|
|**Annotate**|par_name_metadata|clustering_label_1| The name of the metadata slot that will contain the annotations|
|**Annotate**|par_annotate_labels|NULL| A list of cluster labels. There must as many labels as clusters at the defined clustering resolution. Please refrain from using "_" when annotating.|


To run Step 7, use the following command:
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--marker T
```

The resulting output files are deposited into `~/working_directory/step7`. For a description of the outputs see [here](outputs.md). 

 - - - -
#### Marker GSEA without access to the internet
In order to test the cluster marker genes for enrichment across EnrichR libraries (`par_run_enrichR = "yes"`), the HPC must have access to the internet. If your HPC cannot access the internet, it is possible to run the enrichment step on your local machine directly in R. To do so, begin by downloading the `ClusterMarkers.csv` file obtained from running `par_run_find_marker = "yes"` to your computer using the following commands:

```
scp username@beluga.computecanada.ca:~/working_directory/step7/info7/marker/ClusterMarkers.csv ~/Desktop/working_directory
scp username@beluga.computecanada.ca:~/working_directory/step6/objs6/ClusterMarkers.csv  ~/Desktop/working_directory

```
Then run the follwing code in R:

```
install.packages("devtools", dependencies = TRUE)
library(devtools)
devtools::install_github("neurobioinfo/scrnabox/scrnaboxR")
library(scrnaboxR)

PWD="~/Desktop/working_directory"
cluster_marker="PWD/ClusterMarkers.csv"
db=c("Descartes_Cell_Types_and_Tissue_2021","CellMarker_Augmented_2021","Azimuth_Cell_Types_2021")

scrnaboxR::annotation(PWD,cluster_marker,db)

```


 - - - -
### Step 8: Differential gene expression contrasts
This step performs differential gene expression (DGE) analyses according to user-defined contrasts. Contrast can be performed between samples (e.g. case vs control; **sample-sample contrasts**) or between samples, stratified by cell type (e.g. case vs control for excitatory neurons only; **sample-cell contrasts**). In addition, **psuedo-bulk analysis** can be performed. 
 
The following parameters are adjustable for Step 8:

|Parameter|Default|Description (the cluster annotation method associated with the parameter is shown)|
|:--|:--|:--|
|par_run_add_metadata|Yes| Whether or not to add metadata to the Seurat object to facilitate differential gene expression contrasts.|
|par_run_sample_sample_wilcoxon|Yes| Whether or to perform DGE contrasts between samples across all cells using the Wilcoxon method.|
|par_run_sample_cell_wilcoxon|Yes| Whether or to perform DGE contrasts between samples stratified by cell type using the Wilcoxon method.|
|par_run_pseudo_bulk|Yes| Whether or not to perform pseudo-bulk analysis|
|par_save_RNA| No| Whether or not to export an RNA expression matrix|
|par_save_metadata| No|Whether or not to export a metadata dataframe|
|par_metadata|NULL|Path to a csv file defining new metadata that should be added to the Seurat object to facilitate DEG analysis. At least one column should contain "**orig.ident**".|

In addition to the execution parameters, users should fill in the contrast matrices to define their contrasts. <br />
 - - - -
### Sample-sample contrast matrix
To perform sample-sample contrasts, users must fill in the `step8_contrast_genotype.txt` file located in `~/working_directory/job_info/parameters`. The sample-sample contrasts matrix contains the following columns:

1. **contast_name:** An abritrary name for the contrast
2. **meta_data_variable:** The metadata slot containing the Sample IDs defined in group1 and group2 
3. **group1:** A list of sample IDs to be contrasted against the sample IDs listed in group2
4. **group2:**A list of sample IDs to be contrasted against the sample IDs listed in group1

Multiple contrasts can be defined in the same file. In addition, multiple samples can be listed under group1 and group 2. For example: 
```
contrast_name meta_data_variable group1 group2
design1 orig.ident Control1,Control2,Control3 Case1,Case2,Case3
design2 orig.ident Control1 Case1,Case2,Case3
design3 DiseaseStatus HC Case
```

- - - -
### Sample-cell contrast matrix
To perform sample-cell contrasts, users must fill in the `step8_contrast_celltype.txt` file located in `~/working_directory/job_info/parameters`. The sample-cell contrasts matrix contains the following columns:

1. **contast_name:** An abritrary name for the contrast
2. **meta_data_celltype:** The metadata slot containing the cell type annotations
3. **cell_type:** The cell type used for differential gene expression
4. **meta_data_variable:** The metadata slot containing the Sample IDs defined in group1 and group2 
5. **group1:** A list of sample IDs to be contrasted against the sample IDs listed in group2
6. **group2:**A list of sample IDs to be contrasted against the sample IDs listed in group1

Multiple contrasts can be defined in the same file. In addition, multiple samples can be listed under group1 and group 2. For example: 

```
contrast_name meta_data_celltype cell_type meta_data_variable group1 group2
design1 clustering_1 Oligodendrocytes orig.ident Control1,Control2,Control3 Case1,Case2,Case3
design1 clustering_1 Oligodendrocytes orig.ident Control1 Case1,Case2,Case3
design2 clustering_2 Microglia DiseaseStatus HC PD
```
- - - -

### Pseudo-bulk contrast matrix
To perform pseudo-bulk contrasts, users must fill in the `step8_contrast_pseudo_bulk.txt` file located in `~/working_directory/job_info/parameters`. The default pseudo-bulk contrasts matrix contains the following columns:

1. **ContrastName:** An abritrary name for the contrast
2. **CellType:** The metadata slot containing the cell type annotations. Pseudo-bulk DGE analysis will be performed on all cell types defined.
3. **MainContrast:** The metadata slot defining the main variables for the contrast (e.g. Case or Control)
2. **SampleID:** The metadata slot containing the sample IDs.

Only one contrast can can be defined in the same file. Pseudo-bulk analysis will not work without >1 Sample for each group defined in the MainContrast. For example: 

```
ContrastName CellType MainContrast SampleID
Pseudo_design1 clustering_1 DiseaseStatus orig.ident
```
- - - -

To run Step 8, use the following command:
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--dgelist T
```

The resulting output files are deposited into `~/working_directory/step8`. For a description of the outputs see [here](outputs.md). 

