---
layout: post
title: Steps of scRNA pipeline to run Non-hashtag data
description: A short introduction to scRNA pipeline to run Non-hashtag data
date: 2023-06-16
author: Saeid Amiri
published: true
tags: scRNA 
categories: 
comments: false
---
## ScRNAbox pipeline: Standard scRNAseq Analysis Track
## Contents
- [Introduction](#introduction)
  - [Setup](#setup)
  - [Step 1: FASTQ pre-processing](#step-1-fastq-pre-processing)
  - [Step 2: Ambient RNA removal and create Seurat object](#step-2-ambient-rna-removal-and-create-seurat-object)  
  - [Step 3: Quality control and filtering](#step-3-quality-control-and-filtering)
  - [Step 4: Doublet removal](#step-4-doublet-removal)
  - [Step 5: Integration and linear dimensional reduction](#step-5-integration-and-linear-dimensional-reduction)
  - [Step 6: Clustering](#step-6-clustering)   
  - [step 7: Cluster annotation](#step-7-cluster-annotation)    
  - [step 8: Differential gene expression contrasts](#step-8-differential-gene-expression-contrasts)     

 - - - -

## Introduction 
This guide provides an overview of the Analytical Steps that comprise the Standard Analysis Track of the scRNAbox pipeline. The Standard Analysis Track is designed for scRNAseq experiments where each sample is captured and sequenced separately; thus, users should have unique FASTQ files for each of the samples in their experiment. If instead samples were labelled with sample-specific barcodes and pooled prior to sequencing, users should leverage the [Cell Hashtag scRNAseq](HTO.md) Analysis Track.<br /> 

The Analytical Steps involved in the Standard Analysis Track of the scRNAbox pipeline are outlined in the figure below.<br />  
 
 <img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/328966c3-04dc-4e42-91cb-c6181026d1d6" width="550" height="100">


**Note:** This tutorial assumes that `scrnabox.slurm`,`CellRanger`, `R`, and the required R packages have already been installed onto the HPC system. If this is not the case, please visit [Installation](installation.md) to do so before proceeding. If the required packages are installed, you can proceed to [Setup](#setup).

 - - - -

### Setup

Before running the pipeline, create a dedicated folder for the analysis (hereafter referred to as the working directory). Then, export the path to the working directory and the path to `scrnabox.slurm`:
```
mkdir ~/working_directory
cd ~/working_directory

export SCRNABOX_HOME=~/scrnabox.slurm
export SCRNABOX_PWD=~/working_directory
```

Next, run the pipeline initiation Step (`--steps 0`) and define the Standard scRNAseq Analysis Track (`--method SCRNA`) using the following command from the working directory:
```
cd ~/working_directory 

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 0 \
--method SCRNA
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

Upon completing the setup procedures, users can run their analysis using the scRNAbox pipeline. 

 - - - -

### scRNAbox Analytical Steps
Specific Analytical Steps are called using the `--steps` flag. The output of each Analytical Step is deposited into its respective folder within the working directory (e.g. `~/working_directory/step1`). Prior to running each Analytical Step, users are strongly encouraged to modify the execution parameters of the analysis using the adjustable, Step-specific parameters text files. The parameters text files are located in `~/working_directory/job_info/parameters`. To ensure replicability, a summary report file documents the execution parameters for each iteration of each Analytical Step, which is located in `~/working_directory/job_info/summary_report.txt`.


For detailed descriptions of each Analytical Step please see our pre-print manuscript. 
- - - -

### Step 1: FASTQ pre-processing
In this step, feature-barcode expression matrices are generated from FASTQ files using the CellRanger _counts_ pipeline. Prior to running CellRanger, `library.csv` files must be prepared to define the FASTQ files for each sample. ScRNAbox provides an option for automating this process or users may manually prepare the libraries. For more information, please see the the [Library preparation](library_prep.md) tutorial.<br />

The following parameters are adjustable for Step 1:

|Parameter|Default|Description|
|:--|:--|:--|
|REF_DIR_GRCH|NULL|Path to reference genome for FASTQ alignment. 10X Genomics reference genomes are available for download. For more information see their [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ct).|
|R1LENGTH|NULL|Minimum number of bases to retain for R1 sequence of gene expression|
|MEMPERCORE|30|For clusters whose job managers do not support memory requests, it is possible to request memory in the form of cores. This option will scale up the number of threads requested via the __MRO_THREADS__ variable according to how much memory a stage requires when given to the ratio of memory on your nodes.|
|par_automated_library_prep|No| Whether or not to perform automated library prep. Alternatively, you may set this parameter to "no" and manually prepare the libraries.|
|par_fastq_directory|NULL|Path to directory containing the FASTQ files. This directory should only contain FASTQ files for the experiment.|
|par_sample_names|NULL|The sample names used to name the FASTQ files according to CellRanger nomeclature|
|par_rename_samples|Yes| Whether or not you want to rename your samples. These names will be used to identify cells in the Seurat objects|
|par_new_sample_names|NULL| New sample names. Make sure they are defined in the same order as 'par_sample_names'|
|par_paired_end_seq|Yes| Whether or not paired-end sequencing was performed|

Given that CellRanger runs a user interface, it is recommended to run Step 1 in a 'screen'. To run Step 1, use the following command:
```
screen -S run_scrnabox
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 1
```

The resulting output files are deposited into `~/working_directory/step1`


**Note:** If you do not have access to FASTQ files for your experiment, you may intiate the pipeline at which ever Analytical Step takes your data object as input. In the case where FASTQ files are not available, users do not have to create the `samples_info` folder. For more information see [Processed Data](PROC.md).  
- - - -

### Step 2: Ambient RNA removal and create Seurat object 
In this step, the ambient RNA rate is estimated and the gene expression profiles are corrected for RNA contamination (optional) using SoupX (Young et al. 2020). Then, CellRanger- or SoupX-generated feature-barcode expression matrices are transformed into Seurat objects. Genes expressed in less than a minimum number of cells and cells expressing less than a minimum number of genes can be filtered.<br />
<br />
The following parameters are adjustable for Step 2:

|Parameter|Default|Description|
|:--|:--|:--|
|par_save_RNA| No| Whether or not to export an RNA expression matrix|
|par_save_metadata| No|Whether or not to export a metadata dataframe|
|par_ambient_RNA| Yes|Whether or not to correct the feature-barcode expression matrices for ambient RNA contamination|
|par_count_matrices| NULL|If users skipped Step 1, the may provide the path to a directory that contains existing feature-barcode expression matrices to initiate the pipeline at Step 2 |
|par_min.cells_L| 0|Only retain genes expressed in a minimum number of cells|
|par_min.features_L| 0|Only retain cells expressing a minimum number of genes|

To run Step 2, use the following command:
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 2
```

The resulting output files are deposited into `~/working_directory/step2`
 - - - -
### Step 3: Quality control and filtering
Low quality cells are filtered based on the user-defined thresholds for:

- the number of genes detected per cell;
- the number of unique transcripts detected per cell;
- the percentage of mitochondrial-encoded transcripts; 
- the percentage of ribosomal-encoded transcripts.

In addition, mitochondrial- and ribosomal-encoded genes can be filtered out, as well as a custom user-defined list of genes. Finally, normalization and scaling is performed on the individual Seurat objects prior to cell-cycle scoring. <br />
<br />
The following parameters are adjustable for Step 3:

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


To run Step 3, use the following command:
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 3 
```
The resulting output files are deposited into `~/working_directory/step3`

 - - - -
### Step 4: Doublet removal
Doublets are identified and removed from downstream analysis (optional) using the DoubletFinder tool (McGinnis et al. 2019).  <br />
<br />
The following parameters are adjustable for Step 4:

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

To run Step 4, use the following command:
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 4 
```
The resulting output files are deposited into `~/working_directory/step4`
 - - - -
### Step 5: Integration and linear dimensional reduction
Individual Seurat objects are integrated to enable the joint analysis across samples using Seurat's integration algorithm (Stuart et al. 2019); if experiments are limited to a single sequencing run, the integration Step can be bypassed. Normalization, scaling, and linear dimensional reduction is then performed on the resulting Seurat object to inform the optimal parameters for clustering in Step 6.<br />
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
The resulting output files are deposited into `~/working_directory/step5`

 - - - -
### Step 6: Clustering 
Clustering is performed to define groups of cells with similar expression profiles using the graph-based clustering approach implemented in the Seurat framework (Macosko et al. 2015).<br />
<br />
The following parameters are adjustable for Step 6:

|Parameter|Default|Description|
|:--|:--|:--|
|par_save_RNA| No| Whether or not to export an RNA expression matrix|
|par_save_metadata| No|Whether or not to export a metadata dataframe|
|par_skip_step5|No|Whether or not the user skipped integration in Step 5| 
|par_FindNeighbors_dims|30| Number of dimensions from linear dimensional reduction used as input to identify neighbours. Can be informed by the elbow and Jackstraw plots produced in Step 5|
|par_RunUMAP_dims|30| Number of dimensions from linear dimensional reduction used as input to compute UMAP. Can be informed by the elbow and Jackstraw plots produced in Step 5 |
|par_FindNeighbors_k.param|60|Defines k for the k-nearest neighbor algorithm|
|par_FindNeighbors_prune.SNN|1/15|Sets the cutoff for acceptable Jaccard index when computing the neighborhood overlap for the shared nearest-neighbour (SNN) construction
|par_FindClusters_resolution|0.1 to 0.9, in intervals of 0.1|Value of the clustering resolution parameter. You may provide multiple resolution values|
|par_compute_ARI|Yes| Whether or not you want to compute the Adjusted Rand Index (ARI) between clusters at a given clustering resolution|
|par_RI_reps|100|Number of iterations for clustering the data at a given resolution in order to calculate the ARI|

To run Step 6, use the following command:
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 6 
```
The resulting output files are deposited into `~/working_directory/step6`

 - - - -
### Step 7: Cluster annotation
Cell populations, or clusters, are annotated to define cell types by three distinct methods:<br />
<br />
 **Method 1. Cluster marker gene set enrichment analysis (GSEA)**: Seurat's _FindAllMarkers_ function is used to identify differentially expressed marker genes (DEG) by the Wilcoxon rank-sum test (Macosko et al. 2015). DEGs in the positive direction     (Log2 fold-change > 0.00) are then tested for enrichment across user-defined gene set libraries that define cell types using the EnrichR tool (Chen et al. 2013).<br />
 <br />
 **Method 2. Module score**: Seurat’s implementation (_AddModuleScore_) of Tirosh et al.’s algorithm is used to comparatively quantify the expression of gene sets across clusters at the single-cell level (Tirosh et al. 2016). Users must          define their desired gene sets in the parameters file of Step 7. Gene sets should be defined in a csv file, where the column names correspond to the arbitrary name of the gene set and the corresponding rows define the genes in the gene set. For example: <br />
 

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
 
<br /> 
 **Method 3. Reference-based annotation**: Seurat's _FindTransferAnchors_ and _TransferData_ functions are used to leverage cell-type annotations from a reference Seurat object and generate annotation predictions for the query dataset (Macosko et        al. 2015). User's must define the location of their referene Seurat object in the parameters file of Step 7.<br />

 In addition to the three annotation methods, scRNAbox allows users to **visualize the expression of select features** at the cluster or cell level via a user-provided list of gene identifiers in order to further inform the cellular species that make up specific clusters. 

The following parameters are adjustable for Step 7:

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


To identify cluster marker genes (**Method 1**), compute the module score for user-defined gene sets (**Method 2**), and **visualize select features**, use the following command:
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--marker T
```

Then, to test the top cluster marker genes for enrichment across gene set libraries that define cell types (**Method 1**), use the following command:
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--enrich T
```

**Note:** In order to test the cluster marker genes for enrichment across EnrichR libraries, the HPC must have access to the internet. If your HPC cannot access the internet, it is possible to run the enrichment step on your local machine directly in R. To do so, begin by downloading the Step 6 Seurat RDS object (`~/working_directory/step6/objs6/seu_step6.rds `)  and cluster marker RDS object (`~/working_directory/step7/info7/marker/ClusterMarkers.rds`) to your computer using the following commands:

```
scp username@beluga.computecanada.ca:~/working_directory/step7/info7/marker/ClusterMarkers.rds  ~/Desktop/working_directory
scp username@beluga.computecanada.ca:~/working_directory/step6/objs6/seu_step6.rds  ~/Desktop/working_directory

```
Then run the follwing code in R:

```
SCRNABOX_PWD <- "~/Desktop/working_directory"
dir.create("SCRNABOX_PWD/step7/annot")
level_cluster='integrated_snn_res.0.7'
ClusterMarkers='SCRNABOX_PWD/step7/info7/ClusterMarkers.rds'
PWD='SCRNABOX_PWD/step7/annot/'
PSUE='SCRNABOX_PWD/step6/objs6/seu_step6.rds'
top_sel=5
db <- c('Descartes_Cell_Types_and_Tissue_2021','CellMarker_Augmented_2021','Azimuth_Cell_Types_2021')
scrnaboxR::annotation(level_cluster,ClusterMarkers,PWD,PSUE,top_sel,db)
```

To perform reference-based annotation (**Method 3**), run the following command: 

```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--fta T
```
The resulting output files are deposited into `~/working_directory/step7`.

 - - - -
### Step 8: Differential gene expression contrasts
This step performs differential gene expression (DGE) analyses according to user-defined contrasts. Contrast can be performed between samples (e.g. case vs control; **sample-sample contrasts**) or between samples, stratified by cell type (e.g. case vs control for excitatory neurons only; **sample-cell contrasts**).
 
The following parameters are adjustable for Step 8:

|Parameter|Default|Description (the cluster annotation method associated with the parameter is shown)|
|:--|:--|:--|
|par_save_RNA| No| Whether or not to export an RNA expression matrix|
|par_save_metadata| No|Whether or not to export a metadata dataframe|
|par_level_cluster|integrated_snn_res.0.7| The cluster resolution that you used in the cluster annotation module. If you skipped integration in Step 5, use par_level_cluster='RNA_snn_res.0.7', if you want to proceed with a clustering resolution of 0.7|
|par_step8_clus_label|NULL|List of user-currated cluster labels obtained from the annotation module. Make sure to have the same number of labels as clusters at the desired clustering resolution.|
|par_new_genotype|yes|Whether or not you want to add new sample labels to simplify the contrasts. For example, you may wish to set both control1 and control2 as control.|
|par_old_sample_label|NULL|list of old sample labels (i.e. those used to create the samples_info folder in the setup procedures)|
|par_new_sample_label|NULL|list of new sample labels corresponding to the old sample labels defined in the parameter above|


The first step of the DGE contrasts module is to create a  DGEListobject. This step may require alot of RAM, we suggest 3*size(seu_int_clu.rds). Users can adjust the RAM in the `scrnabox_config.ini` file located in `~/working_directory/job_info/configs`. To create a DGEListobject, run the following command:

```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--dgelist T
```

The next step is to perform DGE contrasts. As mentioned, users can perform **sample-sample contrasts** or **sample-cell contrasts**.

To perform **sample-sample contrasts**, users must first modify the `step8_contrast_genotype.txt` contrast matrix located in `~/working_directory/job_info/parameters`. The contrast matrix comprises three columns that define the arbitrary name of the contrast (cont_name), the control sample (control), and the case sample (versus). For example:

```
cont_name control versus
design1 control case1,case2,case3,case4,case5
design2 control case1
design3 control case2
design4 control case3
design5 control case4
design6 control case5
design7 control case1,case2,case3,case4
design8 control case1,case2,case3
design9 control case1,case2
```
To perform **sample-sample contrasts**, run the following command:

```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--genotype T
```

To perform **sample-cell contrasts**, users must first modify the `step8_contrast_celltype.txt` contrast matrix located in `~/working_directory/job_info/parameters`. The contrast matrix comprises four columns that define the arbitrary name of the contrast (cont_name), the cell type (cell), the control sample (control), and the case sample (versus). For example:

```
cont_name cell control versus
design1 excitatory_neuron control case1,case2,case3,case4,case5
design2 excitatory_neuron control case1
design3 excitatory_neuron control case2
design4 excitatory_neuron control case3
design5 excitatory_neuron control case4
design6 excitatory_neuron control case5
design7 excitatory_neuron control case1,case2,case3,case4
design8 excitatory_neuron control case1,case2,case3
design9 excitatory_neuron control case1,case2
```
To perform **sample-cell contrasts**, run the following command:
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--celltype T
```
The resulting output files are deposited into `~/working_directory/step8`.

**Note:** If you have a large number of contrasts to run, it may be more efficient to split them up and submit batch jobs instead.



