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
  - [Step 1: FASTQ pre-processing](#step-1-fastq-pre-processing)
  - [Step 2: Create Seurat object](#step-2-create-seurat-object)  
  - [Step 3: Quality control and filtering](#step-3-quality-control-and-filtering)
  - [Step 4: Demultiplexing and doublet removal](#step-4-demultiplexing-and-doublet-removal)
  - [Step 5: Integration and linear dimensional reduction](#step-5-integration-and-linear-dimensional-reduction)
  - [Step 6: Clustering](#step-6-clustering)   
  - [step 7: Cluster annotation](#step-7-cluster-annotation)    
  - [step 8: Differential gene expression contrasts](#step-8-differential-gene-expression-contrasts)     
- [Integrating Seurat objects](#integrating-seurat-objects)  


## Introduction 
This guide provides instructions for analyzing single-cell RNA sequencing (scRNAseq) data using the Cell Hashtag Analysis Track of the scRNAbox pipeline. The Cell Hashtag Analysis Track is designed for multiplexed scRNAseq experiments, whereby samples are tagged with sample-specific barcodes, pooled, and sequenced together; thus, users should have FASTQ files that contain scRNAseq data from multiple samples. If instead samples were sequenced separately, resulting in unique FASTQ files for each sample, users should leverage the [Standard scRNAseq](SCRNA.md) Analysis Track.<br /> 


The main component of the scRNAbox pipeline is `scrnabox.slurm`, which is an open-source pipeline for scRNAseq analysis that is specifically designed to run on high-performance computing (HPC) systems using the [Slurm Workload Manager](https://slurm.schedmd.com/). `scrnabox.slurm` outlines the Analytical Steps involved in a comprehensive scRNAseq analysis workflow, including FASTQ pre-processing, quality control and filtering, clustering, cluster annotation, and differential gene expression contrasts. The Analytical Steps involved in the Cell Hashtag Analysis Track of the scRNAbox pipeline are outlined in the figure below.<br />  

 <img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/eccddd8e-4ea2-4c1e-9427-8ba40e6418ba" width="550" height="100">

**Note:** This tutorial assumes that `scrnabox.slurm`,`cellranger`, `R`, and the required R packages have already been installed onto the HPC system. If this is not the case, please visit [Installation](installation.md) to do so before proceeding. If the required packages are installed, you can proceed to [Setup](#setup).

### Setup

Before running the pipeline, create a dedicated folder for the analysis (hereafter referred to as the working directory). Then, define the path of the working directory (`SCRNABOX_PWD=`) and the path to `scrnabox.slurm` (`SCRNABOX_HOME=`). For example:
```
mkdir ~/working_directory
cd ~/working_directory 
export SCRNABOX_HOME=~/scrnabox.slurm
export SCRNABOX_PWD=~/working_directory
```
For a description of the options for running `scrnabox.slurm` and to ensure that the path was properly defined, run the following command:

```
bash $SCRNABOX_HOME/launch_scrnabox.sh -h 
```

If the pipeline has been installed properly, the above command should return the folllowing:
```
        mandatory arguments:
                -d  (--dir)  = Working directory (where all the outputs will be printed) (give full path)
                --steps  =  Specify what steps, e.g., 2 to run just step 2, 2-4, run steps 2 through 4)

        optional arguments:
                -h  (--help)  = See helps regarding the pipeline options. 
                --method  = Choose what scRNA method you want to use; use HTO  and SCRNA for for hashtag nad Standard scRNA, respectively. 
                --nFeature_RNA_L  = Lower threshold of number of unique RNA transcripts for each cell, it filters nFeature_RNA > nFeature_RNA_L.  
                --nFeature_RNA_U  = Upper threshold of number of unique RNA transcripts for each cell, it filters --nFeature_RNA_U.  
                --nCount_RNA_L  = Lower threshold for nCount_RNA, it filters nCount_RNA > nCount_RNA_L   
                --nCount_RNA_U  = Upper threshold for  nCount_RNA, it filters nCount_RNA < nCount_RNA_U  
                --mitochondria_percent_L  = Lower threshold for the amount of mitochondrial transcript, it is in percent, mitochondria_percent > mitochondria_percent_L. 
                --mitochondria_percent_U  = Upper threshold for the amount of mitochondrial transcript, it is in percent, mitochondria_percent < mitochondria_percent_U. 
                --log10GenesPerUMI_U  = Upper threshold for the log number of genes per UMI for each cell, it is in percent,log10GenesPerUMI=log10(nFeature_RNA)/log10(nCount_RNA). mitochondria_percent < log10GenesPerUMI_U. 
                --log10GenesPerUMI_L  = Lower threshold for the log number of genes per UMI for each cell, log10GenesPerUMI=log10(nFeature_RNA)/log10(nCount_RNA). mitochondria_percent > log10GenesPerUMI_L.  
                --msd  = you can get the hashtag labels by running the following code 
                --marker  = Find marker. 
                --sinfo  = Do you need sample info? 
                --fta  = FindTransferAnchors 
                --enrich  = Annotation 
                --dgelist  = creates a DGEListobject from a table of counts obtained from seurate objects. 
                --genotype  = Run the genotype contrast. 
                --celltype  = Run the Genotype-cell contrast. 
                --cont  = You can directly call the contrast to the pipeline.  
                --seulist                = You can directly call the list of seurat objects to the pipeline. 
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
├── working_directory
    ├── job_info
        ├── configs
        ├── logs
        ├── parameters
```
- The `configs/` directory contains the `scrnabox_config.ini` file which allows users to specify their job allocations (memory, threads, and walltime) for each Analytical Step using the Slurm Workload Manager; <br /> 
- The `logs/` directory records the events of each Analytical Step; <br />
- The `parameters/` directory contains adjustable, Step-specific text files which allow users to define the execution parameters for each Analytical Step. <br />

Next, navigate to the `scrnabox_config.ini` file in `~/working_directory/job_info/configs` to define the path to the R library (`R_LIB_PATH=`), the version of R (`R_VERSION=`), and the path to CellRanger (`MODULECELLRANGER=`). For example: 

```
MODULECELLRANGER=mugqic/cellranger/5.0.1
R_VERSION=4.2.1
R_LIB_PATH=Path_to_R_library
```

Finally, in preparation for Step 1 (FASTQ pre-processing with CellRanger) users must create `library.csv` and `feature_ref.csv` files for each of their sequencing runs.<br />

#### library.csv
The `library.csv` file defines the necessary information of the FASTQ files for the experiment, including the gene expression and antibody assays. The structure of the `library.csv` file should be: <br />
```
fastqs,sample,library_type
~/fastqs/,CTRL1_GEX,Gene Expression
~/fastqs/,CTRL1_HTO,Antibody Capture
```
- The `fastqs` column defines the path to the directory that contains the FASTQ files for the experiment. <br /> 
- The `sample` column defines the sample name of the corresponding FASTQ file. Please note that FASTQ files must be named according to standard CellRanger nomenclature. For example, "CTRL1_S1_L001_R1_001.fastq". For more information please visit CellRanger's [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/fastq-input). <br />
- The `library_type` column defines the assay type. For the Cell Hashtag Analysis track, each sequencing run should have a "Gene Expression" and "Antibody Capture" assay. For more information, please visit CellRanger's [documentation]("https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis") <br />

For example, if the experiment comprises three sequencing runs the following steps should be taken: <br />

1) Navigate to the working directory and create a `samples_info` folder: <br />
```
cd ~/working_directory
mkdir samples_info
```
2) Navigate to the `samples_info` folder and create a folder for each sequencing run: <br />
```
cd samples_info
mkdir run1
mkdir run2
mkdir run3
```
3) Navigate to the folder for each sequencing and create the `library.csv` file. <br />

After performing steps 1-3 above, the structure of the samples_info folder for an experiment with three sequencing runs should be:
```
├── working_directory
    ├── samples_info
        ├── run1
            ├── library.csv
        ├── run2
            ├── library.csv
        ├── run3
            ├── library.csv
```
#### feature_ref.csv
The `feature_ref.csv` file defines the necessary information for processing the sample-specific barcodes that will eventually be used to demultiplex the pooled samples. For example, if there are four samples pooled together with four unique barcode identifiers, the structure of the `feature_ref.csv` file should be:
```
id,name,read,pattern,sequence,feature_type
Hash1,B0251_TotalSeqB,R2,5PNNNNNNNNNN(BC),GTCAACTCTTTAGCG,Antibody Capture
Hash2,B0252_TotalSeqB,R2,5PNNNNNNNNNN(BC),TGATGGCCTATTGGG,Antibody Capture
Hash3,B0253_TotalSeqB,R2,5PNNNNNNNNNN(BC),TTCCGCCTCTCTTTG,Antibody Capture
Hash4,B0254_TotalSeqB,R2,5PNNNNNNNNNN(BC),AGTAAGTTCAGCGTA,Antibody Capture
```
- The `id` column defines the barcode ID which will be used to track the feature counts. <br /> 
- The `name` column defines the arbitrary name for the barcode identifier. <br /> 
- The `read` column defines which RNA sequencing read contains the barcode sequence. This value Will be either R1 or R2.<br /> 
- The `pattern` column defines the pattern of the barcode identifiers. For more information please visit the 10X Genomics [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis#pattern)<br /> 
- The `sequence` column defines nucleotide sequence associated with the barcode identifier.<br /> 
- The `feature_type` column defines the type of feature used for sample identification. Please ensure that the feature_type in the `feature_ref.csv` file matches a library_type in the `library.csv` file.  <br /> 

For more information regarding the preparation of the `feature_ref.csv`, please see CellRanger's [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis).

`feature_ref.csv` files can be prepared the same way as the `library.csv` files. After producing the `feature_ref.csv` for each sequncing run, the structure of the samples_info folder for an experiment with three sequencing runs should be:
```
├── working_directory
    ├── samples_info
        ├── run1
            ├── library.csv
            ├── feature_ref.csv
        ├── run2
            ├── library.csv
            ├── feature_ref.csv
        ├── run3
            ├── library.csv
            ├── feature_ref.csv
```

Upon completing the above steps to setup the pipeline, users can run their analysis. Please note that if you do not have access to FASTQ files for your experiment, you may intiate the pipeline at which ever Analytical Step takes your data object as input. In the case where FASTQ files are not available, users do not have to create the `samples_info` folder. For more information see [Processed Data](PROC.md). 

### scRNAbox Analytical Steps
Specific Analytical Steps are called using the `--steps` flag. The output of each Analytical Step is deposited into its respective folder within the working directory (e.g. `~/working_directory/step1`). Prior to running each Analytical Step, users are strongly encouraged to modify the execution parameters of the analysis using the adjustable, Step-specific parameters text files. The parameters text files are located in `~/working_directory/job_info/parameters`. For instructions on how to modify these text files please see [FAQ](FAQ.md). To ensure replicability, a summary report file documents the execution parameters for each iteration of each analytical Step, which is located in `~/working_directory/job_info/summary_report.txt`.

For detailed descriptions of each Analytical Step please see our pre-print manuscript. 

### Step 1: FASTQ pre-processing
In this step, feature-barcode expression matrices are generated from FASTQ files using the CellRanger _counts_ pipeline. Given that CellRanger runs a user interface, it is recommended to run Step 1 in a 'screen'. <br />
<br />
The following parameters are adjustable for Step 1:

|Parameter|Default|Description|
|:--|:--|:--|
|REF_DIR_GRCH|NULL|Path to reference genome for FASTQ alignment. 10X Genomics reference genomes are available for download. For more information see their [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ct).|
|R1LENGTH|NULL|Minimum number of bases to retain for R1 sequence of gene expression|
|MEMPERCORE|30|For clusters whose job managers do not support memory requests, it is possible to request memory in the form of cores. This option will scale up the number of threads requested via the __MRO_THREADS__ variable according to how much memory a stage requires when given to the ratio of memory on your nodes.|

To run Step 1, use the following command:
```
screen -S run_scrnabox
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 1
```

The resulting output files are deposited into `~/working_directory/step1`

**Parameter names need to be updated.**

### Step 2: Create Seurat object 
In this step, CellRanger-generated feature-barcode expression matrices are transformed into Seurat objects. Genes expressed in less than a minimum number of cells and cells expressing less than a minimum number of genes can be filtered.<br />
<br />
The following parameters are adjustable for Step 2:

|Parameter|Default|Description|
|:--|:--|:--|
|Save_RNA| No| Whether or not to export an RNA expression matrix|
|Save_metadata| No|Whether or not to export a metadata dataframe|
|count_matrices| NULL|If users skipped Step 1, the may provide the path to a directory that contains existing feature-barcode expression matrices to initiate the pipeline at Step 2 |
|min.cells_L| 0|Only retain genes expressed in a minimum number of cells|
|min.features_L| 0|Only retain cells expressing a minimum number of genes|

To run Step 2, use the following command:
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 2
```

The resulting output files are deposited into `~/working_directory/step2`

**Parameter names need to be updated and the code for ambient RNA detection must be added.**


### Step 3: Quality control and filtering
Low quality cells are filtered based on the user-defined thresholds for the number of genes detected per cell, number of unique transcripts detected per cell, percentage of mitochondrial-encoded transcripts, and percentage of ribosomal-encoded transcripts. Mitochondrial- and ribosomal-encoded genes can be filtered out, as well as a custom user-defined list of genes. Finally, normalization and scaling is performed on the individual Seurat objects prior to cell-cycle scoring. <br />
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


To run Step 3, use the following command:
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 3 
```
The resulting output files are deposited into `~/working_directory/step3`


### Step 4: Demultiplexing and doublet removal
Seurat’s implementation (MULTIseqDemux) of the tag assignment algorithm outlined in Multi-seq is used to demultiplex pooled samples and identify doublets according to the expression matrices of the sample-specific barcodes (McGinnis et al 2019).

The following parameters are adjustable for Step 4:

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

The resulting output files are deposited into `~/working_directory/step4`

### Step 5: Integration and linear dimensional reduction
Individual Seurat objects are integrated to enable the joint analysis across sequencing runs using Seurat's integration algorithm (Stuart et al. 2019); if experiments are limited to a single sequencing run, the integration Step can be bypassed. Normalization, scaling, and Linear dimensional reduction is then performed on the resulting Seurat object to inform the optimal parameters for clustering in Step 6.<br />
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
The resulting output files are deposited into `~/working_directory/step5`

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

### Step 7: Cluster annotation
Cell populations, or clusters, with similar expression profiles are annotated to define cell types by three distinct methods:<br />
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
 <br />
 In addition, scRNAbox allows users to **visualize the expression of select features** at the cluster or cell level via a user-provided list of gene identifiers in order to further inform the cellular species that make up specific clusters. 

The following parameters are adjustable for Step 7:

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

Please note that in order to test the cluster marker genes for enrichment across EnrichR libraries, the HPC must have access to the internet. If your HPC cannot access the internet, it is possible to run the enrichment step on your local machine directly in R. To do so, begin by downloading the Step 6 Seurat RDS object (`~/working_directory/step6/objs6/seu_step6.rds `)  and cluster marker RDS object (`~/working_directory/step7/info7/marker/ClusterMarkers.rds`) to your computer using the following commands:

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

### Step 8: Differential gene expression contrasts
This step performs differetial gene expression (DGE) analyses according to user-defined contrasts. Contrast can be performed between samples (e.g. case vs control; **sample-sample contrasts**) or between samples, stratified by cell type (e.g. case vs control for excitatory neurons only; **sample-cell contrasts**).
 
The following parameters are adjustable for Step 8:

|Parameter|Default|Description (the cluster annotation method associated with the parameter is shown)|
|:--|:--|:--|
|par_save_RNA| No| Whether or not to export an RNA expression matrix|
|par_save_metadata| No|Whether or not to export a metadata dataframe|
|par_level_cluster|integrated_snn_res.0.7| The cluster resolution that you used in the cluster annotation module. If you skipped integration in Step 5, use par_level_cluster='RNA_snn_res.0.7', if you want to proceed with a clustering resolution of 0.7|
|par_step8_clus_label|NULL|List of user-currated cluster labels obtained from the annotation module. Make sure to have the same number of labels as clusters at the desired clustering resolution.|
|par_new_genotype|yes|Whether or not you want to add new sample labels to simplify the contrasts. For example, you may wish to set both control1 and control2 as control.|
|par_old_antibody_label|NULL|list of old sample labels (i.e. those used to create the samples_info folder in the setup procedures)|
|par_new_antibody_label|NULL|list of new sample labels corresponding to the old sample labels defined in the parameter above|


The first step of the DGE contrasts module is to create a  DGEListobject. This step may require alot of RAM, we suggest 3*size(seu_int_clu.rds). Users can adjust the RAM in the `scrnabox_config.ini` file located in `~/working_directory/job_info/configs`. To create a DGEListobject, run the following command:

```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--dgelist T
```

The next step is to perform DGE contrasts. As mentioned, useers can perform **sample-sample contrasts** or **sample-cell contrasts**

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

**Note:** that if you have a large number of contrasts to run, it may be more efficient to split them up and submit batch jobs instead.

## Integrating seurat objects
To combine different seurat objects, you can run the following codes. 
```
LISTOFSEU=~/list.txt

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps  integrate \
--seulist ${LISTOFSEU}
```

LISTOFSEU includes the path of seurate files, put them in different lines. By default, the pipeline standardize the seurat objects before integrating, you can change the default in `${SCRNABOX_PWD}/job_info/parameters/stepint_par.txt`.
