---
layout: post
title: Run pipeline on processed data
description: A short introduction how run scRNA pipeline on processed data
date: 2022-02-23
author: Saeid Amiri
published: true
tags: scRNA 
categories: 
comments: false
---
## Application Case 2: Cell Hashtag scRNAseq Analysis Track of scRNAbox
## Contents

- [Introduction](#introduction)
- [Downloading the midbrain dataset](#downloading-the-midbrain-dataset)
- [Installation](#installation)
    - [scrnabox.slurm installation](#scrnaboxslurm-installation)
    - [CellRanger installation](#cellranger-installation)
    - [R library preparation and R package installation](#r-library-preparation-and-r-package-installation)
- [scRNAbox: Standard Analysis Track](#scrnabox-standard-analysis-track)
    - [Step 0: Pipeline initiation](#step-0-pipeline-initiation)
    - [Step 1: FASTQ to gene expression matrix](#step-1-fastq-to-gene-expression-matrix)  
    - [Step 2: Create Seurat object and remove ambient RNA ](#step-2-create-seurat-object-and-remove-ambient-rna)  
    - [Step 3: Quality control and filtering](#step-3-quality-control-and-filtering)
    - [Step 4: Demultiplexing and doublet detection](#step-4-demultiplexing-and-doublet-detection) 
    - [Step 5: Integration and linear dimensional reduction](#step-5-integration-and-linear-dimensional-reduction)
    - [Step 6: Clustering](#step-6-clustering) 
    - [Step 7: Cluster annotation](#step-7-cluster-annotation) 
        - [Method 1: Cluster marker GSEA](#method-1-cluster-marker-gsea) 
        - [Method 2: Module score](#method-2-module-score) 
        - [Method 3: Reference-based annotation](#method-3-reference-based-annotation) 
        - [Visualizing the expression of select features](#visualizin-the-expression-of-select-features)    
    - [Step 8: Differential gene expression contrasts](#step-8-differential-gene-expression-contrasts)
        - [Create DGEList object](#create-dgelist-object) 
        - [Sample-sample contrasts](#sample-sample-contrasts) 
        - [Sample-cell contrasts](#sample-cell-contrasts)
 - [Job Configurations](#job-configurations)      
 - - - -

## Introduction 
This guide illustrates the steps taken for Application Case 2 in our pre-print manuscript. Here, we are using the Cell Hashtag scRNAseq Analysis Track of scRNAbox to analyze a publicly available scRNAseq dataset produced by [Stoeckius et al.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1603-1). This data set (referred to as the PBMC dataset in the manuscript) describes >20,000 peripheral blood mononuclear cells (PBMC) from eight human donors, which were tagged with sample-specific barcodes, pooled, and sequenced together in a single run. 
 - - - -

## Downloading the PBMC dataset
The scRNAseq data produced by [Stoeckius et al.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1603-1) is publicly available in the Gene Expression Omnibus with accession code [GSE108313](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108313). To download the data, we must first install SRAtoolkit (if this is not already installed on your High-Performance Computing (HPC) system). We will create a directory for our raw data and download SRAtoolkit with the following code:

```
mkdir data_download
cd data_download
wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -vxzf sratoolkit.tar.gz
export PATH=$PATH:$PWD/sratoolkit.3.0.5-ubuntu64/bin
```

For more information regarding the SRAtoolkit, please visit the [documentation](https://github.com/ncbi/sra-tools/wiki).

The Sequence Read Archive (SRA) run identifiers for the RNA and antibody assays are:

|Assay|SRR|
|:--|:--|
|RNA|SRR8281306|
|Antibody|SRR8281307|

To download the FASTQ files for the RNA and antibody assays, run the following code. Please note that this may take a very long time. 

```
export PATH=$PATH:$PWD/sratoolkit.3.0.5-ubuntu64/bin
module load StdEnv/2020 gcc/9.3.0
module load sra-toolkit/3.0.0 

#RNA
prefetch SRR8281306 --max-size 100GB
fasterq-dump SRR8281306

#Antibody
prefetch SRR8281307 --max-size 100GB
fasterq-dump SRR8281307 

```
If the FASTQ files for the RNA and antibody assays have been downloaded properly, the `data_download` folder should contain the following:
```
data_download
├── SRR8281306
│   └── SRR8281306.sra
├── SRR8281306_1.fastq
├── SRR8281306_2.fastq
├── SRR8281307
│   └── SRR8281307.sra
├── SRR8281307_1.fastq
└── SRR8281307_2.fastq
```


Next, we will rename the FASTQ files according to the CellRanger nomenclature and transfer the FASTQ files to a folder named `fastqs`. For more information regarding the nomeclature required by the CellRanger _counts_ pipeline, please visit CellRanger's [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/fastq-input).

**Note**: The `fastqs` folder should only contain FASTQ files for the experiment.

```
mkdir fastqs

# RNA assay
cp ~/data_download/SRR8281306_1.fastq  ~/fastqs/run1GEX_S1_L001_R1_001.fastq
cp ~/data_download/SRR8281306_2.fastq  ~/fastqs/run1GEX_S1_L001_R2_001.fastq

# HTO assay
cp ~/data_download/SRR8281307_1.fastq ~/fastqs/run1HTO_S1_L001_R1_001.fastq
cp ~/data_download/SRR8281307_2.fastq ~/fastqs/run1HTO_S1_L001_R2_001.fastq

```
If the above steps were conducted properly, the `fastqs` folder should contain the following files:
```
fastqs
├── run1GEX_S1_L001_R1_001.fastq
├── run1GEX_S1_L001_R2_001.fastq
├── run1HTO_S1_L001_R1_001.fastq
└── run1HTO_S1_L001_R2_001.fastq
```

 - - - -

## Installation
#### scrnabox.slurm installation
Now that the raw data has been downloaded and organized, we can install the latest version of `scrnabox.slurm` (v0.135):
```
wget https://github.com/neurobioinfo/scrnabox/releases/download/v0.1.35/scrnabox.slurm.zip
unzip scrnabox.slurm.zip
```
For a description of the options for running `scrnabox.slurm` run the following command:
```
export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
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
 - - - -

#### CellRanger installation 

For information regarding the installation of CellRanger, please visit the 10X Genomics [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation). If CellRanger is already installed on your HPC system, you may skip the CellRanger installation procedures.

For our analysis of the midbrain dataset we used the 10XGenomics GRCh38-3.0.0 reference genome and CellRanger v5.0.1. For more information regarding how to prepare reference genomes for the CellRanger _counts_ pipeline, please see the 10X Genomics [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#grch38_3.0.0).

 - - - -

#### R library preparation and R package installation
We must prepapre a common R library where we will load all of the required R packages. If the required R packages are already installed on your HPC system in a common R library, you may skip the following procedures. 
First, we will creat an `R` folder and download our desired R version. The analyses presented in our pre-print manuscript were conducted using R v4.2.1

```
#make common R library
mkdir R_library
cd R_library

#install and open R in the terminal
module load r/4.2.1
R

#set common R library path
R_LIB_PATH="~/R_library"
.libPaths(R_LIB_PATH)

library(Seurat)
library(ggplot2)
library(dplyr)
library(foreach)
library(doParallel)
library(Matrix)
library(cowplot)
library(clustree)
library(xlsx)
library(enrichR)
library(stringi)
library(limma)
library(tidyverse)
library(edgeR)
library(vctrs)
library(RColorBrewer)
library(fossil)
library(openxlsx)
library(stringr)
library(ggpubr)
library(SoupX)
library(MatrixGenerics)
library(BiocGenerics)
library(S4Vectors)
library(IRanges)
library(GenomeInfoDb)
library(GenomicRanges)
library(Biobase)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(DropletUtils)  
library(stringr)
devtools::install_github(“neurobioinfo/scrnabox/scrnaboxR”)
```

**Saeid, is there a way to automatically load these packages when scrnaboxR is installed?**

 - - - -

## scRNAbox: Standard Analysis Track
### Step 0: Pipeline initiation
Now that `scrnabox.slurm`, `CellRanger`, `R`, and the Required R packages have been installed, we can proceed to our analysis with the Standard scRNAseq Analysis Track of the scRNAbox pipeline. We will create a `pipeline` folder designated for the analysis and run the pipeline initiation Step using the following code:
```
mkdir pipeline
cd pipeline

export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=~/pipeline

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 0 \
--method HTO
```
Next, we will navigate to the `scrnabox_config.ini` file in `~/pipeline/job_info/configs` to define the path to the R library (`R_LIB_PATH=`), the version of R (`R_VERSION=`), and the path to CellRanger (`MODULECELLRANGER=`):
```
cd ~/pipeline/job_info/configs
nano scrnabox_config.ini

MODULECELLRANGER=mugqic/cellranger/5.0.1
R_VERSION=4.2.1
R_LIB_PATH=~/R
```
 - - - -

### Step 1: FASTQ to gene expression matrix
In this Step, we will run the CellRanger _counts_ pipeline to generate feature-barcode expression matrices from the FASTQ files. While it is possible to manually prepare the `library.csv` and `feature_ref.csv` files for the sequencing run prior to running Step 1, for this analysis we are going to opt for automated library preparation. For more information regarding the manual prepartion of `library.csv` and `feature_ref.csv` files, please see the the [CellRanger library preparation](library_prep.md) tutorial. <br /> 
<br /> 
For our analysis of the PBMC dataset we set the following execution parameters for Step 1 (`~/pipeline/job_info/parameters/step1_par.txt`):

|Parameter|Value|
|:--|:--|
|par_automated_library_prep|Yes|
|par_fastq_directory|~/fastqs|
|par_RNA_run_names|run1GEX|
|par_HTO_run_names|run1HTO|
|par_seq_run_names|run1|
|par_paired_end_seq|Yes|
|id|Hash1, Hash2, Hash3, Hash4, Hash5, Hash6, Hash7, Hash8|
|name|A_TotalSeqA, B_TotalSeqA, C_TotalSeqA, D_TotalSeqA, E_TotalSeqA, F_TotalSeqA, G_TotalSeqA, H_TotalSeqA|
|read|R2|
|pattern|5P(BC)|
|sequence|AGGACCATCCAA, ACATGTTACCGT, AGCTTACTATCC, TCGATAATGCGA, GAGGCTGAGCTA, GTGTGACGTATT, ACTGTCTAACGG, TATCACATCGGT|
|REF_DIR_GRCH|~/genome/10xGenomics/refdata-cellranger-GRCh38-3.0.0|
|R1LENGTH|NULL|
|MEMPERCORE|30|


**Note:** The parameters file for each Analytical Step is located in `~/pipeline/job_info/parameters`. For a comprehensive description of the execution parameters for each Analytical Step see the [Execution parameters](reference.md) section of the scRNAbox documentation. 

Given that CellRanger runs a user interface and is not submitted as a Job, it is recommended to run Step 1 in a **'screen'** which will allow the the task to keep running if the connection is broken. To run Step 1, use the following command:
```
export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=~/pipeline

screen -S run_stoeckius_application_case
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 1
```
The outputs of the CellRanger _counts_ pipeline are deposited into `~/pipeline/step1`. The expression matrix, features, and barcode files outputed by CellRanger are located in `~/pipeline/step1/sample/ouput_folder/outs/raw_feature_bc_matrix`.

 - - - -

### Step 2: Create Seurat object and remove ambient RNA 
In this Step, we are going to begin by correcting the RNA assay for ambient RNA removal using SoupX (Young et al. 2020). We will then use the the ambient RNA-corrected feature-barcode matrices produced to create a Seurat object.  <br /> 
<br /> 
For our analysis of the PBMC dataset we set the following execution parameters for Step 2 (`~/pipeline/job_info/parameters/step2_par.txt`):

|Parameter|Value|
|:--|:--|
|par_save_RNA| Yes|
|par_save_metadata| Yes|
|par_ambient_RNA| Yes|
|par_count_matrices| NULL|
|par_min.cells_L| 0|
|par_min.features_L|0|

We can run Step 2 using the following code:
```
export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=~/pipeline

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 2
```

Step 2 produces the following outputs:
```
~/pipeline
├── step2
│   ├── figs2
│   │   └── vioplot_seu1.png
│   ├── info2
│   │   ├── meta_info_seu1.txt
│   │   ├── run1_ambient_rna_summary.rds
│   │   ├── sessionInfo.txt
│   │   ├── seu1_RNA.txt
│   │   ├── seu_MetaData1.txt
│   │   └── summary_seu1.txt
│   └── objs2
│       └── seu1.rds
└── step2_ambient 
    └── run1
    ├── barcodes.tsv
    ├── genes.tsv
    └── matrix.mtx   
```

**Note:** For a comprehensive description of the outputs for each Analytical Step, please see the [Outputs](outputs.md) section of the scRNAbox documentation.

<p align="center">
<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/78a382d3-3399-483c-99d3-7aba18c84689" width="350" height="100"> <br /> 
</p>

**Figure 1. Figure produced by Step 2 of the Cell Hashtag Analysis Track.** Violin plots are produced to visualize the distribution of genes per cell (nFeature_RNA), molecules per cell (nCount_RNA), percentage of mitochondrial-encoded genes per cell (percent.mt), and perentage of ribosomal genes per cell (percent.ribo) for each sequencing run.

 - - - -

### Step 3: Quality control and filtering
In this Step, we are going to perform quality control (QC) procedures and filter out low quality cells. We are going to filter out cells expressing <200 features, cells with <200 unqiue molecules, and cells with >10% mitochondrial-encoded genes. In addition, we are going to remove mitochondrial-encoded and ribosomal genes and will perform cell cycle scoring. Prior to performing cell cycle scoring, we must normalize and scale the RNA assay counts matrix. 

For our analysis of the PBMC dataset we set the following execution parameters for Step 3 (`~/pipeline/job_info/parameters/step2_par.txt`):

|Parameter|Value|
|:--|:--|
|par_save_RNA| Yes| 
|par_save_metadata| Yes|
|par_seurat_object| NULL |
|par_nFeature_RNA_L|NULL |
|par_nFeature_RNA_U|NULL |
|par_nCount_RNA_L|1500 |
|par_nCount_RNA_U|NULL |
|par_mitochondria_percent_L|NULL | 
|par_mitochondria_percent_U|10 |
|par_ribosomal_percent_L|NULL |
|par_ribosomal_percent_U|10 |
|par_remove_mitochondrial_genes|Yes| 
|par_remove_ribosomal_genes|Yes| 
|par_remove_genes|NULL|
|par_normalization.method|LogNormalize|
|par_scale.factor|10000|
|par_selection.method|vst|
|par_nfeatures|2500|
|par_top|10|
|par_npcs_pca|30|
|par_cells|500|
|par_dims|12|
|par_dims_umap|10|
|par_n.neighbors|65|

We can run Step 3 using the following code:
```
export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=~/pipeline

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 3
```

Step 3 produces the following outputs.
```
step3
├── figs3
│   ├── cellcycle_seu1.png
│   ├── dimplot_pca_seu1.png
│   ├── dimplot_umap_seu1.png
│   ├── elbowplot_seu1.png
│   ├── QC_vioplot_seu1.png
│   └── VariableFeaturePlot_seu1.png
├── info3
│   ├── MetaDataseu1.txt
│   ├── meta_info_seu1.txt
│   ├── most_variable_genes_seu1.txt
│   ├── sessionInfo.txt
│   ├── seu1_RNA.txt
│   └── summary_seu1.txt
└── objs3
    └── seu1.rds
```
<p align="center">
<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/104f3c00-3a5a-45ef-bca7-b4fa7e0b2e88" width="800" height="120"> <br /> 
</p>

**Figure 2. Figures produced by Step 3 of the Cell hashtag Analysis Track.** **A)** Distribution of QC metrics after filtering according to the user-defined thresholds. **B)** Variable features plot showing the top 2500 most variable features; the top 10 most variable features are labelled. **C)** Elbow plot to visualize the percentage of variance explained by each principal component (PC). **D)** Principal component analysis (PCA) visualizing the first two PCs. **E)** Uniform Manifold Approximation and Projections (UMAP) plot, taking the first ten PCs as input. **F)** Distibution of G2M and S scores across cells.

 - - - -

### Step 4: Demultiplexing and doublet detection
In this Step, we are going to demultiplex the pooled samples and remove doublets (erroneous libraries produced by two or more cells) based on the expression matrices of the sample-specific barcodes (antibody assay). 

If the barcode labels used in the analysis are unkown, the first step is to retrieve them from the Seurat object. To do this, we do not need to modify the execution parameters and can go straight to running the following code:

```
export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=~/pipeline

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 4 \
--msd T 
```

The above code produces the following file: 
```
step4
├── figs4
├── info4
│   └── seu1.rds_old_antibody_label_MULTIseqDemuxHTOcounts.csv
└── objs4
```
Which contains the names of the barcode labels (i.e. **A_TotalSeqA**, **B_TotalSeqA**, **C_TotalSeqA**, **D_TotalSeqA**, **E_TotalSeqA**, **F_TotalSeqA**, **G_TotalSeqA**, **H_TotalSeqA**, **Doublet**, **Negative**).

Now that we know the barcode labels used in the PBMC dataset, we can perform demultiplexing and doublet detection.

For our analysis of the PBMC dataset we set the following execution parameters for Step 4 (`~/pipeline/job_info/parameters/step4_par.txt`):

|Parameter|Value|
|:--|:--|
|par_save_RNA|Yes|
|par_save_metadata|Yes|
|par_normalization.method|CLR|
|par_scale.factor| 1000|
|par_selection.method|vst|
|par_nfeatures|2500|
|par_dims_umap|3|
|par_n.neighbor|65|
|par_dimensionality_reduction|Yes|
|par_npcs_pca|30|
|par_dropDN|Yes|
|par_label_dropDN|Doublet, Negative|
|par_quantile|0.9|
|par_autoThresh|TRUE|
|par_maxiter|5|
|par_RidgePlot_ncol|3|
|par_old_antibody_label|A-TotalSeqA, B-TotalSeqA, C-TotalSeqA, D-TotalSeqA, E-TotalSeqA, F-TotalSeqA, G-TotalSeqA, H-TotalSeqA, Doublet, Negative|
|par_new_antibody_label|sample-A, sample-B, sample-C, sample-D, sample-E, sample-F, sample-G, sample-H, Doublet, Negative|

We can run Step 4 using the following code:

```
export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=~/pipeline

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 4
```

Step 4 produces the following outputs.

```
step4
├── figs4
│   ├── HTO_dimplot_pca_seu1.png
│   ├── HTO_dimplot_umap_seu1.png
│   ├── seu1_DotPlotHTOIDMSD.png
│   ├── seu1_HeatmapHTOIDMSD.png
│   ├── seu1_nCounts_demult_groupsMSD.png
│   └── seu1_RidgeplotHTOMSD.png
├── info4
│   ├── meta_infoseu1.txt
│   ├── sessionInfo.txt
│   ├── seu1filtered_MULTIseqDemuxHTOcounts.csv
│   └── seu1_MULTIseqDemuxHTOcounts.csv
└── objs4
    └── seu1.rds
```

<p align="center">
<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/8eb35514-c77e-4b98-9615-5bd1c3976a91" width="650" height="100"> <br /> 
</p>


**Figure 3. Figures produced by Step 4 of the Cell Hashtag Analysis Track.** **A)** Principal component analysis (PCA) showing the first two principal components (PC) of the antibody assay. **B)** Uniform Manifold Approximation and Projections (UMAP) plot, taking the first three PCs of the antibody assay as input. **C)** Ridgeplot visualizing the enrichment of barcode labels across sample assignments at the sample level. **D)** Dot plot visualizing the enrichment of barcode labels across sample assignments at the sample level. **E)** Heatmap visualizing the enrichment of barcode labels across sample assignments at the cel level. **D)** Violin plot visualizing the distribution of unique molecules identified per cell, startified by sample assignment.

 - - - -

### Step 5: Integration and linear dimensional reduction
In this Step, we would typically integrate the individual Seurat objects for each sequencing run to enable joint analyses. However, because the PBMC dataset only consists of one sequencing run, we are going to skip integration. We will then perform normalization, scaling and linear dimensional reduction on the RNA assay. The outputs from Step 5 will inform the optimal clustering parameters for Step 6. 

For our analysis of the PBMC dataset we set the following execution parameters for Step 5:

|Parameter|Value|
|:--|:--|
|par_save_RNA| Yes| 
|par_save_metadata| Yes|
|par_skip_integration| Yes| 
|par_FindIntegrationAnchors_dim|NULL|
|par_normalization_and_scaling|no **this only controls normalization so adress the parameters if better**|
|par_DefaultAssay|RNA|
|par_normalization.method|LogNormalize|
|par_scale.factor|1000|Scale factor for scaling the data|
|par_selection.method|vst|
|par_nfeatures|2000|
|par_RunUMAP_n.neighbors|65|
|par_RunPCA_npcs|30| 
|par_RunUMAP_dims|30| 
|par_compute_jackstraw |Yes|

We can run Step 5 using the following code:

```
export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=~/pipeline

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 5
```
Step 5 produces the following outputs.
```
step5
├── figs5
│   ├── DimPlot_pca.png
│   ├── DimPlot_umap.png
│   ├── elbowPlot.png
│   └── Jackstraw_plot.png
├── info5
│   ├── meta_info_seu_step5.csv
│   ├── MULTI_IDcount.csv
│   ├── sessionInfo.txt
│   ├── seu_int_MetaData.txt
│   └── seu_int_RNA.txt
└── objs5
    └── seu_step5.rds
```
<p align="center">
<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/81a6e7dd-ae3f-4ebb-a084-3fb2daaf4c73" width="550" height="100"> <br /> 
</p>

**Figure 4. Figures produced by Step 5 of the Cell Hashtag Analysis Track.** **A)** Principal component analysis (PCA) visualizing the first two principal components (PC) of the integrated assay. **B)** Uniform Manifold Approximation and Projections (UMAP) plot of the integrated assay, taking the first ten PCs as input.  **C)** Jackstraw plot to visualize the distribution of p-values for each PC. **D)** Elbow plot to visualize the percentage of variance explained by each PC.

 - - - -

### Step 6: Clustering
In this Step, we will cluster the cells to indentify groups of cells with similar expression profiles. Based on the Elbow and Jackstraw plots produced in Step 5, we are going to use the first 15 PCs for the nearest-neighbour graph construction and to run the Uniform Manifold Approximation and Projection (UMAP) dimensional reduction. We will cluster the cells at a clustering resolution of 0.1 to 1.5, in intervals of 1.5. To determine the stability of clusters, we will run the Louvain clustering algorithm five times for each clustering resolution, while shuffling the order of the nodes in the graph for each iteration, and will compute the Adjusted Rand Index (ARI) between pairs of clusters at a given clustering resolution.

For our analysis of the PBMC dataset we set the following execution parameters for Step 6:

|Parameter|Value|
|:--|:--|
|par_save_RNA| Yes| 
|par_save_metadata| Yes|
|par_skip_step5|No|
|par_FindNeighbors_dims|15| 
|par_RunUMAP_dims|15| 
|par_FindNeighbors_k.param|30|
|par_FindNeighbors_prune.SNN|1/15|
|par_FindClusters_resolution|0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5|
|par_compute_ARI|Yes|
|par_RI_reps|5|

We can run Step 6 using the following code:

```
export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=~/pipeline

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 6
```
Step 6 produces the following outputs.
```
step6
├── ARI
│   ├── ARI.png
│   └── test.xlsx
├── figs6
│   ├── clustree_int.png
│   ├── RNA_snn_res.0.1.png
│   ├── RNA_snn_res.0.2.png
│   ├── RNA_snn_res.0.3.png
│   ├── RNA_snn_res.0.4.png
│   ├── RNA_snn_res.0.5.png
│   ├── RNA_snn_res.0.6.png
│   ├── RNA_snn_res.0.7.png
│   ├── RNA_snn_res.0.8.png
│   ├── RNA_snn_res.0.9.png
│   ├── RNA_snn_res.1.1.png
│   ├── RNA_snn_res.1.2.png
│   ├── RNA_snn_res.1.3.png
│   ├── RNA_snn_res.1.4.png
│   ├── RNA_snn_res.1.5.png
│   └── RNA_snn_res.1.png
├── info6
│   ├── meta_info.csv
│   ├── sessionInfo.txt
│   ├── seu_MetaData.txt
│   └── seu_RNA.txt
└── objs6
    └── seu_step6.rds
```

<p align="center">
<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/099fd2e4-c48a-4b7d-a2b4-230ce4dfc819" width="750" height="100"> <br /> 
</p>


**Figure 5. Figures produced by Step 6 of the Cell Hashtag Analysis Track.** **A)** ClustTree plot to visualize inter-cluster dynamics at varying cluster resolutions. **B)** Mean (top panel) and standard deviation (sd; middle panel) of the Adjusted RNA Index (ARI) between clustering pairs at each user-defined clustering resolution. The bottom panel shows the number of clusters at each user-defined clustering resolution. **C)** Uniform Manifold Approximation and Projections (UMAP) plot at a clustering resolution of 1.4. 

 - - - -

### Step 7: Cluster annotation
In this Step, we are going to annotate the clusters identified in Step 6 to define the cellular species in the midbrain dataset. scRNAbox provides three distinct methods for cluster annotations

- **Method 1:** Cluster marker gene set enrichment analysis (GSEA)
- **Method 2:** Module score
- **Method 2:** Reference-based annotation

In addition to these three Methods, we can **visualize the expression of select features** to further inform the cellular species in the dataset.

For comprehensive description of each cluster annotation Method, please see the [Cell Hashtag scRNAseq Analysis Track](HTO.md) section of the scRNAbox documentation or our pre-print manuscript. 

- - - -
#### Method 1: Cluster marker GSEA
Using Method 1, we are first going to identify differentially expressed marker genes for each cluster. We must define the number of marker genes for each cluster that we want scRNAbox to report and selecte a clustering resolution that we want to annotate. In this case we will report the top five marker genes for each cluster at a clustering resolution of 1.5.

To identify the marker genes for each cluster, we set the following execution parameters for Step 7 (`step7_par.txt`):

|Parameter|Value|
|:--|:--|
|par_save_RNA| Yes| 
|par_save_metadata| Yes|
|par_level_cluster| RNA_snn_res.1.4|
|par_level_genotype|MULTI_ID_Lables|
|par_top_sel|5|
|par_db|NULL|
|par_compute_module_score|No|
|par_module_score|NULL|
|par_reference|NULL|
|par_level_celltype|NULL|
|par_FindTransferAnchors_dim|NULL|
|par_futureglobalsmaxSize|NULL|
|par_visualize_select_features|No|
|par_select_features|NULL|

We can identify the marker genes for each cluster using the following code:

```
export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=~/pipeline

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--marker T
```

The above code produces the following outputs:

```
step7
├── figs7
│   ├── marker
│   │   └── heatmap.pdf
│   ├── module_score
│   ├── umap.pdf
│   ├── umap_splitted.pdf
│   └── visualize_select_features
├── info7
│   ├── marker
│   │   ├── cluster_just_genes.xlsx
│   │   ├── ClusterMarkers.csv
│   │   ├── ClusterMarkers.rds
│   │   ├── cluster_whole.xlsx
│   │   └── top_sel.csv
│   ├── module_score
│   └── sessionInfo_marker.txt
└── objs7
```

In addition to identifying the marker genes for each cluster, the above code produces UMAP plots at the user-defined clustering resolution (1.4) to visualize the clustering landscape across all cells in the dataset and at the sample level. 

<p align="center">
<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/9b619b29-24f5-4630-9641-f47a618ec1b5" width="750" height="100"> <br /> 
</p>

**Figure 6. Uniform Manifold Approximation and Projections (UMAP) plots at the user-defined clustering resolution.** **A)** The clustering landscape at the user-defined clustering resolution across all cells in the dataset. **B**) The clustering landscape at the user-defined clustering resolution, stratified by sample. 

Now that we have identified the marker genes for each cluster, we will perform a **gene set enrichment analysis (GSEA)**; we will test the differentially expressed genes (DEG) in the positive direction (Log2 fold-change > 0.00) for enrichment across gene set libraries that define cell types using the EnrichR tool. For this analysis, we will leverage the following libraries:

- Descartes_Cell_Types_and_Tissue_2021; 
- CellMarker_Augmented_2021; 
- Azimuth_Cell_Types_2021 cell type libraries.

To perform GSEA, we set the following execution parameters for Step 7 (`step7_par.txt`):

|Parameter|Value|
|:--|:--|
|par_save_RNA| Yes| 
|par_save_metadata| Yes|
|par_level_cluster| integrated_snn_res.1.5|
|par_top_sel|5|
|par_db|Descartes_Cell_Types_and_Tissue_2021,<br /> CellMarker_Augmented_2021,<br />Azimuth_Cell_Types_2021|
|par_compute_module_score|No|
|par_module_score|NULL|
|par_reference|NULL|
|par_level_celltype|NULL|
|par_FindTransferAnchors_dim|NULL|
|par_futureglobalsmaxSize|NULL|
|par_visualize_select_features|No|
|par_select_features|NULL|

If your HPC **allows access to the internet**, we can perform GSEA using the following code:

```
export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=~/pipeline

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--enrich T
```

**Note:** If your HPC **does not allow access to the internet**, you will have to run GSEA locally. For more information, please see the Standard scRNAseq [documentation](SCRNA.md) under the Step 7 section.

The above code produces the following outputs. As an example, we are only showing the outputs for cluster 0.
```
step7
└── annot_enrich
    ├── cluster0
    │   ├── Er.genes.1.csv
    │   ├── Er.genes.2.csv
    │   ├── Er.genes.3.csv
    │   ├── plotenrich1.pdf
    │   ├── plotenrich2.pdf
    │   └── plotenrich3.pdf
    ├── cluster1
    ├── cluster2
    ├── cluster3
    ├── cluster4
    ├── cluster5
    ├── cluster6
    ├── cluster7
    └── cluster8

```
After performing cluster marker GSEA and curating the results, we can produce our first iteration of the cluster annotations.<br />

**Note:** Visualizing intermediate cluster annotations is not incorporated into the scRNAbox pipeline; however we provide the code to do so below. Once users are satisfied with their final cluster annotations, they can provide the curated results in the parameters file for [Step 8](), which will be discussed below. 

```
## load and open R
module load r/4.2.1
R

## load parameters
# path to common R library
r_lib_path = "~/R"
# path to pipeline directory
output_dir = "~/pipeline"
# clustering resolution to cluster
clustering_resolution = "RNA_snn_res.1.5"
# intermediate cluster annotations
intermediate_cluster_labels = c("CD4 T", "CD4 T", "CD16 mono", "CD14 mono", "CD4 T", "B", "CD16 mono", "CD4 T/CD8 T", "NK/DC" )


## load library
.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr','stringi','limma','tidyverse','edgeR')
lapply(packages, library, character.only = TRUE)

## load Step 6 Seurat RDS object
sample_name<-list.files(path = paste(output_dir, "/step6/objs6",sep=""),pattern = "*.rds")
seu.int.c<-readRDS(paste(output_dir,'/step6/objs6/',sample_name, sep=''))

## set cluster annotations obtained from cluster annotations
cluster.ids<-intermediate_cluster_labels

## set cluster resolution and rename cluster identities
seu.int.c <- SetIdent(seu.int.c, value = clustering_resolution)
names(cluster.ids) <- levels(seu.int.c)    
seu.int.c <- RenameIdents(seu.int.c, cluster.ids) 

##plot UMAP
DimPlot(seu.int.c, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(file = paste(output_dir,'/step7/figs7','/intermediate_cluster_annotation.pdf', sep=''))
```

<p align="center">
<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/1c0d4342-3c8b-454d-bed8-fe71f31128b0" width="500" height="100"> <br /> 
</p>

**Figure 7. Figures produced by Method 1 (Cluster Marker GSEA) of the scRNAbox cluster annotation module.** **A)** The top expressional markers that define each cluster are visualized through a heatmap showing the expression  across cells, stratified by cluster. **B)** Differentially expressed marker genes in the positive direction (Log2 fold-change > 0.00) can be tested for enrichment across gene-set libraries that define cell types using the EnrichR tool. The enrichment results are visualized through a bar plot which displays the 20 most enriched terms for a particular cluster. As an example, we show the enrichment results of cluster 0 using the Azimuth_Cell_Types_2021 cell type library. **C)** Uniform Manifold Approximation and Projections (UMAP) plot showing the intermediate cluster annotations. 

- - - -
#### Visualizing the expression of select features
Now that we have broadly defined the cellular species that comprise our clusters, we are going to explore the expression of the marker genes used by [Smajic et al.](https://academic.oup.com/brain/article/145/3/964/6469020) to define their clusters:

|Cell tyep|Gene|
|:--|:--|
|Oligodendrocytes|_MOBP_|
|OPC|_VCAN_|
|Astrocytes|_AQP4_|
|Ependymal|_FOXJ1_|
|Microglia|_CD74_|
|Endothelial|_CLDN5_|
|Pericytes|_PDGFRB_|
|Excitatory neurons|_SLC17A6_|
|Inhibitory neurons|_GAD2_|
|GABAergic neurons|_GAD2_, _GRIK1_|
|Dopaminergic neurons (DaN)|_TH_|
|Degenerating DaN|_CADPS2_|

To visualize these features, we set the following execution parameters for Step 7 (`step7_par.txt`):

|Parameter|Value|
|:--|:--|
|par_save_RNA| Yes| 
|par_save_metadata| Yes|
|par_level_cluster| integrated_snn_res.1.5|
|par_top_sel|5|
|par_db|Descartes_Cell_Types_and_Tissue_2021,<br /> CellMarker_Augmented_2021,<br />Azimuth_Cell_Types_2021|
|par_compute_module_score|No|
|par_module_score|NULL|
|par_reference|NULL|
|par_level_celltype|NULL|
|par_FindTransferAnchors_dim|NULL|
|par_futureglobalsmaxSize|NULL|
|par_visualize_select_features|Yes|
|par_select_features|MOBP, VCAN, AQP4,FOXJ1, CD74, CLDN5, PDGFRB, SLC17A6, GAD2, GRIK1, TH, CADPS2|

We can visualize the expression of these features using the following code:

```
export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=~/pipeline

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--marker T
```

The above code produces the following outputs:
```
step7
└── figs7
    └── visualize_select_features
        ├──select_feature_dot_plot.pdf 
        ├──select_feature_feature_plot.pdf
        └──select_feature_violin_plot.pdf
```


<p align="center">
<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/f6150d60-fcf4-4426-9e56-7e60f94eedd1" width="900" height="100"> <br /> 
</p>

**Figure 8. Figures produced by the "visualize select features" option in the scRNAbox cluster annotation module.**  The expression of select features can be visualized at the cluster level via **A)** a dot plot and **B)** violin plots. **C)** The expression of select features can be visualized at the cell level via feature plots.

Based on the results of the above analyses, we can re-visit our cluster annotations using the same intermediate annotation code as above and visualize the annotations via a UMAP.

<p align="center">
<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/d682e588-0e51-4c8b-832e-384dc428a467" width="400" height="100"> <br /> 
</p>

**Figure 9. Intermediate cluster annotations.** Uniform Manifold Approximation and Projections (UMAP) plot showing the intermediate cluster annotations after leveraging the "visualize select feature" function of scRNAbox's cluster annotation module.

- - - -
#### Method 2: Module score
Using Method 2, we are going to comparatively quantify the expression of gene sets across clusters at the single-cell level. We will first define the gene set that we want to explore in a csv file; as an example, we are going to explore the expression of some well-known marker genes for the cellular species of interest. 

We will first produce a csv file with the following structure. This csv can be found **HERE**

|da_neurons|	NPC_orStemLike|	mature_neurons|	excitatory_neurons|	inhbitory_neurons	|astrocytes|	oligodendrocytes |	radial_glia|	epithelial| 	microglia| 
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
|TH|	DCX	|RBFOX3	|GRIA2|GAD1	|GFAP|	MBP|	PTPRC|	HES1 |	IBA1| 
|SLC6A3|	NEUROD1|	SYP|	GRIA1|	GAD2|	S100B|	MOG	|AIF1|	HES5|	P2RY12| 
|SLC18A2|	TBR1|	VAMP1	|GRIA4|	GAT1|	AQP4|	OLIG1|	ADGRE1|	SOX2|	P2RY13| 
|SOX6	|PCNA	|VAMP2|	GRIN1|	PVALB	|APOE	|OLIG2|	VIM|SOX10	|TREM119| 
|NDNF	|MKI67|	TUBB3	|GRIN2B|	GABR2|	SOX9|	SOX10|	TNC|	NES	|GPR34| 
|SNCG	|SOX2|	SYT1|	GRIN2A|	GABR1|	SLC1A3|	 |	PTPRZ1|	CDH1|	SIGLECH| 
|ALDH1A1|	NES|	BSN	|GRIN3A	|GBRR1|		|	|FAM107A	|NOTCH1	|TREM2| 
|CALB1|	PAX6|	HOMER1|	GRIN3|	GABRB2|	|	|	HOPX|		|CX3CR1| 
|TACR2|		|SLC17A6|	GRIP1	|GABRB1	|		| |LIFR	|	|FCRLS| 
|SLC17A6|	|		|CAMK2A	|GABRB3	|	|	|ITGB5	|	|OLFML3| 
|SLC32A1|	|		|	|GABRA6		|	| |IL6ST	|	|HEXB| 
|OTX2	|		|	| |GABRA1	|		| |SLC1A3	|	|TGFBR1| 
|GRP	|		|	| |GABRA4	|				| | | |SALL1| 
|LPL	|		|	| |TRAK2	|				| | | |MERTK| 
|CCK	|		|	|	|				| | | | |PROS1| 
|VIP	|		|	|	|				| | | | | 

We can then define the location of our csv file in the execution parameters for Step 7 (`step7_par.txt`):

|Parameter|Value|
|:--|:--|
|par_save_RNA| Yes| 
|par_save_metadata| Yes|
|par_level_cluster| integrated_snn_res.1.5|
|par_top_sel|5|
|par_db|NULL|
|par_compute_module_score|Yes|
|par_module_score|~/pipeline/module_score.csv|
|par_reference|NULL|
|par_level_celltype|NULL|
|par_FindTransferAnchors_dim|NULL|
|par_futureglobalsmaxSize|NULL|
|par_visualize_select_features|No|
|par_select_features|NULL|

We can compute the module score for our gene sets using the following code:

```
export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=~/pipeline

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--marker T
```
The above code produces the following outputs:
```
step7
├── figs7
│   └── module_score
│       ├──module_score_astrocytes.png
│       ├──module_score_epithelial.png
│       ├──module_score_da_neurons.png
│       ├──module_score_excitatory_neurons.png
│       ├──module_score_inhbitory_neurons.png
│       ├──module_score_mature_neurons.png
│       ├──module_score_microglia.png
│       ├──module_score_NPC_orStemLike.png
│       ├──module_score_oligodendrocytes.png
│       └──module_score_radial_glia.png
└── info7
    └──module_score
       └──geneset_by_cluster.csv

```

<p align="center">
<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/d42626d0-9c41-49e8-a3a9-8f38c0f91345" width="900" height="100"> <br /> 
</p>

**Figure 10. Figures produced by Method 2 (Module score) of the scRNAbox cluster annotation module.** Uniform Manifold Approximation and Projections (UMAP) plots showing the module score across established marker genes for **A)** dopaminergic neurons, **B)** neural progenitor cells, **C)** mature neurons, **D)** excitatory neurons, **E)** inhibitory neurons, **F)** astrocytes, **G)** oligodendrocytes, **H)** radial glia, **I)** epithelial cells, and **J)** microglia.

- - - -
#### Method 3: Reference-based annotation
Using Method 3, we are going to leverage the cell-type annotations from a reference Seurat object and generate annotation predictions for the query dataset. For reference-based annotation we must define the path to a our reference Seurat object and the column of the reference Seurat object's metadata that contains the cell type annotations. For the midbrain dataset, we are going to use a reference Seurat object from [Kamath et al.](https://www.nature.com/articles/s41593-022-01061-1).

To perform reference-based annotations, we set the following execution parameters for Step 7 (`step7_par.txt`):

|Parameter|Value|
|:--|:--|
|par_save_RNA| Yes| 
|par_save_metadata| Yes|
|par_level_cluster| integrated_snn_res.1.5|
|par_top_sel|NULL|
|par_db|NULL|
|par_compute_module_score|No|
|par_module_score|NULL|
|par_reference|~/reference_seurat_object.rds|
|par_level_celltype|Cell_Type|
|par_FindTransferAnchors_dim|10|
|par_futureglobalsmaxSize|60000 * 1024^2|
|par_visualize_select_features|No|
|par_select_features|NULL|

We can perform reference-based annotations using the following code:
```
export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=~/pipeline

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--fta T
```
The above code produces the following outputs:
```
step7
├── figs7
│   └── reference_based_annotation
│       └──UMAP_transferred_labels.pdf
└── objs7
    └──seu_step7.rds
```

<p align="center">
<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/00f64588-79d8-4d77-909f-2251f4381860" width="400" height="100"> <br /> 
</p>


**Figure 11. Figure produced by Method 3 (reference-based annotation) of the scRNAbox cluster annotation module.** Uniform Manifold Approximation and Projections (UMAP) plots showing the cluster annotations from the reference Seurat object projected onto the query Seurat object.

- - - -

### Step 8: Differential gene expression contrasts
In this step we are going to perform differential gene expression (DGE) analysis between our samples. ScRNAbox faciltates DGE contrasts between samples (**sample-sample contrasts**) and between samples, stratified by cell type (**sample-cell contrasts**). The DGE contrasts module contains three components:

1) **Create DGEList object** <br />
2) **Sample-sample contrasts**<br />
3) **Sample-cell contrasts**<br />
- - - -
#### Create DGEList object
First, we are going to create a DGElist object. Before doing so, we must define our desired clustering resolution and the final cluster annotations informed by Step 7. We are also going to rename our samples in order to faciliate DGE contrasts.

To create a DGElist object, we set the following execution parameters for Step 8 (`step8_par.txt`):

|Parameter|value|
|:--|:--|
|par_save_RNA| Yes| 
|par_save_metadata| Yes|
|par_level_cluster|integrated_snn_res.1.5|
|par_step8_clus_label|Oligodendrocytes, Oligodendrocytes, Excitatory_Neurons, Oligodendrocytes, Oligodendrocytes, Oligodendrocytes, Oligodendrocytes, Oligodendrocytes, OPC, Endothelial_cells, Microglia,  Oligodendrocytes, Astrocytes, Excitatory_Neurons, Astrocytes, Pericytes, Endothelial_cells, Ependymal, GABAergic_neurons, Oligodendrocytes, Inhibitory_neurons, Oligodendrocytes, OPC|
|par_new_genotype|yes|
|par_old_sample_label|Control1, Control2, Control3, Control4, Control5, Control6, Parkinson1, Parkinson2, Parkinson3, Parkinson4, Parkinson5|
|par_new_sample_label|Control, Control, Control, Control, Control, Control, Parkinson, Parkinson, Parkinson, Parkinson, Parkinson|

**Note:** Cell names and sample names cannot have spaces. For example, do not write "Endothelial cells", instead write "Endothelial_cells".

We can create the DGElist object using the following code:

```
export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=~/pipeline

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--dgelist T
```

The above code produces the following outputs:

```
step8
├── figs8
│   └── final_cluster_annotation.pdf
├── info8
│   ├── de_genes.rds
│   ├── dge.rds
│   ├── meta_info_de_genes.txt
│   ├── meta_info_dge.txt
│   └── meta_info_seu_step8.txt
└── objs8
    └── seu_step8.rds
```


<p align="center">
<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/13a7dc2f-56b2-4a0d-b034-fb50089582a9" width="400" height="100"> <br /> 
</p>


**Figure 12. Final cluster annotations used for differential gene expression (DGE) contrasts .** Uniform Manifold Approximation and Projections (UMAP) plots showing the final cluster annotation obtained by curating the results from scRNAbox's cluster annotation module (Step 7). The final cluster annotations will be used throughout the DGE contrasts module. 

- - - -
#### Sample-sample contrasts
Now that we have our DGElist object, we can perform DGE contrasts between samples (sample-sample contrasts). As an example, we will test for DGE between Parkinson's disease samples and controls. We must first define our contrast matrix in the sample-sample contrasts parameters file (`step8_contrast_genotype.txt`):

```
cont_name control versus
design1 Control Parkinson
```
We can perform sample-sample DGE contrasts using the following code:
```
export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=~/pipeline

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--genotype T
```
The above code produces the following outputs:

```
step8
└── cont_genotype
    └── design1.csv 
```

- - - -
#### Sample-cell contrasts
We can also perform DGE contrasts between samples, stratified by cell type (sample-cell contrasts). As an example, we will test for DGE between microglia from Parkinson's disease and controls. We must first define our contrast matrix in the sample-cell contrasts parameters file (`step8_contrast_celltype.txt`):

```
cont_name cell control versus
design1_cell Microglia Control Parkinson
```

We can perform sample-cel DGE contrasts using the following code:
```
export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=~/pipeline

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--celltype T
```
The above code produces the following outputs:

```
step8
└── cont_celltype
    └── design1.csv 
```
- - - -

### Step 9: Enrichment analysis

- - - -

### Job Configurations
The following job configurations were used for our analysis of the midbrain dataset. Job Configurations can be modified for each Analytical Step in the `scrnabox_config.ini` file in `~/pipeline/job_info/configs`

|Step |THREADS_ARRAY|MEM_ARRAY|WALLTIME_ARRAY|
|:--|:--|:--|:--|
|Step2|4|16g|00-05:00|
|Step3|4|16g|00-05:00|
|Step4|4|45g|00-05:00|
|Step5|4|45g|00-05:00|
|Step6|4|16g|00-05:00|
|Step7 marker|4|40g|00-01:00|
|Step7 fta|4|150g|00-09:00|
|Step8 dgelist|4|40g|00-12:00|
|Step8 cont|10|40g|00-12:00|

