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
## HTO analysis track: PBMC dataset
## Contents

- [Introduction](#introduction)
- [Downloading the pbmc dataset](#downloading-the-pbmc-dataset)
- [Installation](#installation)
    - [scrnabox.slurm installation](#scrnaboxslurm-installation)
    - [CellRanger installation](#cellranger-installation)
    - [R library preparation and R package installation](#r-library-preparation-and-r-package-installation)
- [scRNAbox: HTO Analysis Track](#scrnaboxpipeline)
    - [Step 0: Set up](#step-0-set-up)
    - [Step 1: FASTQ to gene expression matrix](#step-1-fastq-to-gene-expression-matrix)  
    - [Step 2: Create Seurat object and remove ambient RNA ](#step-2-create-seurat-object-and-remove-ambient-rna)  
    - [Step 3: Quality control and filtering](#step-3-quality-control-and-filtering)
    - [Step 4: Demultiplexing and doublet detection](#step-4-demultiplexing-and-doublet-detection)      
- [Publication-ready figures](#publication-ready-figures)
- [Job Configurations](#job-configurations) 
 - - - -

## Introduction 
This guide illustrates the steps taken for our analysis of the PBMC dataset in our [pre-print manuscript](https://www.biorxiv.org/content/10.1101/2023.11.13.566851v1). Here, we are using the HTO analysis track of scRNAbox to analyze a publicly available scRNAseq dataset produced by [Stoeckius et al.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1603-1). This data set describes peripheral blood mononuclear cells (PBMC) from eight human donors, which were tagged with sample-specific barcodes, pooled, and sequenced together in a single run. 
 
- - - -
## Downloading the PBMC dataset
In you want to use the PBMC dataset to test the scRNAbox pipeline, please see [here](pbmc_download.md) for detialed instructions on how to download the publicly available data.
 - - - -

## Installation
#### scrnabox.slurm installation
To download the latest version of `scrnabox.slurm` (v0.1.52) run the following command: 
```
wget https://github.com/neurobioinfo/scrnabox/releases/download/v0.1.52/scrnabox.slurm.zip
unzip scrnabox.slurm.zip
```
For a description of the options for running `scrnabox.slurm` run the following command:
```
bash /pathway/to/scrnabox.slurm/launch_scrnabox.sh -h 
```

If the `scrnabox.slurm` has been installed properly, the above command should return the folllowing:
```
scrnabox pipeline version 0.1.52
------------------- 
        mandatory arguments:
                -d  (--dir)  = Working directory (where all the outputs will be printed) (give full path)
                --steps  =  Specify what steps, e.g., 2 to run just step 2, 2-4, run steps 2 through 4)

        optional arguments:
                -h  (--help)  = See helps regarding the pipeline arguments. 
                --method  = Select your preferred method: HTO and SCRNA for hashtag, and Standard scRNA, respectively. 
                --msd  = You can get the hashtag labels by running the following code (HTO Step 4). 
                --markergsea  = Identify marker genes for each cluster and run marker gene set enrichment analysis (GSEA) using EnrichR libraries (Step 7). 
                --knownmarkers  = Profile the individual or aggregated expression of known marker genes. 
                --referenceannotation  = Generate annotation predictions based on the annotations of a reference Seurat object (Step 7). 
                --annotate  = Add clustering annotations to Seurat object metadata (Step 7). 
                --addmeta  = Add metadata columns to the Seurat object (Step 8). 
                --rundge  = Perform differential gene expression contrasts (Step 8). 
                --seulist  = You can directly call the list of Seurat objects to the pipeline.  
 
                --rcheck  = You can identify which libraries are not installed.  
 
 ------------------- 
 For a comprehensive help, visit  https://neurobioinfo.github.io/scrnabox/site/ for documentation.
```
 - - - -

#### CellRanger installation 

For information regarding the installation of CellRanger, please visit the 10X Genomics [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation). If CellRanger is already installed on your HPC system, you may skip the CellRanger installation procedures.

For our analysis of the midbrain dataset we used the 10XGenomics GRCh38-3.0.0 reference genome and CellRanger v5.0.1. For more information regarding how to prepare reference genomes for the CellRanger _counts_ pipeline, please see the 10X Genomics [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#grch38_3.0.0).

 - - - -

#### R library preparation and R package installation
We must prepapre a common R library where we will load all of the required R packages. If the required R packages are already installed on your HPC system in a common R library, you may skip the following procedures. 
<br />

We will first install `R`. The analyses presented in our [pre-print manuscript]() were conducted using v4.2.1.

```
# install R
module load r/4.2.1
```

Then, we will run the installation code, which creates a directory where the R packages will be loaded and will install the required R packages:

```
# Folder for R packages 
R_PATH=~/path/to/R/library
mkdir -p $R_PATH

# Install package
Rscript ./scrnabox.slurm/soft/R/install_packages_scrnabox.R $R_PATH
```
 - - - -

## scRNAbox pipeline
### Step 0: Set up
Now that `scrnabox.slurm`, `CellRanger`, `R`, and the required R packages have been installed, we can proceed to our analysis with the scRNAbox pipeline. We will create a `pipeline` folder designated for the analysis and run Step 0, selecting the HTO analysis track (`--method HTO`), using the following code:
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
Next, we will navigate to the `scrnabox_config.ini` file in `~/pipeline/job_info/configs` to define the HPC account holder (**ACCOUNT**), the path to the environmental module (**MODULEUSE**), the path to CellRanger from the environmental module directory (**CELLRANGER**), CellRanger version (**CELLRANGER_VERSION**), R version (**R_VERSION**), and the path to the R library (**R_LIB_PATH**):

```
cd ~/pipeline/job_info/configs
nano scrnabox_config.ini

ACCOUNT=account-name
MODULEUSE=/path/to/environmental/module 
CELLRANGER=/path/to/cellranger/from/module/directory 
CELLRANGER_VERSION=5.0.1
R_VERSION=4.2.1
R_LIB_PATH=/path/to/R/library
```

Next, we can check to see if all of the required R packages have been properly installed using the following command:

```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 0 \
--rcheck 
```
 - - - -

### Step 1: FASTQ to gene expression matrix
In Step 1, we will run the CellRanger _counts_ pipeline to generate feature-barcode expression matrices from the FASTQ files. While it is possible to manually prepare the `library.csv` and `feature_ref.csv` files for the sequencing run prior to running Step 1, for this analysis we are going to opt for automated library preparation. For more information regarding the manual prepartion of `library.csv` and `feature_ref.csv` files, please see the the [CellRanger library preparation](library_prep.md) tutorial. <br /> 
<br /> 
For our analysis of the PBMC dataset we set the following execution parameters for Step 1 (`~/pipeline/job_info/parameters/step1_par.txt`):

|Parameter|Value|
|:--|:--|
|par_automated_library_prep|Yes|
|par_fastq_directory|/path/to/directory/contaning/fastqs|
|par_RNA_run_names|run1GEX|
|par_HTO_run_names|run1HTO|
|par_seq_run_names|run1|
|par_paired_end_seq|Yes|
|par_id|Hash1, Hash2, Hash3, Hash4, Hash5, Hash6, Hash7, Hash8|
|par_name|A_TotalSeqA, B_TotalSeqA, C_TotalSeqA, D_TotalSeqA, E_TotalSeqA, F_TotalSeqA, G_TotalSeqA, H_TotalSeqA|
|par_read|R2|
|par_pattern|5P(BC)|
|par_sequence|AGGACCATCCAA, ACATGTTACCGT, AGCTTACTATCC, TCGATAATGCGA, GAGGCTGAGCTA, GTGTGACGTATT, ACTGTCTAACGG, TATCACATCGGT|
|par_ref_dir_grch|~/genome/10xGenomics/refdata-cellranger-GRCh38-3.0.0|
|par_r1_length|NULL (commented out)|
|par_r2_length|NULL (commented out)|
|par_mempercode|30|
|par_include_introns|NULL (commented out)|
|par_no_target_umi_filter|NULL (commented out)|
|par_expect_cells|NULL (commented out)|
|par_force_cells|NULL (commented out)| 
|par_no_bam|NULL (commented out)| 


**Note:** The parameters file for each step is located in `~/pipeline/job_info/parameters`. For a comprehensive description of the execution parameters for each step see [here](reference.md). 

Given that CellRanger runs a user interface and is not submitted as a Job, it is recommended to run Step 1 in a **'screen'** which will allow the the task to keep running if the connection is broken. To run Step 1, use the following command:
```
export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=~/pipeline

screen -S run_PBMC_application_case
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 1
```
The outputs of the CellRanger _counts_ pipeline are deposited into `~/pipeline/step1`.

 - - - -

### Step 2: Create Seurat object and remove ambient RNA 
In Step 2, we are going to begin by correcting the RNA assay for ambient RNA removal using SoupX ([Young et al. 2020](https://academic.oup.com/gigascience/article/9/12/giaa151/6049831)). We will then use the the ambient RNA-corrected feature-barcode matrices to create a Seurat object.  <br /> 
<br /> 
For our analysis of the PBMC dataset we set the following execution parameters for Step 2 (`~/pipeline/job_info/parameters/step2_par.txt`):

|Parameter|Value|
|:--|:--|
|par_save_RNA| Yes|
|par_save_metadata| Yes|
|par_ambient_RNA| Yes|
|par_normalization.method|LogNormalize|
|par_scale.factor|10000|
|par_selection.method|vst|
par_nfeatures|2500|

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
step2
├── figs2
│   ├── ambient_RNA_estimation_run1.pdf
│   ├── ambient_RNA_markers_run1.pdf
│   ├── cell_cyle_dim_plot_run1.pdf
│   ├── vioplot_run1.pdf
│   └── zoomed_in_vioplot_run1.pdf
├── info2
│   ├── estimated_ambient_RNA_run1.txt
│   ├── MetaData_1.txt
│   ├── meta_info_1.txt
│   ├── run1_ambient_rna_summary.rds
│   ├── sessionInfo.txt
│   ├── seu1_RNA.txt
│   └── summary_seu1.txt
├── objs2
│   └── run1.rds
└── step2_ambient
    └── run1
        ├── barcodes.tsv
        ├── genes.tsv
        └── matrix.mtxs 
```

**Note:** For a comprehensive description of the outputs for each analytical step, please see the [Outputs](outputs.md) section of the scRNAbox documentation.

<p align="center">
<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/784c866c-7ee2-46e6-a106-763c0693554f"> <br /> 
</p>

**Figure 1. Figures produced by Step 2 of the scRNAbox pipeline.** **A)** Estimated ambient RNA contamination rate (Rho) by SoupX. Estimates of the RNA contamination rate using various estimators are visualized via a frequency distribution; the true contamination rate is assigned as the most frequent estimate (red line; 8.7%). **B)** Log10 ratios of observed counts to expected counts for marker genes from each cluster. Clusters are defined by the CellRanger counts pipeline. The red line displays the estimated RNA contamination rate if the estimation was based entirely on the corresponding gene. **C)** Principal component analysis (PCA) of Seurat S and G2M cell cycle reference genes. **D)** Violin plots showing the distribution of cells according to quality control metrics calculated in Step 2. **E)** Zoomed in violin plots, from the minimum to the mean, showing the distribution of cells according to quality control metrics calculated in Step 2.

 - - - -

### Step 3: Quality control and filtering
In Step 3, we are going to perform quality control procedures and filter out low quality cells. We are going to filter out cells with < 50 unique RNA transcripts, > 6000 unique RNA transcripts, < 200 total RNA transcripts, > 7000 total RNA transcripts, and > 50% mitochondria.

For our analysis of the PBMC dataset we set the following execution parameters for Step 3 (`~/pipeline/job_info/parameters/step2_par.txt`):

|Parameter|Value|
|:--|:--|
|par_save_RNA| Yes| 
|par_save_metadata| Yes|
|par_seurat_object| NULL |
|par_nFeature_RNA_L|50 |
|par_nFeature_RNA_U|6000 |
|par_nCount_RNA_L|200 |
|par_nCount_RNA_U|7000 |
|par_mitochondria_percent_L|0 | 
|par_mitochondria_percent_U|50 |
|par_ribosomal_percent_L|0 |
|par_ribosomal_percent_U|100 |
|par_remove_mitochondrial_genes|No| 
|par_remove_ribosomal_genes|No| 
|par_remove_genes|NULL|
|par_regress_cell_cycle_genes|Yes|
|par_normalization.method|LogNormalize|
|par_scale.factor|10000|
|par_selection.method|vst|
|par_nfeatures|2500|
|par_top|10|
|par_npcs_pca|30|

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
│   ├── dimplot_pca_run1.pdf
│   ├── elbowplot_run1.pdf
│   ├── filtered_QC_vioplot_run1.pdf
│   └── VariableFeaturePlot_run1.pdf
├── info3
│   ├── MetaData_run1.txt
│   ├── meta_info_run1.txt
│   ├── most_variable_genes_run1.txt
│   ├── run1_RNA.txt
│   ├── sessionInfo.txt
│   └── summary_run1.txt
└── objs3
    └── run1.rds
```
<p align="center">
<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/76148dee-c567-4203-9e26-05baa3a71704"> <br /> 
</p>

**Figure 2. Figures produced by Step 3 of the scRNAbox pipeline.** **A)** Violin plots showing the distribution of cells according to quality control metrics  after filtering by user-defined thresholds. **B)** Scatter plot showing the top 2500 most variable features; the top 10 most variable features are labelled. **C)** Principal component analysis (PCA) visualizing the first two principal component (PC). **D)** Elbow plot to visualize the percentage of variance explained by each PC.  

 - - - -

### Step 4: Demultiplexing and doublet detection
In Step 4, we are going to demultiplex the pooled samples and remove doublets (erroneous libraries produced by two or more cells) based on the expression of the sample-specific barcodes (antibody assay). 

If the barcode labels used in the analysis are unknown, the first step is to retrieve them from the Seurat object. To do this, we do not need to modify the execution parameters and can go straight to running the following code:

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
|par_scale.factor| 10000|
|par_selection.method|vst|
|par_nfeatures|2500|
|par_dimensionality_reduction|Yes|
|par_npcs_pca|30|
|par_dims_umap|3|
|par_n.neighbor|65|
|par_dropDN|Yes|
|par_label_dropDN|Doublet, Negative|
|par_quantile|0.9|
|par_autoThresh|TRUE|
|par_maxiter|5|
|par_RidgePlot_ncol|3|
|par_old_antibody_label|A-TotalSeqA, B-TotalSeqA, C-TotalSeqA, D-TotalSeqA, E-TotalSeqA, F-TotalSeqA, G-TotalSeqA, H-TotalSeqA, Doublet|
|par_new_antibody_label|sample-A, sample-B, sample-C, sample-D, sample-E, sample-F, sample-G, sample-H, Doublet|

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
│   ├── run1_DotPlot_HTO_MSD.pdf
│   ├── run1_Heatmap_HTO_MSD.pdf
│   ├── run1_HTO_dimplot_pca.pdf
│   ├── run1_HTO_dimplot_umap.pdf
│   ├── run1_nCounts_RNA_MSD.pdf
│   └── run1_Ridgeplot_HTO_MSD.pdf
├── info4
│   ├── run1_filtered_MULTIseqDemuxHTOcounts.csv
│   ├── run1_MetaData.txt
│   ├── run1_meta_info_.txt
│   ├── run1_MULTIseqDemuxHTOcounts.csv
│   ├── run1_RNA.txt
│   └── sessionInfo.txt
└── objs4
    └── run1.rds
```

<p align="center">
<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/693cfc0e-39e7-4421-adcc-2e666d1996bd"> <br /> 
</p>


**Figure 3. Figures produced by Step 4 of the Cell Hashtag Analysis Track.** **A)** Uniform Manifold Approximation and Projections (UMAP) plot, taking the first three pricipal components (PC) of the antibody assay as input. **B)** Principal component analysis (PCA) showing the first two PCs of the antibody assay.  **C)** Ridgeplot visualizing the enrichment of barcode labels across sample assignments at the sample level. **D)** Dot plot visualizing the enrichment of barcode labels across sample assignments at the sample level. **E)** Heatmap visualizing the enrichment of barcode labels across sample assignments at the cel level. **D)** Violin plot visualizing the distribution of the number of total RNA transcripts identified per cell, startified by sample assignment.

 - - - -

### Publication-ready figures
The code used to produce the publication-ready figures used in our [pre-print manuscript]() is avaliable here [here](https://github.com/neurobioinfo/scrnabox/tree/main/tutorial/Midbrain_dataset_figures).  

- - - -

### Job Configurations
The following job configurations were used for our analysis of the PBMC dataset. Job Configurations can be modified for each analytical step in the `scrnabox_config.ini` file in `~/pipeline/job_info/configs`

|Step |THREADS_ARRAY|MEM_ARRAY|WALLTIME_ARRAY|
|:--|:--|:--|:--|
|Step2|4|16g|00-05:00|
|Step3|4|16g|00-05:00|
|Step4|4|16g|00-05:00|



