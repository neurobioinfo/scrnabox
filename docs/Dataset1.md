---
layout: post
title: Run pipeline on processed data
description: A short introduction how run scRNA pipeline on processed data
date: 2022-02-23
author: Michael Fiorini
published: true
tags: scRNA 
categories: 
comments: false
---
## ScRNAbox analysis of the midbrain dataset
## Contents

- [Introduction](#introduction)
- [Downloading the midbrain dataset](#downloading-the-midbrain-dataset)
- [Installation](#installation)
    - [scrnabox.slurm installation](#scrnaboxslurm-installation)
    - [CellRanger installation](#cellranger-installation)
    - [R library preparation and R package installation](#r-library-preparation-and-r-package-installation)
- [scRNAbox: Standard Analysis Track](#scrnabox-standard-analysis-track)
    - [Step 0: Set up](#step-0-set-up)
    - [Step 1: FASTQ to gene expression matrix](#step-1-fastq-to-gene-expression-matrix)  
    - [Step 2: Create Seurat object and remove ambient RNA ](#step-2-create-seurat-object-and-remove-ambient-rna)  
    - [Step 3: Quality control and filtering](#step-3-quality-control-and-filtering)
    - [Step 4: Doublet removal](#step-4-doublet-removal) 
    - [Step 5: Integration](#step-5-integration)
    - [Step 6: Clustering](#step-6-clustering) 
    - [Step 7: Cluster annotation](#step-7-cluster-annotation) 
        - [Tool 1: Cluster marker GSEA](#tool-1-cluster-marker-gsea) 
        - [Tool 2: Expression profiling of known marker genes](#tool-2-expression-profiling-of-known-marker-genes) 
        - [Tool 3: Reference-based annotation](#tool-3-reference-based-annotation) 
        - [Annotate](#annotate)    
    - [Step 8: Differential gene expression](#step-8-differential-gene-expression-analysis)
        - [Add metadata](#add-metadata) 
        - [Cell-based DGE using all cells](#cell-based-dge-using-all-cells)
        - [Cell-based DGE using cell type groups](#cell-based-dge-using-cell-type-groups)
        - [Sample-based DGE using all cells](#sample-based-dge-using-all-cells)
        - [Sample-based DGE using cell type groups](#sample-based-dge-using-cell-type-groups)
 - [Analysis of differentially expressed genes](#analysis-of-differentially-expressed-genes)
 - [Publication-ready figures](#publication-ready-figures)
 - [Job Configurations](#job-configurations) 

 - - - -

## Introduction 
This guide illustrates the steps taken to analyze the midbrain dataset ([Smajic et al. 2022](https://academic.oup.com/brain/article/145/3/964/6469020)) that was presented in our [pre-print manuscript](). This dataset describes single-nuclei transcriptomes from the post-mortem midbrains of five individuals with Parkinson’s disease (PD) and six controls sequenced separately. 

 - - - -
## Downlaoding the midbrain dataset
In you want to use the midbrain dataset to test the scRNAbox pipeline, please see [here](midbrain_download.md) for detialed instructions on how to download the publicly available data.
 - - - -

## Installation
#### scrnabox.slurm installation
To download the latest version of `scrnabox.slurm` (v0.1.35) run the following command: 
```
wget https://github.com/neurobioinfo/scrnabox/releases/download/v0.1.35/scrnabox.slurm.zip
unzip scrnabox.slurm.zip
```
For a description of the options for running `scrnabox.slurm` run the following command:
```
bash /pathway/to/scrnabox.slurm/launch_scrnabox.sh -h 
```

If the `scrnabox.slurm` has been installed properly, the above command should return the folllowing:
```
        mandatory arguments:
                -d  (--dir)  = Working directory (where all the outputs will be printed) (give full path)
                --steps  =  Specify what steps, e.g., 2 to run just step 2, 2-4, run steps 2 through 4)

        optional arguments:
                -h  (--help)  = See helps regarding the pipeline arguments. 
                --method  = Select your preferred method: HTO and SCRNA for hashtag, and Standard scRNA, respectively. 
                --msd  = You can get the hashtag labels by running the following code 
                --markergsea  = Identify marker genes for each cluster and run marker gene set enrichment analysis (GSEA) using EnrichR libraries. 
                --knownmarkers  = Run module score and visualize the expression of known cell type marker genes. 
                --referenceannotation  = Run module score and visualize the expression of known cell type marker genes. 
                --annotate  = Run module score and visualize the expression of known cell type marker genes. 
                --addmeta  = Add metadata columns to the Seurat object 
                --rundge  = Perform differential gene expression contrasts 
                --seulist  = You can directly call the list of seurat objects to the pipeline.  
 
 ------------------- 
 For a comprehensive help, visit https://github.com/neurobioinfo/scrnabox for documentation.
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
Now that `scrnabox.slurm`, `CellRanger`, `R`, and the required R packages have been installed, we can proceed to our analysis with the scRNAbox pipeline. We will create a `pipeline` folder designated for the analysis and run Step 0, selecting the standard analysis track (`--method SCRNA`), using the following code:
```
mkdir pipeline
cd pipeline

export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=~/pipeline

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 0 \
--method SCRNA
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
 - - - -

### Step 1: FASTQ to gene expression matrix
In Step 1, we will run the CellRanger _counts_ pipeline to generate feature-barcode expression matrices from the FASTQ files. While it is possible to manually prepare the `library.csv` files for each of the 11 samples in the experiment prior to running Step 1, we are going to opt for automated library preparation. For more information regarding the manual prepartion of `library.csv` files, please see the the [CellRanger library preparation](library_prep.md) tutorial. <br /> 
<br /> 
For our analysis of the midbrain dataset we set the following execution parameters for Step 1 (`~/pipeline/job_info/parameters/step1_par.txt`):

|Parameter|Value|
|:--|:--|
|par_automated_library_prep|Yes|
|par_fastq_directory|/path/to/directory/containing/fastqs|
|par_sample_names|PD1, PD2, PD3, PD4, PD5, CTRL1, CTRL2, CTRL3, CTRL4, CTRL5, CTRL6|
|par_rename_samples|Yes|
|par_new_sample_names|Parkinson1, Parkinson2, Parkinson3, Parkinson4, Parkinson5, Control1, Control2, Control3, Control4, Control5, Control6|
|par_paired_end_seq|TRUE|
|par_ref_dir_grch|~/genome/10xGenomics/refdata-cellranger-GRCh38-3.0.0|
|par_r1_length|NULL (commented out)|
|par_r2_length|NULL (commented out)|
|par_mempercode|30|
|par_include_introns|Yes|
|par_no_target_umi_filter|NULL (commented out)|
|par_expect_cells|NULL (commented out)|
|par_force_cells|NULL (commented out)| 
|par_no_bam|NULL (commented out)| 

**Note:** The parameters file for each step is located in `~/pipeline/job_info/parameters`. For a comprehensive description of the execution parameters for each step see [here](reference.md). 

Given that CellRanger runs a user interface and is not submitted as a job, it is recommended to run Step 1 in a **'screen'**, which will allow the the task to keep running if the connection is broken. To run Step 1, use the following command:
```
export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=~/pipeline

screen -S run_smajic_application_case
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 1
```
The outputs of the CellRanger _counts_ pipeline are deposited into `~/pipeline/step1`.

 - - - -

### Step 2: Create Seurat object and remove ambient RNA 
In Step 2, we are going to use the CellRanger-generated feature-barcode matrices to produce unique Seurat objects for each of the 11 samples. In this step, we have the option to correct the expression matrices for ambient RNA contamination; however, because [Smajic et al.](https://academic.oup.com/brain/article/145/3/964/6469020) did not perform this analytical procedure we will skip it. In addition, we will perform cell cycle scoring. Prior to performing cell cycle scoring, we must normalize and scale the counts matrix. 
<br /> 

For our analysis of the midbrain dataset we set the following execution parameters for Step 2 (`~/pipeline/job_info/parameters/step2_par.txt`):

|Parameter|Value|
|:--|:--|
|par_save_RNA| Yes| 
|par_save_metadata| Yes|
|par_count_matrices| NULL|
|par_ambient_RNA| No|
|par_min.cells_L| 1|
|par_normalization.method|LogNormalize|
|par_scale.factor|10000|
|par_selection.method|vst|
|par_nfeatures|2500|


We can run Step 2 using the following code:
```
export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=~/pipeline

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 2
```

Step 2 produces the following outputs for each sample. As an example we show the outputs for sample _Control1_:
```
step2
├── figs2
│   ├── ambient_RNA_estimation_Control1.pdf
│   ├── ambient_RNA_markers_Control1.pdf
│   ├── cell_cyle_dim_plot_Control1.pdf
│   ├── vioplot_Control1.pdf
│   └── zoomed_in_vioplot_Control1.pdf
├── info2
│   ├── estimated_ambient_RNA_Control1.txt
│   ├── MetaData_Control1.txt
│   ├── meta_info_Control1.txt
│   ├── Control1_ambient_rna_summary.rds
│   ├── Control1_RNA.txt
│   ├── sessionInfo.txt
│   └── summary_Control1.txt
└── objs2
    └── Control1.rds
```

**Note:** For a comprehensive description of the outputs for each analytical step, please see the [Outputs](outputs.md) section of the scRNAbox documentation.

<p align="center">
<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/9cbbf8ea-8328-4915-9ebc-0a28022fe435"> <br /> 
</p>

**Figure 1. Figures produced by Step 2 of the scRNAbox pipeline.** The figures for the _Control1_ sample are shown as an example. **A)** Estimated ambient RNA contamination rate (Rho) by SoupX. Estimates of the RNA contamination rate using various estimators are visualized via a frequency distribution; the true contamination rate is assigned as the most frequent estimate (red line; 5.1%). **B)** Log10 ratios of observed counts to expected counts for marker genes from each cluster. Clusters are defined by the CellRanger counts pipeline. The red line displays the estimated RNA contamination rate if the estimation was based entirely on the corresponding gene. **C)** Principal component analysis (PCA) of Seurat S and G2M cell cycle reference genes. **D)** Violin plots showing the distribution of cells according to quality control metrics calculated in Step 2. **E)** Zoomed in violin plots, from the minimum to the mean, showing the distribution of cells according to quality control metrics calculated in Step 2.

 - - - -

### Step 3: Quality control and filtering
In Step 3, we are going to perform quality control procedures and filter out low quality cells. We are going to filter out cells with < 1000 unique RNA transcripts, < 1500 total RNA transcripts, and > 10% mitochondria and ribosomal RNA. In addition, we are going to remove mitochondrial-encoded and ribosomal genes.

For our analysis of the midbrain dataset we set the following execution parameters for Step 3 (`~/pipeline/job_info/parameters/step2_par.txt`):

|Parameter|Value|
|:--|:--|
|par_save_RNA| Yes| 
|par_save_metadata| Yes|
|par_seurat_object| NULL |
|par_nFeature_RNA_L|1000 |
|par_nFeature_RNA_U|NULL |
|par_nCount_RNA_L|1500 |
|par_nCount_RNA_U|NULL |
|par_mitochondria_percent_L|0 | 
|par_mitochondria_percent_U|10 |
|par_ribosomal_percent_L|0 |
|par_ribosomal_percent_U|10 |
|par_remove_mitochondrial_genes|Yes| 
|par_remove_ribosomal_genes|Yes| 
|par_remove_genes|NULL|
|par_regress_cell_cycle_genes|No|
|par_regress_custom_genes|No|
|par_regress_genes|NULL|
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

Step 3 produces the following outputs for each sample. As an example we show the outputs for sample _Control1_:
```
step3
├── figs3
│   ├── dimplot_pca_Control1.pdf
│   ├── elbowplot_Control1.pdf
│   ├── filtered_QC_vioplot_Control1.pdf
│   └── VariableFeaturePlot_Control1.pdf
├── info3
│   ├── MetaData_Control1.txt
│   ├── meta_info_Control1.txt
│   ├── most_variable_genes_Control1.txt
│   ├── Control1_RNA.txt
│   ├── sessionInfo.txt
│   └── summary_Control1.txt
└── objs3
    └── Control1.rds

```
<p align="center">
<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/9864c056-a1ea-4d77-af6f-dd1c72f7a59b"> <br /> 
</p>

**Figure 2. Figures produced by Step 3 of the scRNAbox pipeline.** The figures for the _Control1_ sample are shown as an example. **A)** Violin plots showing the distribution of cells according to quality control metrics  after filtering by user-defined thresholds. **B)** Scatter plot showing the top 2500 most variable features; the top 10 most variable features are labelled. **C)** Principal component analysis (PCA) visualizing the first two principal component (PC). **D)** Elbow plot to visualize the percentage of variance explained by each PC.  

 - - - -

### Step 4: Doublet removal
In this Step, we are going to identify doublets (erroneous barcodes produced by two or more cells) and remove them from downstream analyses using the [DoubletFinder](https://www.cell.com/cell-systems/pdfExtended/S2405-4712(19)30073-0) tool (McGinnis et al. 2019). For optimal performance, DoubletFinder requires the user to define the following parameters:

 - The number of statistically significant PCs (par_PCs)
 - The number of artificial doublets to generate (par_pN)
 - The expected doublet rate for each sample (par_expected_doublet_rate)

 The **number of statistically significant PCs** can be informed by the elbow plots produced in Step 3; in this case the top 25 PCs should maintain a robust compression of the data across samples. DoubletFinder is largely invariant to the **number of artifical doublets generated**, therefore we will maintain the default parameter of 0.25. The **expected doublet rate** can be informed by the number of recovered cells (~8% for ~10,000 cells recovered). The number of recovered cells can be informed by the `barcodes.tsv.gz` file produced by the CellRanger _counts_ pipeline, which is located in `~/pipeline/step1/<sample>/output_folder/outs/filtered_feature_bc_matrix`. 

The expected doublet rates are approximations obtained from the 10X Genomics Next GEM Single Cell 3' v3.1 [documentation](https://kb.10xgenomics.com/hc/en-us/articles/360059124751-Why-is-the-multiplet-rate-different-for-the-Next-GEM-Single-Cell-3-LT-v3-1-assay-compared-to-other-single-cell-applications-), which was used by [Smajic et al.](https://academic.oup.com/brain/article/145/3/964/6469020) for library preparation. 

For our analysis of the midbrain dataset we set the following execution parameters for Step 4 (`~/pipeline/job_info/parameters/step4_par.txt`):

|Parameter|Value|
|:--|:--|
|par_save_RNA| Yes| 
|par_save_metadata|Yes|
|par_seurat_object| NULL |
|par_RunUMAP_dims|25|
|par_RunUMAP_n.neighbors|65|
|par_dropDN| Yes| 
|par_PCs|25| 
|par_pN|0.25|
|par_sct|FALSE|
|par_sample_names|Control1, Control2, Control3, Control4, Control5, Control6, Parkinson1, Parkinson2, Parkinson3, Parkinson4, Parkinson5|
|par_expected_doublet_rate|0.042, 0.042, 0.023, 0.04, 0.027, 0.053, 0.023, 0.053, 0.034, 0.023, 0.05|

We can run Step 4 using the following code:

```
export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=~/pipeline

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 4
```
Step 4 produces the following outputs for each sample. As an example we show the outputs for sample _Control1_:
```
step4
├── figs4
│   ├── Control1_DF.classifications.pdf
│   └── Control1_doublet_summary.pdf
├── info4
│   ├── MetaData_Control1.txt
│   ├── meta_info_Control1.txt
│   ├── n_predicted_doublets_Control1.txt
│   ├── Control1_RNA.txt
│   └── sessionInfo.txt
└── objs4
    └── Control1.rds

```
<p align="center">
<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/6e16710f-9d0b-446f-9d07-0dcb4955bdd3"> <br /> 
</p>


**Figure 3. Figures produced by Step 4 of the scRNAbox pipeline standard track.** The figures for the _Control1_ sample are shown as an example. **A)** Results of doublet detection analysis with DoubletFinder. Left: violin plot displaying the distribution of the proportion of artificial nearest neighbours (pANN) across singlets and doublets. Right: a bar plot of the number of predicted singlets and doublets. **B)** Uniform Manifold Approximation Projection (UMAP) plots coloured by droplet assignments (singlet or doublet).

 - - - -

### Step 5: Integration
In Step 5, we are going to integrate the individual Seurat objects to enable joint analyses across all 11 samples. We will then perform normalization, scaling and linear dimensional reduction on the integrated assay. The outputs from Step 5 will inform the optimal clustering parameters for Step 6. 

For our analysis of the midbrain dataset we set the following execution parameters for Step 5:

|Parameter|Value|
|:--|:--|
|par_save_RNA| Yes| 
|par_save_metadata| Yes|
|par_seurat_object| NULL |
|par_one_seurat| No|
|par_integrate_seurat| Yes| 
|par_merge_seurat| No| 
|par_DefaultAssay|RNA|
|par_normalization.method|LogNormalize|
|par_scale.factor|10000|
|par_selection.method|vst|
|par_nfeatures|4000|
|par_FindIntegrationAnchors_dim|25|
|par_RunPCA_npcs|30|
|par_RunUMAP_dims|25|
|par_RunUMAP_n.neighbors|65|
|par_compute_jackstraw |yes|

We can run Step 5 using the following code:

```
export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=~/pipeline

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 5
```
Step 5 produces the following outputs:
```
step5
├── figs5
│   ├── integrated_DimPlot_pca.pdf
│   ├── integrated_DimPlot_umap.pdf
│   ├── integrated_elbow.pdf
│   └── integrated_Jackstraw_plot.pdf
├── info5
│   ├── int_meta_info_seu_step5.csv
│   ├── sessionInfo.txt
│   ├── seu_int_RNA.txt
│   └── seu_int_MetaData.txt
└── objs5
    └── seu_step5.rds
```
<p align="center">
<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/69cbb1c9-0090-4e51-aad0-ce145af9502a"> <br /> 
</p>


**Figure 4. Figures produced by Step 5 of the scRNAbox pipeline standard track.** **A)** Principal component analysis (PCA) visualizing the first two principal components (PC) of the integrated assay, colour coded by sample. **B)** Uniform Manifold Approximation and Projections (UMAP) plot of the integrated assay, colour coded by sample. **C)** Elbow plot to visualize the percentage of variance explained by each PC. **D)** Jackstraw plot to visualize the distribution of p-values for each PC. 

 - - - -

### Step 6: Clustering
In Step 6, we will cluster the cells to indentify groups with similar expression profiles. Based on the Elbow and Jackstraw plots produced in Step 5, we are going to use the first 25 PCs as input. We will cluster the cells at clustering resolutions ranging from 0.0 to 2.0. To determine the stability of clusters at each clustering resolution, we will run the Louvain clustering algorithm 25 times for each resolution, while shuffling the order of the nodes in the graph for each iteration. We will then compute the Adjusted Rand Index (ARI) between pairs of clusters at a given clustering resolution.

For our analysis of the midbrain dataset we set the following execution parameters for Step 6:

|Parameter|Value|
|:--|:--|
|par_save_RNA| Yes| 
|par_save_metadata| Yes|
|par_seurat_object| NULL |
|par_skip_integration|No|
|par_FindNeighbors_dims|25| 
|par_RunUMAP_dims|25| 
|par_FindNeighbors_k.param|30|
|par_FindNeighbors_prune.SNN|1/15|
|par_FindClusters_resolution|0, 0.05, 0.2, 0.6, 0.8, 1.0, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0|
|par_compute_ARI|Yes| 
|par_RI_reps|25|
 
We can run Step 6 using the following code:

```
export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=~/pipeline

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 6
```
Step 6 produces the following outputs:
```
step6
├── ARI
│   ├── ARI.pdf
│   └── clustering_ARI.xlsx
├── figs6
│   ├── clustree_int.pdf
│   ├── integrated_snn_res.0.05.pdf
│   ├── integrated_snn_res.0.2.pdf
│   ├── integrated_snn_res.0.6.pdf
│   ├── integrated_snn_res.0.8.pdf
│   ├── integrated_snn_res.0.pdf
│   ├── integrated_snn_res.1.2.pdf
│   ├── integrated_snn_res.1.4.pdf
│   ├── integrated_snn_res.1.5.pdf
│   ├── integrated_snn_res.1.6.pdf
│   ├── integrated_snn_res.1.8.pdf
│   ├── integrated_snn_res.1.pdf
│   └── integrated_snn_res.2.pdf
├── info6
│   ├── meta_info.csv
│   ├── sessionInfo.txt
│   ├── seu_MetaData.txt
│   └── seu_RNA.txt
└── objs6
    └── seu_step6.rds
```

<p align="center">
<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/7b4e633e-4995-4097-a994-ac6fbffa332f" width="750" height="100"> <br /> 
</p>


**Figure 5. Figures produced by Step 6 of the scRNAbox pipeline.** **A)** Uniform Manifold Approximation and Projections (UMAP) plot, coloured according to the clusters identified at a resolution of 1.5. UMAP plots are produced for each user-defined clustering resolution. **B)** Mean (top panel) and standard deviation (sd; middle panel) of the Adjusted RNA Index (ARI) between clustering pairs at each user-defined clustering resolution. The bottom panel shows the number of clusters at each user-defined clustering resolution. **C)** ClustTree plot to visualize inter-cluster dynamics at varying cluster resolutions.

 - - - -

### Step 7: Cluster annotation
In Step 7, we are going to annotate the clusters identified in Step 6 to identify the cell types comprising the midbrain dataset. scRNAbox provides three distinct tools for cluster annotations:

- **Tool 1:** Cluster marker and gene set enrichment analysis (GSEA)
- **Tool 2:** Expression profiling of known marker genes
- **Tool 3:** Reference-based annotation

Additionally, users can add cluster annotations to the Seurat object.

For comprehensive description of each cluster annotation tool, please see the [Step 7:Cluster annotation](Step7.md) section of the scRNAbox documentation or our [pre-print manuscript](). 

- - - -
#### Tool 1: Cluster marker GSEA
Using Tool 1, we are first going to identify differentially expressed marker genes for each cluster. We must define the number of marker genes for each cluster that we want scRNAbox to report and select a clustering resolution that we want to annotate. In this case we will report the top two marker genes for each cluster at a clustering resolution of 1.5.

To identify the marker genes for each cluster, we set the following execution parameters for Step 7 (`step7_par.txt`):

|Annotation tool|Parameter|Value|
|:--|:--|:--|
|**General**|par_save_RNA| No|
|**General**|par_save_metadata| No|
|**General**|par_seurat_object| NULL |
|**General**|par_level_cluster| integrated_snn_res.1.5| 
|**Tool 1**|par_run_find_marker|Yes|
|**Tool 1**|par_run_enrichR|No|
|**Tool 1**|par_top_sel|2|
|**Tool 1**|par_db|NULL|

We can identify the marker genes for each cluster using the following code:

```
export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=~/pipeline

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--markergsea T
```

The above code produces the following outputs:

```
step7
├── figs7
│   └── marker
│       └── heatmap.pdf
├── info7
│   ├── marker
│   │   ├── cluster_just_genes.xlsx
│   │   ├── ClusterMarkers.csv
│   │   ├── ClusterMarkers.rds
│   │   ├── cluster_whole.xlsx
│   │   └── top_sel.csv
│   └── sessionInfo_find_marker.txt
└── objs7
    └── seu_step7.rds
```

<p align="center">
<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/e7034065-0efa-48ab-8416-8e4792f86cc7"> <br /> 
</p>

**Figure 6. Figure produced by Tool 1 (marker gene gene set enrichment analysis) of scRNAbox's cluster annotation module (Step 7).** A heatmap is produced to visualize the expression of the top markers genes at the cell level, stratified by cluster. The top marker genes for each cluster identified at a clustering resolution of 1.5 is shown. 
- - - -

Now that we have identified the marker genes for each cluster, we will perform a **gene set enrichment analysis (GSEA)**; we will test the differentially expressed genes (DEG) in the positive direction (Log2 fold-change > 0.00) for enrichment across gene set libraries that define cell types using the EnrichR tool. For this analysis, we will use the following libraries:

- Descartes_Cell_Types_and_Tissue_2021; 
- CellMarker_Augmented_2021; 
- Azimuth_Cell_Types_2021 cell type libraries.

To perform GSEA, we set the following execution parameters for Step 7 (`step7_par.txt`):

|Annotation tool|Parameter|Value|
|:--|:--|:--|
|**General**|par_save_RNA| No|
|**General**|par_save_metadata| No|
|**General**|par_seurat_object| NULL |
|**General**|par_level_cluster| integrated_snn_res.1.5| 
|**Tool 1**|par_run_find_marker|No|
|**Tool 1**|par_run_enrichR|Yes|
|**Tool 1**|par_top_sel|2|
|**Tool 1**|par_db|NULL|

If your HPC **allows access to the internet**, we can perform GSEA using the following code:

```
export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=~/pipeline

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--markergsea T
```

**Note:** If your HPC **does not allow access to the internet**, you will have to run GSEA locally. For more information, please see the [Step 7:Cluster annotation](Step7.md) section of the scRNAbox documentation.

The above code produces the following outputs. As an example, we are only showing the outputs for cluster 0:
```
step7
└── annot_enrich
    └── cluster0
        ├── Er_genes_clust_0_Azimuth_Cell_Types_2021.csv
        ├── Er_genes_clust_0_CellMarker_Augmented_2021.csv
        ├── Er_genes_clust_0_Descartes_Cell_Types_and_Tissue_2021.csv
        ├── plotenrich_clust_0_1.pdf
        ├── plotenrich_clust_0_2.pdf
        └── plotenrich_clust_0_3.pdf

```
<p align="center">
<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/a29d66a3-2d1e-47ec-be8d-c7421ac10da6"> <br /> 
</p>

**Figure 7. Figures produced by Tool 1 (marker gene gene set enrichment analysis) of scRNAbox's cluster annotation module (Step 7).** Upon identifying the top marker genes for each cluster at the user-defined clustering resolution, users can perfrom a gene set enrichment analysis (GSEA) using the EnrichR tool for all marker genes in the positive direction (Log2 fold-change > 0.00). Bar plots are produced to visualize the most enriched terms for each cluster. As an example, the top enrichment results across **A)** Azimuth Cell Types 2021, **B)** Descartes Cell Types and Tissue 2021, and  **C)** CellMarker Augmented 2021 cell type libraries for cluster 5 are shown.  

**Note:** It is possible to identify cluster-specific marker genes and perform GSEA at the same time by setting both `par_run_find_marker= "Yes"` and `par_run_enrichR= "Yes`.

- - - -
#### Tool 2: Expression profiling of known marker genes
Using Tool 2, we are going to profile the midbrain dataset for known cell type marker genes. First, we are going to visualize the expression of the marker genes used by [Smajic et al.](https://academic.oup.com/brain/article/145/3/964/6469020) to define their clusters:

|Cell type|Gene|
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

|Annotation tool|Parameter|Value|
|:--|:--|:--|
|**General**|par_save_RNA| No|
|**General**|par_save_metadata| No|
|**General**|par_seurat_object| NULL |
|**General**|par_level_cluster| integrated_snn_res.1.5| 
|**Tool 2**|par_run_module_score|No|
|**Tool 2**|par_run_visualize_markers|Yes|
|**Tool 2**|par_module_score|NULL|
|**Tool 2**|par_select_features_list|MOBP, VCAN, AQP4, FOXJ1, CD74, CLDN5, GFRB, SLC17A6, GAD2, GRIK1, TH, CADPS2|
|**Tool 2**|par_select_features_csv|NULL|

We can visualize the expression of these features using the following code:

```
export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=~/pipeline

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--knownmarkers T
```

The above code produces the following outputs:
```
step7
├── figs7
    └── visualize_features
        ├── list_dot_plot.pdf
        ├── list_feature_plot.pdf
        └── list_violin_plot.pdf
```


<p align="center">
<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/5c87b8e1-1a2a-4ae5-8877-1486f6d2fa0d"> <br /> 
</p>

**Figure 8. Figures produced by Tool 2 (profiling the expression of known marker genes) of scRNAbox's cluster annotation module (Step 7).**  ScRNAbox allows users to visualize the individual and aggregated expression of known marker genes. For the midbrain dataset, we visualized the individual expression of the marker genes used by [Smajic et al.](https://academic.oup.com/brain/article/145/3/964/6469020) to define their clusters. Feature plots visualizing the expression of each individual gene at the cell level are shown. Additionally, violin plots and dot plots are produced to visualize the expression of individual genes at the cluster level. 

- - - -
Next, we are going to use a csv file to defineline multiple lists of cell type marker genes:

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

**Note**: This is csv file is available [here](https://github.com/neurobioinfo/scrnabox/blob/main/tutorial/Midbrain_dataset_example_files/module_score.csv).

We will use this csv file to visualize the individual expression of each gene and calculate the module score, which will allow us to visualize the aggregated expression of each gene set.

To visualize the unique and aggregated expression of these features, we set the following execution parameters for Step 7 (`step7_par.txt`):

|Annotation tool|Parameter|Value|
|:--|:--|:--|
|**General**|par_save_RNA| No|
|**General**|par_save_metadata| No|
|**General**|par_seurat_object| NULL |
|**General**|par_level_cluster| integrated_snn_res.1.5| 
|**Tool 2**|par_run_module_score|Yes|
|**Tool 2**|par_run_visualize_markers|Yes|
|**Tool 2**|par_module_score|~/scrnabox/tutorial/Midbrain_dataset_example_files/module_score.csv|
|**Tool 2**|par_select_features_list|NULL|
|**Tool 2**|par_select_features_csv|~/scrnabox/tutorial/Midbrain_dataset_example_files/module_score.csv|

We can visualize the individual and aggregated expression of these features using the following code:

```
export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=~/pipeline

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--knownmarkers T
```

The above code produces the following outputs:
```
step7
├── figs7
│   ├── module_score
│   │   ├── module_score_astrocytes.pdf
│   │   ├── module_score_da_neurons.pdf
│   │   ├── module_score_epithelial.pdf
│   │   ├── module_score_excitatory_neurons.pdf
│   │   ├── module_score_inhbitory_neurons.pdf
│   │   ├── module_score_mature_neurons.pdf
│   │   ├── module_score_microglia.pdf
│   │   ├── module_score_NPC_orStemLike.pdf
│   │   ├── module_score_oligodendrocytes.pdf
│   │   └── module_score_radial_glia.pdf
│   └── visualize_features
│       ├── astrocytes_dot_plot.pdf
│       ├── astrocytes_feature_plot.pdf
│       ├── astrocytes_violin_plot.pdf
│       ├── da_neurons_dot_plot.pdf
│       ├── da_neurons_feature_plot.pdf
│       ├── da_neurons_violin_plot.pdf
│       ├── epithelial_dot_plot.pdf
│       ├── epithelial_feature_plot.pdf
│       ├── epithelial_violin_plot.pdf
│       ├── excitatory_neurons_dot_plot.pdf
│       ├── excitatory_neurons_feature_plot.pdf
│       ├── excitatory_neurons_violin_plot.pdf
│       ├── inhbitory_neurons_dot_plot.pdf
│       ├── inhbitory_neurons_feature_plot.pdf
│       ├── inhbitory_neurons_violin_plot.pdf
│       ├── mature_neurons_dot_plot.pdf
│       ├── mature_neurons_feature_plot.pdf
│       ├── mature_neurons_violin_plot.pdf
│       ├── microglia_dot_plot.pdf
│       ├── microglia_feature_plot.pdf
│       ├── microglia_violin_plot.pdf
│       ├── NPC_orStemLike_dot_plot.pdf
│       ├── NPC_orStemLike_feature_plot.pdf
│       ├── NPC_orStemLike_violin_plot.pdf
│       ├── oligodendrocytes_dot_plot.pdf
│       ├── oligodendrocytes_feature_plot.pdf
│       ├── oligodendrocytes_violin_plot.pdf
│       ├── radial_glia_dot_plot.pdf
│       ├── radial_glia_feature_plot.pdf
│       └── radial_glia_violin_plot.pdf
└── info7
    └── module_score
        └── geneset_by_cluster.csv


```

<p align="center">
<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/9a612ece-ce9c-4000-a8ba-426f92eece20"> <br /> 
</p>

**Figure 9. Figures produced by Tool 2 (profiling the expression of known marker genes) of scRNAbox's cluster annotation module (Step 7).** ScRNAbox allows users to visualize the individual and aggregated expression of known marker genes. For the midbrain dataset, we used gene sets of well known marker genes for cell types in the human midbrain. The microglia gene set is shown as an example. **A)** Dot plot visualizing the individual expression of marker genes in the microglia gene set. **B)** Uniform manifold approximation and projection (UMAP) plot visualizing the module score of the microglia gene set at the cell level. The module score allows for profiling the aggregated expression of multiple genes. 

- - - -
#### Tool 3: Reference-based annotation
Using Method 3, we are going to leverage the cell-type annotations from a reference Seurat object and generate annotation predictions for the query dataset. For reference-based annotation we must define the path to a our reference Seurat object and the column of the reference Seurat object's metadata that contains the cell type annotations. For the midbrain dataset, we are going to use a reference Seurat object from [Kamath et al.](https://www.nature.com/articles/s41593-022-01061-1). The reference Seurat object was obtain from the [Broad Institute Single Cell Portal](https://singlecell.broadinstitute.org)

To perform reference-based annotations, we set the following execution parameters for Step 7 (`step7_par.txt`):

|Annotation tool|Parameter|Value|
|:--|:--|:--|
|**General**|par_save_RNA| Yes| 
|**General**|par_save_metadata| Yes|
|**General**|par_seurat_object| NULL |
|**General**|par_level_cluster| integrated_snn_res.1.5|
|**Tool 3**|par_reference|/path/to/Kamath/reference/seurat/object|
|**Tool 3**|par_reference_name|Kamath| 
|**Tool 3**|par_level_celltype|Cell_Type|
|**Tool 3**|par_FindTransferAnchors_dim|10| 
|**Tool 3**|par_futureglobalsmaxSize|60000 * 1024^2|


We can perform reference-based annotations using the following code:
```
export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=~/pipeline

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--referenceannotation T
```

The above code produces the following outputs:
```
step7
├── figs7
│   └── reference_based_annotation
│       └── Kammath_UMAP_transferred_labels.pdf
├── info7
│   ├── reference_based_annotation
│   │   └── Kammath_prediction_summary.xlsx
│   ├── sessionInfo_annotate.txt
│   ├── seu_MetaData.txt
│   └── seu_RNA.txt
└── objs7
    └── seu_step7.rds
```

<p align="center">
<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/ae1d99f8-f914-4f04-9737-afe8d35463c6"> <br /> 
</p>


**Figure 10. Figure produced by Tool 3 (reference-based annotation) of scRNAbox's cluster annotation module (Step 7).** ScRNAbox allows users to predict cell type annotations of their query dataset based on the predictions of a reference dataset. Left: Uniform Manifold Approximation and Projections (UMAP) plot of clustered and annotated reference Seurat object; single-nucleus RNA sequencing (scRRNAseq) of midbrain tissue produced by [Kamath et al.](https://www.nature.com/articles/s41593-022-01061-1), coloured by cell type. Right: UMAP of the label transfer predictions for the midbrain dataset, coloured by predicted cell type. Abbreviations: astro, astrocytes; da, dopaminergic neurons; endo, endothelial cells; nonda, non-dopaminergic neurons; mg, microglia; olig, oligodendrocytes; opc, oligodendrocyte precursor cells.

- - - -
#### Annotate
Now that we have used the three cluster annotation tools available in the scRNAbox pipeline, we are going to manually curate the results and add cluster annotations to the midbrain dataset.

To add annotations, we set the following execution parameters for Step 7 (`step7_par.txt`):

|Annotation tool|Parameter|Value|
|:--|:--|:--|
|**General**|par_save_RNA| Yes| 
|**General**|par_save_metadata| Yes|
|**General**|par_seurat_object| NULL |
|**General**|par_level_cluster| integrated_snn_res.1.5|
|**Annotate**|par_annotate_resolution|integrated_snn_res.1.5| 
|**Annotate**|par_name_metadata|clustering_annotation| 
|**Annotate**|par_annotate_labels|Oligodendrocyte, Oligodendrocyte, Excitatory, Oligodendrocyte, Oligodendrocyte, Microglia, OPC, Oligodendrocyte, Oligodendrocyte , Inhibitory, Atrocyte, Endothlial, Oligodendrocyte, Microglia, Oligodendrocyte, Astrocyte, Oligodendrocyte, Oligodendrocyte, Atrocyte, Oligodendrocyte, Pericyte, Ependymal, OPC, Pericyte, GABA, Inhibitory, Oligodendrocyte, Microglia, Astrocyte, Endothelial, DaN, CADPS2, Microglia 


We can add our annotations using the following code:
```
export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=~/pipeline

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--annotate T
```

The above code produces the following outputs:
```
step7
├── figs7
│   └── annotate
│       ├── clustering_annotation_cluster_annotation.pdf
│       └── clustering_annotation_split_cluster_annotation.pdf
├── info7
│   ├── meta_info_seu_step7.txt
│   ├── sessionInfo_annotate.txt
│   ├── seu_MetaData.txt
│   └── seu_RNA.txt
└── objs7
    └── seu_step7.rds
```

<p align="center">
<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/a2662fa6-227b-4b88-82d0-bc93c1c847db" width="400" height="100"> <br /> 
</p>


**Figure 11. Annotated Seurat object produced in Step 7 of the scRNAbox pipeline.** Uniform Manifold Approximation and Projections (UMAP) plots showing the final cluster annotations for the midbrain dataset. Abbreviations: astro, astrocytes; DaN, dopaminergic neurons; endo, endothelial cells; excit, excitatory neurons; GABA, GABAergic neurons; inhib, inhibitory neurons; mg, microglia; olig, oligodendrocytes; opc, oligodendrocyte precursor cells; peri, pericytes.

- - - -

### Step 8: Differential gene expression  analysis
In Step 8, we are going to perform differential gene expression (DGE) analysis between PD and controls. Prior to computing DGE, we are going to add additional metadata to the Seurat object.
- - - -
#### Add metadata
We are going to add metadata to the Seurat object using a csv file containing the additional metadata:

|Sample_ID|DiseaseStatus|
|:--|:--|
|Parkinson1| PD| 
|Parkinson2| PD| 
|Parkinson3| PD| 
|Parkinson4| PD| 
|Parkinson5| PD| 
|Control1| HC| 
|Control2| HC| 
|Control3| HC| 
|Control4| HC| 
|Control5| HC| 
|Control6| HC| 

Note: This is csv file is available [here](https://github.com/neurobioinfo/scrnabox/blob/main/tutorial/Midbrain_dataset_example_files/metadata.csv).

To add metadata, we set the following execution parameters for Step 8 (`step8_par.txt`):

|DGE method|Parameter|Value|
|:--|:--|:--|
|**General**|par_save_RNA| Yes| 
|**General**|par_save_metadata| Yes|
|**General**|par_seurat_object| NULL |
|**Add metadata**|par_merge_meta|Sample_ID|
|**Add metadata**|par_metadata|~/scrnabox/tutorial/Midbrain_dataset_example_files/metadata.csv|


We can add metadata to the Seurat object using the following code:

```
export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=~/pipeline

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--addmeta T
```

The above code produces the following outputs:
```
step8
├── info8
│   ├── sessionInfo_add_metadata.txt
│   ├── seu_MetaData.txt
│   └── seu_RNA.txt
└── objs8
    └── seu_step8.rds
```
- - - -

Next, we can perform DGE analysis. ScRNAbox can compute DGE between two condition conditions (in our case PD vs control) using **all cells** or **cell type groups**. Furthermore, scRNAbox provides two frameworks for computing DGE:

**1) Cell-based DGE**<br />
Cells are used as replicates and DGE is computed using the Seurat _FindMarkers_ ([Macosko et al. 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4481139/)). While _FindMarkers_ supports several statistical frameworks to compute DGE, we set the default method in our implementation to MAST, which is tailored for scRNAseq data ([Finak et al. 2015](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0844-5))<br />

**2) Sample-based DGE**<br />
Samples are used as replicates by applying a pseudo-bulk analysis. The Seurat _AggregateExpression_ function is used to compute the sum of RNA counts for each gene across all cells from a particular sample ([Cao et al. 2022](https://academic.oup.com/nar/article/50/21/e121/6709246)). The DESq2 statistical framework is then used to compute DGE between conditions using the aggregated counts. ([Love et al. 2014](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8))

- - - -
#### Cell-based DGE using all cells
To perform **cell-based DGE using all cells** we must begin by preparing the contrast matrix (`step8_contrast_cell_based_all_cells.txt`):

```
contrast_name meta_data_variable group1 group2
HCvPD DiseaseStatus HC PD
```
To perform **cell-based DGE using all cells**, we set the following execution parameters for Step 8 (`step8_par.txt`):

|DGE method|Parameter|Value|
|:--|:--|:--|
|**General**|par_save_RNA| No|
|**General**|par_save_metadata| No|
|**General**|par_seurat_object| NULL |
|**Cell-based DGE with all cells**|par_run_cell_based_all_cells|Yes|
|**Cell-based DGE**|par_statistical_method|MAST|

We can perform **cell-based DGE using all cells** using the following code:

```
export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=~/pipeline

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--rundge T
```

The above code produces the following outputs:

```
step8
└── Cell_based_all_cells
    └── HCvPD
        ├── HCvPD_DEG.csv
        └── HCvPD_volcano_plot.pdf
```

- - - -
#### Cell-based DGE using cell type groups 
To perform **cell-based DGE using cell type groups** we must begin by preparing the contrast matrix (`step8_contrast_sample_based_celltype_groups.txt`):

```
contrast_name meta_data_celltype cell_type meta_data_variable group1 group2
OligocendrocytesPDvHC clustering_annotation Oligodendrocyte DiseaseStatus HC PD
MicrogliaPDvHC clustering_annotation Microglia DiseaseStatus HC PD
ExcitatoryPDvHC clustering_annotation Excitatory DiseaseStatus HC PD
InhibitoryPDvHC clustering_annotation Inhibitory DiseaseStatus HC PD
GABAPDvHC clustering_annotation GABA DiseaseStatus HC PD
OPCPDvHC clustering_annotation OPC DiseaseStatus HC PD
AstrocytesPDvHC clustering_annotation Atrocyte DiseaseStatus HC PD
EndothelialPDvHC clustering_annotation Endothelial DiseaseStatus HC PD
PericytesPDvHC clustering_annotation Pericyte DiseaseStatus HC PD
EpendymalPDvHC clustering_annotation Ependymal DiseaseStatus HC PD
DaNPDvHC clustering_annotation DaN DiseaseStatus HC PD
```
To perform **cell-based DGE using all cells**, we set the following execution parameters for Step 8 (`step8_par.txt`):

|DGE method|Parameter|Value|
|:--|:--|:--|
|**General**|par_save_RNA| No|
|**General**|par_save_metadata| No|
|**General**|par_seurat_object| NULL |
|**Cell-based DGE with cell type groups**|par_run_cell_based_cell_type_groups|Yes|
|**Cell-based DGE**|par_statistical_method|MAST|

We can perform **cell-based DGE using all cells** using the following code:

```
export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=~/pipeline

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--rundge T
```

The above code produces the following outputs:

```
step8
└── Cell_based_celltype_groups
    ├── AstrocytesPDvHC
    │   ├── AstrocytesPDvHC_DEG.csv
    │   └── AstrocytesPDvHC_volcano_plot.pdf
    ├── DaNPDvHC
    │   ├── DaNPDvHC_DEG.csv
    │   └── DaNPDvHC_volcano_plot.pdf
    ├── EndothelialPDvHC
    │   ├── EndothelialPDvHC_DEG.csv
    │   └── EndothelialPDvHC_volcano_plot.pdf
    ├── EpendymalPDvHC
    │   ├── EpendymalPDvHC_DEG.csv
    │   └── EpendymalPDvHC_volcano_plot.pdf
    ├── ExcitatoryPDvHC
    │   ├── ExcitatoryPDvHC_DEG.csv
    │   └── ExcitatoryPDvHC_volcano_plot.pdf
    ├── GABAPDvHC
    │   ├── GABAPDvHC_DEG.csv
    │   └── GABAPDvHC_volcano_plot.pdf
    ├── InhibitoryPDvHC
    │   ├── InhibitoryPDvHC_DEG.csv
    │   └── InhibitoryPDvHC_volcano_plot.pdf
    ├── MicrogliaPDvHC
    │   ├── MicrogliaPDvHC_DEG.csv
    │   └── MicrogliaPDvHC_volcano_plot.pdf
    ├── OligocendrocytesPDvHC
    │   ├── OligocendrocytesPDvHC_DEG.csv
    │   └── OligocendrocytesPDvHC_volcano_plot.pdf
    ├── OPCPDvHC
    │   ├── OPCPDvHC_DEG.csv
    │   └── OPCPDvHC_volcano_plot.pdf
    └── PericytesPDvHC
        ├── PericytesPDvHC_DEG.csv
        └── PericytesPDvHC_volcano_plot.pdf

```
- - - -

#### Sample-based DGE using all cells
To perform **sample-based DGE using all cells** we must begin by preparing the contrast matrix (`step8_contrast_sample_based_all_cells.txt`):

```
ContrastName MainContrast SampleID
PDvControl DiseaseStatus Sample_ID
```
To perform **sample-based DGE using all cells**, we set the following execution parameters for Step 8 (`step8_par.txt`):

|DGE method|Parameter|Value|
|:--|:--|:--|
|**General**|par_save_RNA| No|
|**General**|par_save_metadata| No|
|**General**|par_seurat_object| NULL |
|**Sample-based DGE with all cells**|par_run_sample_based_all_cells|Yes|

We can perform **sample-based DGE using all cells** using the following code:

```
export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=~/pipeline

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--rundge T
```

The above code produces the following outputs:

```
step8
└── Sample_based_all_cells
    └── PDvControl
        ├── Aggregated_expression_summary.csv
        ├── DGE_AllCellsMainContrast HC vs PD.csv
        ├── DGE_AllCellsMainContrast HC vs PD.pdf
        └── SampleBased_DGEsummarytable.csv
```

- - - -

#### Sample-based DGE using cell type groups 
To perform **sample-based DGE using cell type groups** we must begin by preparing the contrast matrix (`step8_contrast_sample_based_celltype_groups.txt`):

```
ContrastName CellType MainContrast SampleID
PDvControlbulk clustering_annotation DiseaseStatus Sample_ID
```
To perform **sample-based DGE using cell type groups**, we set the following execution parameters for Step 8 (`step8_par.txt`):

|DGE method|Parameter|Value|
|:--|:--|:--|
|**General**|par_save_RNA| No|
|**General**|par_save_metadata| No|
|**General**|par_seurat_object| NULL |
|**Sample-based DGE with cell type groups**|par_run_sample_based_cell_type_groups|Yes|


We can perform **sample-based DGE using cell type groups** using the following code:

```
export SCRNABOX_HOME=~/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=~/pipeline

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--rundge T
```

The above code produces the following outputs:

```
step8
└── sample_based_celltype_groups
    └── PDvControlbulk
        ├── Aggregated_expression_summary.csv
        ├── figs
        │   ├── DGE_AstrocytesMainContrast HC vs PD.pdf
        │   ├── DGE_DaNMainContrast HC vs PD.pdf
        │   ├── DGE_EndothelialMainContrast HC vs PD.pdf
        │   ├── DGE_EpendymalMainContrast HC vs PD.pdf
        │   ├── DGE_ExcitatoryMainContrast HC vs PD.pdf
        │   ├── DGE_GABAMainContrast HC vs PD.pdf
        │   ├── DGE_InhibitoryMainContrast HC vs PD.pdf
        │   ├── DGE_MicrogliaMainContrast HC vs PD.pdf
        │   ├── DGE_OligodendrocytesMainContrast HC vs PD.pdf
        │   ├── DGE_OPCMainContrast HC vs PD.pdf
        │   └── DGE_PericytesMainContrast HC vs PD.pdf
        ├── info
        │   ├── DGE_AstrocytesMainContrast HC vs PD.csv
        │   ├── DGE_DaNMainContrast HC vs PD.csv
        │   ├── DGE_EndothelialMainContrast HC vs PD.csv
        │   ├── DGE_EpendymalMainContrast HC vs PD.csv
        │   ├── DGE_ExcitatoryMainContrast HC vs PD.csv
        │   ├── DGE_GABAMainContrast HC vs PD.csv
        │   ├── DGE_InhibitoryMainContrast HC vs PD.csv
        │   ├── DGE_MicrogliaMainContrast HC vs PD.csv
        │   ├── DGE_OligodendrocytesMainContrast HC vs PD.csv
        │   ├── DGE_OPCMainContrast HC vs PD.csv
        │   └── DGE_PericytesMainContrast HC vs PD.csv
        └── SampleBased_DGEsummarytable.csv
```

- - - -

<p align="center">
<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/ac770707-560e-41b5-be74-a8771d97a89e" width="800" height="800"> <br /> 
</p>

**Figure 12. Figures produced by Step 8 of the scRNAbox pipeline.** Step 8 of the scRNAbox pipeline allows users to compute DGE between groups by two methods: 1) using cells as replicates with MAST (cell-based) and 2) using samples as replicates with DESeq2 (sample-based). Volcano plots are produced for each user-defined DGE contrast. For the midbrain dataset, we leveraged both methods to compute DGE between Parkinson's disease (PD) subjects and controls across **A)** all cells, **B)** astrocytes, **C)** dopaminergic neurons (DaN), **D)** endothelial cells, **E)** Ependymal cells, **F)** excitatory neurons, **G)** GABAergic neurons, **H)** Inhibitory neurons, **I)** microglia, **J)** oligodendrocytes, **K)** oligodendrocyte precursor cells (OPC), and **H)** Pericytes. 
- - - -

### Analysis of differentially expressed genes 

For the code used to perform downstream analysis of the differentially expressed genes presented in our [pre-print manuscript]() see [here](DEG.md).  

- - - -

### Publication-ready figures
The code used to produce the publication-ready figures used in our [pre-print manuscript]() is avaliable here [here]().  

- - - -

### Job Configurations
The following job configurations were used for our analysis of the midbrain dataset. Job Configurations can be modified for each Analytical Step in the `scrnabox_config.ini` file in `~/pipeline/job_info/configs`

|Step |THREADS_ARRAY|MEM_ARRAY|WALLTIME_ARRAY|
|:--|:--|:--|:--|
|Step2|4|40g|00-05:00|
|Step3|4|40g|00-05:00|
|Step4|4|40g|00-05:00|
|Step5|4|100g|00-05:00|
|Step6|4|100g|00-05:00|
|Step7 MarkerGSEA|4|40g|00-05:00|
|Step7 KnownMarkers|4|40g|00-02:00|
|Step7 ReferenceAnnotation|4|200g|00-12:00|
|Step7 Annotate|4|40g|00-01:00|
|Step8 AddMeta|4|40g|00-02:00|
|Step8 RunDGE|4|40g|00-12:00|




