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
## Application Case 1: Standard scRNAseq Analysis Track of scRNAbox
## Contents

- [Introduction](#introduction)
- [Downloading the public dataset](#downloading-the-public-dataset)
- [scRNAbox preparation](#scrnabox-preparation)
  - [From step 2](#from-step-2)
  - [From step 3](#from-step-3)  
  - [From step 4](#from-step-4)
  - [From Step 5: integration](#from-step-5-integration)
  - [From Step 6: Clustering](#from-step-6-clustering)   

## Introduction 
This guide illustrates the steps taken for Application Case 1 in our pre-print manuscript. Here, we are using the Standard scRNAseq Analysis Track of scRNAbox to analyze a publicly available scRNAseq dataset produced by [Smajic et al.](https://academic.oup.com/brain/article/145/3/964/6469020). This data set (referred to as the midbrain dataset in the manuscript) describes >41,000 single-nuclei transcriptomes from the post-mortem midbrains of five individuals with Parkinson’s disease (PD) and six controls sequenced separately. 


## Downloading the public dataset
The scRNAseq data produced by Smajic et al. is publicly available in the Gene Expression Omnibus with accession code [GSE157783](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157783). To download the data, we must first install SRAtoolkit (if this is not already installed on your High-Performance Computing (HPC) system). We will create a directory for our raw data and download SRAtoolkit with the following code:

```
mkdir data_download
cd data_download
wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -vxzf sratoolkit.tar.gz
export PATH=$PATH:$PWD/sratoolkit.3.0.5-ubuntu64/bin
```

For more information regarding the SRAtoolkit, please visit the [documentation](https://github.com/ncbi/sra-tools/wiki).

The Sequence Read Archive (SRA) run identifier for each of the 11 samples in the midbrain dataset are:

|Sample|SRR|
|:--|:--|
|PD1|SRR12621862|
|PD2|SRR12621863|
|PD3|SRR12621864|
|PD4|SRR12621865|
|PD5|SRR12621866|
|CTRL1|SRR12621867|
|CTRL2|SRR12621868|
|CTRL3|SRR12621869|
|CTRL4|SRR12621870|
|CTRL5|SRR12621871|
|CTRL6|SRR12621872|

**Note**: If you simply want to test scRNAbox's Standard scRNAseq Analysis Track, it may be best to only incorportate a subset of samples in a test run, as using all 11 samples will take substantially longer. In this case, we suggest including at least one PD sample and one control to facilitate differential gene expression (DGE) contrasts in Step 8. 

To download the FASTQ files for all 11 samples, run the following code. Please note that this may take a very long time. 

```
export PATH=$PATH:$PWD/sratoolkit.3.0.5-ubuntu64/bin
module load StdEnv/2020 gcc/9.3.0
module load sra-toolkit/3.0.0 

#PD1
prefetch SRR12621862 
fasterq-dump SRR12621862 

#PD2
prefetch SRR12621863 
fasterq-dump SRR12621863  

#PD3
prefetch SRR12621864 
fasterq-dump SRR12621864  

#PD4
prefetch SRR12621865 
fasterq-dump SRR12621865  

#PD5
prefetch SRR12621866 
fasterq-dump SRR12621866 

#CTRL1
prefetch SRR12621867 
fasterq-dump SRR12621867 

#CTRL2
prefetch SRR12621868 
fasterq-dump SRR12621868  

#CTRL3
prefetch SRR12621869 
fasterq-dump SRR12621869  

#CTRL4
prefetch SRR12621870 
fasterq-dump SRR12621870  

#CTRL5
prefetch SRR12621871 
fasterq-dump SRR12621871  

#CTRL6
prefetch SRR12621872
fasterq-dump SRR12621872 
```
If the FASTQ files for all 11 samples have been downloaded properly, the `data_download` folder should contain the following:
```
├── SRR12621862
│   └── SRR12621862.sra
├── SRR12621862_1.fastq
├── SRR12621862_2.fastq
├── SRR12621863
│   └── SRR12621863.sra
├── SRR12621863_1.fastq
├── SRR12621863_2.fastq
├── SRR12621864
│   └── SRR12621864.sra
├── SRR12621864_1.fastq
├── SRR12621864_2.fastq
├── SRR12621865
│   └── SRR12621865.sra
├── SRR12621865_1.fastq
├── SRR12621865_2.fastq
├── SRR12621866
│   └── SRR12621866.sra
├── SRR12621866_1.fastq
├── SRR12621866_2.fastq
├── SRR12621867
│   └── SRR12621867.sra
├── SRR12621867_1.fastq
├── SRR12621867_2.fastq
├── SRR12621868
│   └── SRR12621868.sra
├── SRR12621868_1.fastq
├── SRR12621868_2.fastq
├── SRR12621869
│   └── SRR12621869.sra
├── SRR12621869_1.fastq
├── SRR12621869_2.fastq
├── SRR12621870
│   └── SRR12621870.sra
├── SRR12621870_1.fastq
├── SRR12621870_2.fastq
├── SRR12621871
│   └── SRR12621871.sra
├── SRR12621871_1.fastq
├── SRR12621871_2.fastq
├── SRR12621872
│   └── SRR12621872.sra
├── SRR12621872_1.fastq
└── SRR12621872_2.fastq
```


Next, we will rename the FASTQ files according to the CellRanger nomenclature and transfer the FASTQ files to a folder named `fastqs`. For more information regarding the nomeclature required by the CellRanger _counts_ pipeline, please visit CellRanger's [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/fastq-input).

**Note**: The `fastqs` folder should only contain FASTQ files for the experiment.

```
mkdir fastqs

#PD1
cp ~/data_download/SRR12621862_1.fastq ~/fastqs/PD1_S1_L001_R1_001.fastq
cp ~/data_download/SRR12621862_2.fastq ~/fastqs/PD1_S1_L001_R2_001.fastq

#PD2
cp ~/data_download/SRR12621863_1.fastq ~/fastqs/PD2_S1_L001_R1_001.fastq
cp ~/data_download/SRR12621863_2.fastq ~/fastqs/PD2_S1_L001_R2_001.fastq

#PD3
cp ~/data_download/SRR12621864_1.fastq ~/fastqs/PD3_S1_L001_R1_001.fastq
cp ~/data_download/SRR12621864_2.fastq ~/fastqs/PD3_S1_L001_R2_001.fastq

#PD4
cp ~/data_download/SRR12621865_1.fastq ~/fastqs/PD4_S1_L001_R1_001.fastq
cp ~/data_download/SRR12621865_2.fastq ~/fastqs/PD4_S1_L001_R2_001.fastq

#PD5
cp ~/data_download/SRR12621866_1.fastq ~/fastqs/PD5_S1_L001_R1_001.fastq
cp ~/data_download/SRR12621866_2.fastq ~/fastqs/PD5_S1_L001_R2_001.fastq

#Ctrl1
cp ~/data_download/SRR12621867_1.fastq ~/fastqs/CTRL1_S1_L001_R1_001.fastq
cp ~/data_download/SRR12621867_2.fastq ~/fastqs/CTRL1_S1_L001_R2_001.fastq

#Ctrl2
cp ~/data_download/SRR12621868_1.fastq ~/fastqs/CTRL2_S1_L001_R1_001.fastq
cp ~/data_download/SRR12621868_2.fastq ~/fastqs/CTRL2_S1_L001_R2_001.fastq

#Ctrl3
cp ~/data_download/SRR12621869_1.fastq ~/fastqs/CTRL3_S1_L001_R1_001.fastq
cp ~/data_download/SRR12621869_2.fastq ~/fastqs/CTRL3_S1_L001_R2_001.fastq

#Ctrl4
cp ~/data_download/SRR12621870_1.fastq ~/fastqs/CTRL4_S1_L001_R1_001.fastq
cp ~/data_download/SRR12621870_2.fastq ~/fastqs/CTRL4_S1_L001_R2_001.fastq

#Ctrl5
cp ~/data_download/SRR12621871_1.fastq ~/fastqs/CTRL5_S1_L001_R1_001.fastq
cp ~/data_download/SRR12621871_2.fastq ~/fastqs/CTRL5_S1_L001_R2_001.fastq

#Ctrl6
cp ~/data_download/SRR12621872_1.fastq ~/fastqs/CTRL6_S1_L001_R1_001.fastq
cp ~/data_download/SRR12621872_2.fastq ~/fastqs/CTRL6_S1_L001_R2_001.fastq
```
If the above steps were conducted properly, the `fastqs` folder should contain the following files:
```
├── CTRL1_S1_L001_R1_001.fastq
├── CTRL1_S1_L001_R2_001.fastq
├── CTRL2_S1_L001_R1_001.fastq
├── CTRL2_S1_L001_R2_001.fastq
├── CTRL3_S1_L001_R1_001.fastq
├── CTRL3_S1_L001_R2_001.fastq
├── CTRL4_S1_L001_R1_001.fastq
├── CTRL4_S1_L001_R2_001.fastq
├── CTRL5_S1_L001_R1_001.fastq
├── CTRL5_S1_L001_R2_001.fastq
├── CTRL6_S1_L001_R1_001.fastq
├── CTRL6_S1_L001_R2_001.fastq
├── PD1_S1_L001_R1_001.fastq
├── PD1_S1_L001_R2_001.fastq
├── PD2_S1_L001_R1_001.fastq
├── PD2_S1_L001_R2_001.fastq
├── PD3_S1_L001_R1_001.fastq
├── PD3_S1_L001_R2_001.fastq
├── PD4_S1_L001_R1_001.fastq
├── PD4_S1_L001_R2_001.fastq
├── PD5_S1_L001_R1_001.fastq
└── PD5_S1_L001_R2_001.fastq
```

## Installation
#### scrnabox.slurm installation
Now that the raw data has been downloaded and organized, we can install the latest version of `scrnabox.slurm` (v0.135):
```
wget https://github.com/neurobioinfo/scrnabox/releases/download/v0.1.35/scrnabox.slurm.zip
unzip scrnabox.slurm.zip
```
For a description of the options for running `scrnabox.slurm` run the following command:
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
#### CellRanger installation 

**Note:** If CellRnager is already installed on your HPC system you may skip the CellRanger installation procedures

For our analysis of the midbrain dataset we used the 10XGenomics GRCh38-3.0.0 reference genome. For more information regarding how to prepare reference genomes for the CellRanger _counts_ pipeline, please see the 10X Genomics [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#grch38_3.0.0).


#### R library and R package installation
Next, we can prepare a common R library where we will load all of the required R packages. If the required R packages are already installed on your HPC system in a common R library, you may skip the following steps. 
First, we will creat an `R` folder and download our desired R version. The analyses presented in our pre-print manuscript were conducted using R v4.2.1

```
mkdir R
cd R

```


## scRNAbox: Standard Analysis Track
#### Step 0: Pipeline initiation
Now that `scrnabox.slurm`, `CellRanger`, `R`, and the Required R packages have been installed, we can proceed to our analysis with the Standard scRNAseq Analysis Track of the scRNAbox pipeline. We will create a `pipeline` folder designated for the analysis and run the pipeline initiation Step using the following code:
```
mkdir pipeline
cd pipeline

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 0 \
--method SCRNA
```
Next, we will navigate to the `scrnabox_config.ini` file in `~/pipeline/job_info/configs` to define the path to the R library (`R_LIB_PATH=`), the version of R (`R_VERSION=`), and the path to CellRanger (`MODULECELLRANGER=`):
```
cd ~/pipeline/job_info/configs
nano scrnabox_config.ini

MODULECELLRANGER=mugqic/cellranger/5.0.1
R_VERSION=4.2.1
R_LIB_PATH=~/R
```

#### Step 1: FASTQ pre-processing
In this Step, we will run the CellRanger _counts_ pipeline to generate feature-barcode expression matrices from the FASTQ files. While it is possible to manually prepare the `library.csv` files for each of the 11 samples in the experiment prior to running Step 1, for this analysis we are going to opt for automated library preparation. For more information regarding the manual prepartion of `library.csv` files please see the Standard scRNAseq [documentation](SCRNA.md) under the **Setup** section. <br /> 
<br /> 
For our analysis of the midbrain dataset we set the following execution parameters for Step 1:

|Parameter|Value|
|:--|:--|
|par_automated_library_prep|yes|
|par_fastq_directory|~/fastqs|
|par_sample_names|PD1, PD2, PD3, PD4, PD5, CTRL1, CTRL2, CTRL3, CTRL4, CTRL5, CTRL6|
|par_rename_samples|Yes|
|par_new_sample_names|Parkinson1, Parkinson2, Parkinson3, Parkinson4, Parkinson5, Control1, Control2, Control3, Control4, Control5, Control6|
|par_paired_end_seq|TRUE|
|REF_DIR_GRCH|~/genome/10xGenomics/refdata-cellranger-GRCh38-3.0.0|
|R1LENGTH|NULL|
|MEMPERCORE|30|

**Note:** For a comprehensive description of the execution parameters for each Analytical Step see the [Reference](reference.md) section of the scRNAbox documentation. 

#### Step 2: Ambient RNA removal and create Seurat object




