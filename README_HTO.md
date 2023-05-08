---
layout: post
title:  A Guide to Analyzing HTO with the Scrnabox Pipeline
description: A short introduction to  Hashtag oligonucleotide analyzing using scrnabox pipeline
date: 2022-05-08
author: Saeid Amiri
published: true
tags: scRNA HTO
categories: 
comments: false
---
## ScRNA pipeline to run HTO   
## Contents

- [Introduction](#introduction)
  - [Setup](#setup)
  - [Step 1: cellranger](#step-1-cellranger)
  - [Step 2: Seurat object](#step-2-seurat-object)  
  - [Step 3: QC and filter](#step-3-qc-and-filter)
  - [Step 4: demuplixing](#step-4-demuplixing)
  - [Step 5: integration](#step-5-integration)
  - [Step 6: Clustering](#step-6-clustering)   
  - [step 7: Celltype annotation](#step-7-celltype-annotation)    
  - [step 8: DEG contrast](#step-8-DEG-contrast)     
- [Integrating seurat objects](#integrating-seurat-objects)  
- [References](#references)

## Introduction 
This guide provides a brief introduction to analyzing hashtag oligonucleotide (HTO) data using the Scrnabox pipeline, `scrnabox.slurm` is an open-source pipeline for scRNA analysis that includes a job scheduler for HPC system. It outlines the steps involved in processing and analyzing HTO data, including quality control, cell filtering, clustering, and contrast analysis. By following these steps, researchers can gain insights into gene expression patterns in single cells and understand the underlying cellular heterogeneity in their samples.
 
### Setup
Before running the pipeline, create a dedicated folder for the analysis and export the pipeline to this folder
```
mkdir -p  ~/scratch/des
export SCRNABOX_HOME=~/scrnabox.slurm
export SCRNABOX_PWD=~/scratch/des
```

After defining the 'SCRNABOX_PWD' variable, create a folder named samples_info and prepare two files - library.csv and features_ref.csv - containing necessary information about the samples. An example format for these files can be found at ![link](https://github.com/neurobioinfo/scrnabox/tree/main/test_code/LaunchSampleHTO).

```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 0 
```

When the pipeline is executed, it generates several files and folders under ${SCRNABOX_PWD}, including ./job_output/configs/scrnabox.config.ini (which contains the configuration arguments), ./job_output/expected.done.files.txt (which records the completed steps), ./job_output/logs (which contains logs for the submitted jobs), and ./job_output/parameters/ (which includes the arguments and parameters used in running the job and can be modified if necessary)."
 
To view a brief overview of the pipeline, run the following code
```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh --help 
```


### Step 1: cellranger
In this step, Cell Ranger is executed, and the resulting output is saved under ${SCRNABOX_PWD}/step1, which generates a count matrix. As Cell Ranger runs a user interface, it is recommended to run this step in a 'screen'.

```
screen -S run_scrnabox_HTO
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 1
```

### Step 2: Seurat object 
This step involves creating Seurat objects, which are a standard format for data generated using the 10x Genomics platform. The resulting objects are saved under 
`${SCRNABOX_PWD}/step2

```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 2
```

### Step 3: QC and filter
In this step, quality control is performed on the data, and the resulting output is saved under `${SCRNABOX_PWD}/step3`. The data is filtered based on the following criteria: `nFeature_RNA > 1000 & nCount_RNA < 65000 & mitochondria_percent < 25`.  
```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 3 \
--nFeature_RNA_L 1000 \
--nCount_RNA_U 65000 \
--mitochondria_percent_U 25
```
The parameters of this step are: 
- nFeatures_RNA: this parameter presents the number of unique RNA transcripts for each cell.  If the value is less than 1000, it is considered as debris or dead cells and removed.  `--nFeature_RNA_L` and `--nFeature_RNA_U` options set  the lower and upper thresholds for nFeatures_RNA, respectively.
-nCount_RNA: Cells with an abnormally high number of RNA transcripts may be doublets. To remove doublets, the nCount_RNA parameter is used instead. The `--nCount_RNA_L` and `--nCount_RNA_U` options set the lower and upper thresholds for nCount_RNA, respectively.
-mitochondria_percent: High mitochondrial transcript levels relative to total RNA transcripts may indicate dead or dying cells and can affect clustering performance. By default, cells with mitochondrial transcript percentages greater than 25% are removed. The `--mitochondria_percent_U` option sets the upper threshold for mitochondrial transcript percentage.


### Step 4: Demultiplexing
If you are using hashtag, you need to choose the right label (for the hashtags), you can get the hashtag labels by running the following code 
```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 4 \
--msd T 
```

Add the old label and its new corresponding label in '${SCRNABOX_PWD}/job_output/parameters/step4_par.txt'. Run the following to run the demuplixing  
```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 4 
```

This step can be used to remove the 'Doublet'; the default is to remove the doublet, if you want to keep them, just change 'yes' to 'no' in '${SCRNABOX_PWD}/job_output/parameters/step4_par.txt'. 

<!--
**To-do**  
1. We need to justify why we are using the DoubletFinder tool as opposed to other doublet detecting tools  
-we will have to decide as a group how we want to proceed with doublet removal. We can i) justify DoubletFinder ii) Determine which doublet tool is the best from the literature iii) Incorporate multiple doublet tools into the pipeline and let the user choose which tool they'd like to use.  
2. For this step we may want to include an option to output a .csv filtered expression matrix that now contains the original hashtag label.    

|  | Gene X | Gene y |
| - | - | - |
| Cell 1 | Counts | Counts |
| Cell 2 | Counts | Counts |
| Cell n | Counts | Counts |
-->

### Step 5: Integration 
```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 5 
```
<!--
**To-do**  
1. One of the outputs is **mata_info.csv** we should change this to **meta_info.csv**  
2. For this step we may want to include an option to output a .csv expression matrix for the integrated assay  

|  | Run ID | Gene X | Gene y |
| - | - | - | -|
| Cell 1 | 1 | Counts | Counts |
| Cell 2 | 1 | Counts | Counts |
| Cell n | 2 | Counts | Counts |
-->

### Step 6: Clustering 
```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 6 
```

<!--
**To-do**  
1. One of the outputs is **mata_info.csv** we should change this to **meta_info.csv**  
2. For this step we may want to include an option to output a .csv expression matrix for the integrated assay that shows cluster information 

|  | Run ID | Phenotype | Cluster | Gene X | Gene y |
| - | - | - | - | - | - |
| Cell 1 | 1 | Case | 1 | Counts | Counts |
| Cell 2 | 1 | Case | 1 | Counts | Counts |
| Cell n | 2 | Control | 1 | Counts | Counts | 
-->

### step 7: Cluster annotation
In This step, you should find the cluster annotation to use in the Step 8. 
#### Marker 
Finds markers (differentially expressed genes) for each of cluster
```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--marker T
```

#### FindTransferAnchors
Find a set of anchors between a reference and query object and add it to query object `predictions`
```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--fta T
```

#### EnrichR 
If your HPC allows the access to internet under batch submission, run the following codes
```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--enrich T
```

Otherwise, run the enrichment directly:
```
level_cluster='integrated_snn_res.0.7'
PWD='/SCRNABOX_PWD/step7/annot/'
PSUE='/SCRNABOX_PWD/step6/objs/seu_step6.rds'
top_sel=5
db <- c('Descartes_Cell_Types_and_Tissue_2021','CellMarker_Augmented_2021','Azimuth_Cell_Types_2021')
scrnaboxR::annotation(level_cluster,PWD,PSUE,top_sel,db)
```

To run enrichment on your PC, first copy it to your PC, 
```
scp -r usrid@beluga.computecanada.ca:${SCRNABOX_PWD}/step6 ~/Desktop/annot/
```

Then run the following codes
```
level_cluster='integrated_snn_res.0.7'
PWD='~/Desktop/annot/'
PSUE='~/Desktop/annot/step6/objs/seu_step6.rds'
top_sel=5
db <- c('Descartes_Cell_Types_and_Tissue_2021','CellMarker_Augmented_2021','Azimuth_Cell_Types_2021')
scrnaboxR::annotation(level_cluster,PWD,PSUE,top_sel,db)
```

### step 8: DGE contrast
This step involves running the differential gene expression (DEG) analysis. To do this, you will first need to add the labels obtained from Step 7 to `/job_info/parameters/step8_ clus_label.txt`. These labels typically represent cell type or condition information and should be separated by commas. Once the labels have been added, you can use a variety of methods to perform DEG analysis, such as the popular `edgeR` or `DESeq2` packages, this pipeline uses `edgeR`.

 
#### DGEList
This step involves creating a DGEListobject from a table of counts obtained from seurate objects.it is recommended to allocate at least 3 times the size of the seu_int_clu.rds file in RAM, 3*size(seu_int_clu.rds).

```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--dgelist T
```

#### DGE contrasts
At this stage, it is possible to perform a contrast analysis on the clustered results, either at the genotype level or at the genotype-cell level

###### Genotype 
You can find a file called `${SCRNABOX_PWD}/job_info/parameters/step8_contrast_main.txt`, which contains columns for `cont_name`, `control`, `ex_control`, and `all`. You can specify genotype contrasts in this file and then use the `--genotype T` option to run the genotype contrast

```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--genotype T
```

###### Genotype-cell
To perform an interaction analysis between cell type and genotype, specify your contrast in the file `${SCRNABOX_PWD}/job_info/parameters/step8_contrast_inte.txt`. To run Step 8 using the interaction contrast, execute the following command which use the `-celltype T` option to run the Genotype-cell contrast. 


```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--celltype T
```

You can integrate the contrast directly into the pipeline by calling it, 
```
CONTINT=~/des/step7_contrast_inte.txt
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--celltype T \
--cont ${CONTINT}
```

```
CONTMAIN=~/des/step7_contrast_main.txt
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--genotype T \
--cont ${CONTMAIN}
```

Note: If you have a large number of contrasts to run, it may be more efficient to split them up and submit batch jobs instead.

## Integrating seurat objects
To combine different seurat objects, you can run the following codes. 
```
LISTOFSEU=~/list.txt

sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps  integrate \
--seulist ${LISTOFSEU} \
```

LISTOFSEU includes the path of seurate files, put them in different lines. By default, the pipeline standardize the seurat objects before integrating, you can change the default in `${SCRNABOX_PWD}/job_output/parameters/stepint_par.txt`.

------------------------------------
<!-- This is commented out.
## Run multiple steps
To run multiple steps just specify the range of steps, separate using hyphen. 
```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 2-4 \
--nFRNAl 300 \
--nFRNAu 6500 \
--pmt 25


sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 5-6 


sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 3-6 \
--nFRNAl 300 \
--nFRNAu 6500 \
--pmt 25
```

To run all steps, write 'ALL' infront `steps`.   
```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps ALL \
```

## References  -->

