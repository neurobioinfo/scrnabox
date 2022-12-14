---
layout: post
title: Steps of scRNA hashtags pipline
description: A short introduction to scRNA hashtags pipline
date: 2022-12-12
author: Saeid Amiri
published: true
tags: scRNA 
categories: 
comments: false
---
## ScRNA  pipline  
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
scrnabox.svn is an open-source pipeline for scRNA analysis that includes a job scheduler for HPC system.  

### Setup
In order to run the pipeline, first create a folder to do the analysis and export the pipeline

```
mkdir -p  ~/scratch/des
export SCRNABOX_HOME=~/scrnabox.svn
export SCRNABOX_PWD=~/scratch/des
```

Once its 'SCRNABOX_PWD' is defined, you need to create a folder entitled `samples_info` and write samples's `library.csv` and `features_ref.csv`. Then run the following code to setup pipelione, 

```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 0 
```

The pipeline creates files\folders under `${SCRNABOX_PWD}`: `./job_output/configs/scrnabox.config.ini` (include the configure arguments),  `./job_output/expected.done.files.txt` (recorder the done steps), `./job_output/logs` (submitted job would save under this folder), `./job_output/parameters/` (include the arguments and parameter that would use in running job, you can change them). 

### Step 1: cellranger
This step runs cellranger and save the results under `${SCRNABOX_PWD}/step1`. Since cellranger run UI as well, run this step in a `screen`. 
```
screen 
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 1
```

### Step 2: Seurat object 
This step creates seurat's objects and save results under `${SCRNABOX_PWD}/step2`

```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 2
```

### Step 3: QC and filter
This step run QC and save the results under `${SCRNABOX_PWD}/step3`. The following code filter with these criteria: `nFeature_RNA > 300 & nCount_RNA < 6500 & percent.mt < 25`.  
- nFeatures_RNA is the number of unique RNA transcripts for each cell.  If less than 300 we remove these cells as they might be debris or dead cells
- Sometimes cells with too many RNA transcripts are dublexs.  It is better to us nCount_RNA to remove dublets. 
- Cells with a high amount of mitochondrial transcript compared to total RNA transcripts might be dead or dying and can add noise to the data making a clustering performance poor. We remove cells setting a default threshold of 25% (which is very high)

```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 3 \
--nFRNA 300 \
--nCRNA 6500 \
--pmt 25
```

### Step 4: Demuplixing 
In this step, you need to choose the right label (for the hashtags), you can get the hashtag labels by running the following code 

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

### Step 5: Integration 
This step can be done with\without removing the 'Doublet', 'Negative'; the defalut is to remove them, if you want to keep them, just change 'yes' to 'no' in '${SCRNABOX_PWD}/job_output/parameters/step5_par.txt'. 

```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 5 
```

### Step 6: Clustering 
```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 6 
```

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
Since most of HPC can not connect to the internet when the job is submitted, so you should run enrichR locally; first copy it to your PC, 
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
This step run Differetial gene expression (DEG), first add the labels obtained from Step 7 to `/job_output/parameters/step8_ clus_label.txt`. 
 
A) DGEList
This step creates a DGEListobject from a table of counts obtained from seurate objects. This step might need alot of RAM, we suggest 3*size(seu_int_clu.rds)

```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--dgelist T
```

B) DGE contrasts

In this step, one can run the contrast on clustered result, which can be done on genotype and genotype-cell. 

#### genotype 
There is a file ${SCRNABOX_PWD}/job_output/parameters/step8_contrast_main.txt, with columns of cont_name,control,ex_control,all, you can write the genotype contrast here, then select `--genotype T` to run the genotype contrast. 
```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--genotype T
```

#### genotype-cell
To run interact between celltype and genptype, write your contrast in `${SCRNABOX_PWD}/job_output/parameters/step8_contrast_inte.txt`. To run Step 8 on interact contrast, run the following command. Select `-celltype T` to run the main contrast. 
```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--celltype T
```

You can directly call the contrast to the pipeline, 
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

Note: If you have many contrasts, it is better to split them and submit differently.  

## Integrating seurat objects
In order to combine different serurat objects, you can the following codes. By default, the pipline standardize the seurat objects
before integrating, you can change the default in `${SCRNABOX_PWD}/job_output/parameters/step_addseu_par.txt`. 


```
CONTINT2=/home/samamiri/NB043_dge/comm/Dark_Genome/analysis_DarkGenome6weeks_test/ADDSEU0/seu_int.rds
CONTINT1=/home/samamiri/NB043_dge/comm/Dark_Genome/analysis_DarkGenome6weeks_test/ADDSEU0/seu.int.c.rds
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps  addseu \
--seu1 ${CONTINT1} \
--seu2 ${CONTINT2}
```

The integrated seurat object `seu_inetgrated.rds` would save under working directory. 


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

