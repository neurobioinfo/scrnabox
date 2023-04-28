---
layout: post
title: Steps of analyzing HTO via scrnabox pipeline
description: A short introduction to  Hashtag oligonucleotide analyzing using scrnabox pipeline
date: 2022-03-07
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
`scrnabox.slurm` is an open-source pipeline for scRNA analysis that includes a job scheduler for HPC system.  

### Setup
In order to run the pipeline, first create a folder to do the analysis and export the pipeline
```
mkdir -p  ~/scratch/des
export SCRNABOX_HOME=~/scrnabox.slurm
export SCRNABOX_PWD=~/scratch/des
```

Once its 'SCRNABOX_PWD' is defined, you need to create a folder entitled `samples_info` and write samples's `library.csv` and `features_ref.csv`. Then run the following code to setup pipeline, 

```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 0 
```

The pipeline creates few files\folders under `${SCRNABOX_PWD}`: `./job_output/configs/scrnabox.config.ini` (include the configure arguments),  `./job_output/expected.done.files.txt` (recorder the done steps), `./job_output/logs` (submitted job would save under this folder), `./job_output/parameters/` (include the arguments and parameter that would use in running job, you can change them). 

**To-do**  
-We must make it clear how to prepare library.csv and features_ref.csv

### Step 1: cellranger
This step runs cellranger and save the results under `${SCRNABOX_PWD}/step1`, which generates a count matrix. Since cellranger run UI as well, run this step in a `screen`. 
```
screen -S run_scrnabox_HTO
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 1
```

### Step 2: Seurat object 
This step creates the seurat's objects, an standard object for data generated using the 10x Genomics platform, and save the results under `${SCRNABOX_PWD}/step2`
```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 2
```

**To-do**  
-For this step we may want to include an option to output a .csv unfiltered expression matrix


### Step 3: QC and filter
This step run QC and save the results under `${SCRNABOX_PWD}/step3`. The following code filters the data with these criteria: `nFeature_RNA > 300 & nCount_RNA < 6500 & percent.mt < 25`.  
- nFeatures_RNA is the number of unique RNA transcripts for each cell.  If less than 300 we remove these cells as they might be debris or dead cells.  `--nFRNAL` and `--nFRNAU` are the upper and lower thresholds for nFeatures_RNA, respectively.
- Sometimes cells with too many RNA transcripts are dublexs.  It is better to us nCount_RNA to remove dublets. `nCRNAL` and `nCRNAU` are  the upper and lower threshold  for nCount_RNA, respectively. 
- Cells with a high amount of mitochondrial transcript compared to total RNA transcripts might be dead or dying and can add noise to the data making a clustering performance poor. We remove cells setting a default threshold of 25% (which is very high), `--pmtU` is the upper threshold  for the amount of mitochondrial transcript. 
```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 3 \
--nFRNAL 300 \
--nCRNAU 6500 \
--pmtU 25
```
### Step 4: Demuplixing
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

### Step 5: Integration 
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
This step runs the Differetial gene expression (DEG); first add the labels obtained from Step 7 to `/job_output/parameters/step8_ clus_label.txt`. 
 
#### DGEList
This step creates a DGEListobject from a table of counts obtained from seurate objects. It need alot of RAM, we suggest 3*size(seu_int_clu.rds)
```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--dgelist T
```

#### DGE contrasts
In this step, one can run the contrast on clustered result, which can be done on genotype and genotype-cell. 

###### Genotype 
There is a file ${SCRNABOX_PWD}/job_output/parameters/step8_contrast_main.txt, with columns of cont_name,control,ex_control,all, you can write the genotype contrast here, then select `--genotype T` to run the genotype contrast. 
```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--genotype T
```

###### Genotype-cell
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

Note: If you have many contrasts, it is better to split them and submit batch jobs.

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

