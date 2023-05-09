---
layout: post
title: Steps of scRNA pipeline to run Non-hashtag data
description: A short introduction to scRNA pipeline to run Non-hashtag data
date: 2023-05-09
author: Saeid Amiri
published: true
tags: scRNA 
categories: 
comments: false
---
## ScRNA pipeline to run standard scRNA
## Contents
- [Introduction](#introduction)
  - [Setup](#setup)
  - [Step 1: cellranger](#step-1-cellranger)
  - [Step 2: Seurat object](#step-2-seurat-object)  
  - [Step 3: QC and filter](#step-3-qc-and-filter)
  - [Step 4: demuplixing\non-demuplexing](#step-4-demuplixing\nondemuplexing)
  - [Step 5: integration](#step-5-integration)
  - [Step 6: Clustering](#step-6-clustering)   
  - [step 7: Celltype annotation](#step-7-celltype-annotation)    
  - [step 8: DEG contrast](#step-8-DEG-contrast)     
- [Integrating seurat objects](#integrating-seurat-objects)  
- [References](#references)

## Introduction 
This guide provides a brief introduction to analyzing standard data using the Scrnabox pipeline, scrnabox.slurm is an open-source pipeline for scRNA analysis that includes a job scheduler for HPC system. It outlines the steps involved in processing and analyzing HTO data, including quality control, cell filtering, clustering, and contrast analysis. By following these steps, researchers can gain insights into gene expression patterns in single cells and understand the underlying cellular heterogeneity in their samples.
<br />
<br />
<kbd>
![Steps of Standard scRNA-seq ](https://raw.githubusercontent.com/neurobioinfo/scrnabox/main/figs/scrna.png)
</kbd>

### Setup
Before running the pipeline, create a dedicated folder for the analysis and export the pipeline to this folder
```
mkdir -p  ~/scratch/des
export SCRNABOX_HOME=~/scrnabox.slurm
export SCRNABOX_PWD=~/scratch/des
```
After defining the 'SCRNABOX_PWD' variable, create a folder named samples_info and prepare two files - library.csv and features_ref.csv - containing necessary information about the samples. An example format for these files can be found at [link](https://github.com/neurobioinfo/scrnabox/tree/main/test_code/LaunchSampleHTO).  
Then run the following code to setup pipeline:
```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 0 \
--method SCRNA
```

When the pipeline setup is executed, it generates several files and folders under ${SCRNABOX_PWD}, including ./job_output/configs/scrnabox.config.ini (which contains the configuration arguments), ./job_output/expected.done.files.txt (which records the completed steps), ./job_output/logs (which contains logs for the submitted jobs), and ./job_output/parameters/ (which includes the arguments and parameters used in running the job and can be modified if necessary)."  Since you are using hashtags, you need to choose `--method SCRNA`. 


### Step 1: cellranger
In this step, Cell Ranger is executed, and the resulting output is saved under ${SCRNABOX_PWD}/step1, which generates a count matrix. As Cell Ranger runs a user interface, it is recommended to run this step in a 'screen'.
```
screen -S run_scrnabox
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
 - nCount_RNA: Cells with an abnormally high number of RNA transcripts may be doublets. To remove doublets, the nCount_RNA parameter is used instead. The `--nCount_RNA_L` and `--nCount_RNA_U` options set the lower and upper thresholds for nCount_RNA, respectively.
 - mitochondria_percent: High mitochondrial transcript levels relative to total RNA transcripts may indicate dead or dying cells and can affect clustering performance. By default, cells with mitochondrial transcript percentages greater than 25% are removed. The `--mitochondria_percent_U` option sets the upper threshold for mitochondrial transcript percentage.

### Step 4: Double
This step can be used to remove the 'Doublet'. By default, the pipeline removes the doublet, if you want to keep them, just change 'yes' to 'no' in '${SCRNABOX_PWD}/job_info/parameters/step4_par.txt'. 
```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 4 
```

### Step 5: Integration 
In this step, the pipeline  combines multiple single-cell RNA-seq datasets.
```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 5 
```

### Step 6: Clustering 
In this step, the pipeline runs clustering on the integrated dataset to group cells with similar gene expression patterns together based on a k-nearest neighbor graph. 

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
If Your HPC allows the access to internet under batch submission, run the following codes
```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--enrich T
```

Otherwise, run the enrichment directively:
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
This step creates a DGEListobject from a table of counts obtained from seurate objects. It might need alot of RAM, we suggest 3*size(seu_int_clu.rds)
```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--dgelist T
```

#### DGE contrasts
In this step, one can run the contrast on clustered result, which can be done on genotype and genotype-cell. 

#### Genotype 
There is a file ${SCRNABOX_PWD}/job_output/parameters/step8_contrast_main.txt, with columns of cont_name,control,ex_control,all, you can write the genotype contrast here, then select `--genotype T` to run the genotype contrast. 
```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--genotype T
```

#### Genotype-cell
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

