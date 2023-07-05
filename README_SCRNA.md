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
## ScRNA pipeline to run standard scRNA
## Contents
- [Introduction](#introduction)
  - [Setup](#setup)
  - [Step 1: cellranger](#step-1-cellranger)
  - [Step 2: Seurat object](#step-2-seurat-object)  
  - [Step 3: QC and filter](#step-3-qc-and-filter)
  - [Step 4: Remove doublet](#step-4-remove-doublet)
  - [Step 5: integration](#step-5-integration)
  - [Step 6: Clustering](#step-6-clustering)   
  - [step 7: Cluster annotation](#step-7-cluster-annotation)    
  - [step 8: DEG contrast](#step-8-dge-contrast)     
- [Integrating seurat objects](#integrating-seurat-objects)  

## Introduction 
This guide provides a brief introduction to analyzing standard data using the Scrnabox pipeline, `scrnabox.slurm` is an open-source pipeline for scRNA analysis that includes a job scheduler for HPC system. It outlines the steps involved in processing and analyzing data, including quality control, cell filtering, clustering, and contrast analysis. By following these steps, researchers can gain insights into gene expression patterns in single cells and understand the underlying cellular heterogeneity in their samples. The following figure illustrates the steps involved in analyzing standard scRNA-seq data with the scrnabox pipeline
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
To obtain a brief guidance on the parameters, please execute the following code.

```
bash $SCRNABOX_HOME/launch_scrnabox.sh -h 
```

After defining the 'SCRNABOX_PWD' variable, create a folder named `samples_info`, then add subfolder for each sample in the `samples_info`, then prepare two files - library.csv and features_ref.csv - containing necessary information about the sample and save then inside the subfolder (note, the pipeline uses the name of subfolder as sample name under `orig.ident` in Seurat object). An example format for these files can be found at [link](https://github.com/neurobioinfo/scrnabox/tree/main/test_code/LaunchSamplescRNA); create a CSV file named library.csv with three columns: `fastq`, `sample`, and `library_type`. In the fastq column, provide the path to the file. In the sample column, write the first of sample name, e.g., write `Sample1GEXD01_MPS12347745_C12_S1_R1_001.fastq.gz` as "Sample1GEXD01_MPS12347745_C12". In the library_type column, specify the type of the library.  
Then run the following code to setup pipeline:
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 0 \
--method SCRNA
```
In the code, `-d ${SCRNABOX_PWD}`, `-steps 0`, and `--method SCRNA` specify certain parameters for the pipelines. The -d ${SCRNABOX_PWD} flag sets the working directory, -steps 0 determines the specific step to be executed (in this case, we choose the setup step), and `--method SCRNA` indicates that the pipeline is running the SCRNA method. When the pipeline setup is executed, it generates several files and folders under ${SCRNABOX_PWD}, including `./job_output/configs/scrnabox.config.ini` (which contains the configuration arguments), `./job_output/expected.done.files.txt` (which records the completed steps), `./job_output/logs` (which contains logs for the submitted jobs), and `./job_output/parameters/` (which includes the arguments and parameters used in running the job and can be modified if necessary)."  Since you are using hashtags, you need to choose `--method SCRNA`. 


### Step 1: cellranger
In this step, Cell Ranger is executed, and the resulting output is saved under ${SCRNABOX_PWD}/step1, which generates a count matrix. As Cell Ranger runs a user interface, it is recommended to run this step in a 'screen'.
```
screen -S run_scrnabox
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 1
```

### Step 2: Seurat object 
This step involves creating Seurat objects, which are a standard format for data generated using the 10x Genomics platform. The resulting objects are saved under  `${SCRNABOX_PWD}/step2`
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 2
```

### Step 3: QC and filter
In this step, quality control is performed on the data, and the resulting output is saved under `${SCRNABOX_PWD}/step3`. The data is filtered based on the following criteria: `nFeature_RNA > 1000 & nCount_RNA < 65000 & mitochondria_percent < 25`.  
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
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

### Step 4: Remove doublet
This step can be used to remove the 'Doublet'. By default, the pipeline removes the doublet, if you want to keep them, just change 'yes' to 'no' in '${SCRNABOX_PWD}/job_info/parameters/step4_par.txt'.   Also you can add the current label and its corresponding new label in the file '${SCRNABOX_PWD}/job_output/parameters/step4_par.txt'. Once you have added the labels, run the following command to run this step.

```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 4 
```

### Step 5: Integration 
In this step, one combines multiple single-cell RNA-seq datasets.
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 5 
```

### Step 6: Clustering 
In this step, the pipeline runs clustering on the integrated dataset to group cells with similar gene expression patterns together based on a shared nearest neighbor network. 

```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 6 
```

### step 7: Cluster annotation
In this step, you will use various methods to identify and annotate cell clusters. This may involve creating t-SNE or UMAP plots to visualize cell clusters, or using marker gene analysis or gene set enrichment analysis to identify cell types. The goal is to determine the cell types or states present in your dataset and assign them to each cluster. The cluster annotation information will be used in Step 8 for downstream analysis and interpretation
#### Marker 
In this step, the pipeline finds differentially expressed genes for each cluster, which can be used to identify cluster-specific markers. This is done using the FindAllMarkers function in Seurat. The function compares gene expression in each cluster to the expression in all other clusters and identifies genes that are differentially expressed with a significant p-value. The output includes the top differentially expressed genes for each cluster and their corresponding p-values and fold changes. These markers can be used for downstream analysis such as cell type identification and functional annotation. The results are saved under ${SCRNABOX_PWD}/step7. This step produces `./step7/info7/top_sel.csv`, `./step7/info7/cluster_just_genes.xlsx`, `./step7/info7/cluster_whole.xlsx`, `./step7/info7/ClusterMarkers.rds`,  `./step7/figs7/heatmap.pdf`, `./step7/figs7/umap.pdf`, `./step7/figs7/umap_splitted.pdf`

```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--marker T
```

#### FindTransferAnchors
In this step, the pipeline uses the `FindTransferAnchors` function in Seurat identifies anchors between a reference and query object and add it to query object `predictions`. This step produces `./step7/objs7/seu_step7.rds`.  

```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--fta T
```

#### EnrichR 
In this step, the pipeline uses the EnrichR tool to perform gene set enrichment analysis on the differentially expressed genes identified in the previous step. EnrichR compares the list of genes against a large collection of gene set libraries, including pathways, gene ontology terms, and transcription factor targets, to identify enriched sets of genes that are related to biological functions, pathways, or processes. The results are saved as pdfs and csvs file.
 
If your HPC allows access to the internet during batch submission, you can run the following codes
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--enrich T
```

Otherwise, run the enrichment directly:
```
level_cluster='integrated_snn_res.0.7'
ClusterMarkers='/SCRNABOX_PWD/step7/info7/ClusterMarkers.rds'
PWD='/SCRNABOX_PWD/step7/annot/'
PSUE='/SCRNABOX_PWD/step6/objs/seu_step6.rds'
top_sel=5
db <- c('Descartes_Cell_Types_and_Tissue_2021','CellMarker_Augmented_2021','Azimuth_Cell_Types_2021')
scrnaboxR::annotation(level_cluster,ClusterMarkers,PWD,PSUE,top_sel,db)
```

To run enrichment on your PC, first copy it to your PC, 
```
scp -r usrid@beluga.computecanada.ca:${SCRNABOX_PWD}/step7 ~/Desktop/annot/
scp -r usrid@beluga.computecanada.ca:${SCRNABOX_PWD}/step6 ~/Desktop/annot/
```

Then run the following codes
```
level_cluster='integrated_snn_res.0.7'
ClusterMarkers='~/Desktop/annot/step7/info7/ClusterMarkers.rds'
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
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--dgelist T
```


#### DGE contrasts
At this stage, it is possible to perform a contrast analysis on the clustered results, either at the genotype level or at the genotype-cell level

###### Genotype 
You can find a file called `${SCRNABOX_PWD}/job_info/parameters/step8_contrast_main.txt`, which contains columns for `cont_name`, `control`, `ex_control`, and `all`. You can specify genotype contrasts in this file and then use the `--genotype T` option to run the genotype contrast

```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--genotype T
```

###### Genotype-cell
To perform an interaction analysis between cell type and genotype, specify your contrast in the file `${SCRNABOX_PWD}/job_info/parameters/step8_contrast_inte.txt`. To run Step 8 using the interaction contrast, execute the following command which use the `-celltype T` option to run the Genotype-cell contrast. 
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--celltype T
```

You can integrate the contrast directly into the pipeline by calling it, 
```
CONTINT=~/des/step7_contrast_inte.txt
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--celltype T \
--cont ${CONTINT}
```

```
CONTMAIN=~/des/step7_contrast_main.txt
bash $SCRNABOX_HOME/launch_scrnabox.sh \
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

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps  integrate \
--seulist ${LISTOFSEU} \
```

LISTOFSEU includes the path of seurate files, put them in different lines. By default, the pipeline standardize the seurat objects before integrating, you can change the default in `${SCRNABOX_PWD}/job_output/parameters/stepint_par.txt`.

