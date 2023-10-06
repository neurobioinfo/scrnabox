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
## scRNAbox with processed data
## Contents

- [Introduction](#introduction)
  - [From Step 2: Create Seurat object and remove ambient RNA](#from-step-2-create-seurat-object-and-remove-ambient-rna)  
  - [From Step 3: Quality control and filtering](#from-step-3quality-control-and-filtering)
  - [From Step 4: Demultiplexing and doublet removal](#from-step-4-demultiplexing-and-doublet-removal)
  - [From Step 5: Integration and linear dimensional reduction](#from-step-5-integration-and-linear-dimensional-reduction)  
  - [From Step 6: Clustering](#from-step-6-clustering)
  - [From Step 7: Cluster annotation](#from-step-7-cluster-annotation) 
  - [From Step 8: Differential gene expression contrasts](#from-step-8-differential-gene-expression-contrasts)    

 - - - -

## Introduction 
This guide demonstrates how to initiate the scRNAbox pipeline at different Analytical Steps using processed data. The procedures are the same for both the [Standard](SCRNA.md) and [Cell Hashtag](HTO.md) Analysis Tracks.
 - - - -
### From Step 2: Create Seurat object and remove ambient RNA

**I have to return to this.**

 - - - -

### From Step 3: Quality control and filtering
If you have a Seurat object(s) and want to initiate the pipeline at Step 3, create a `/step3/objs3` folder in the working directory, copy the Seurat object(s) to this folder, and run Step 3.

```
export SCRNABOX_PWD=/path/to/working/directory
cd ${SCRNABOX_PWD}
mkdir step3
cd step3
mkdir objs3
cp /path/to/Seurat/object ${SCRNABOX_PWD}/step3/objs3

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 3
```
- - - -


### From Step 4: Demultiplexing and doublet removal
If you have a Seurat object(s) and want to initiate the pipeline at Step 4, create a `/step4/objs4` folder in the working directory, copy the Seurat object(s) to this folder, and run Step 4.

```
export SCRNABOX_PWD=/path/to/working/directory
cd ${SCRNABOX_PWD}
mkdir step4
cd step4
mkdir objs4
cp /path/to/Seurat/object ${SCRNABOX_PWD}/step4/objs4

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 4
```
- - - -
### From Step 5: Integration and linear dimensional reduction
If you have a Seurat object(s) and want to initiate the pipeline at Step 5, create a `/step5/objs5` folder in the working directory, copy the Seurat object(s) to this folder, and run Step 5.

```
export SCRNABOX_PWD=/path/to/working/directory
cd ${SCRNABOX_PWD}
mkdir step5
cd step5
mkdir objs5
cp /path/to/Seurat/object ${SCRNABOX_PWD}/step5/objs5

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 5
```
- - - -

### From Step 6: Clustering 
If you have a Seurat object and want to initiate the pipeline at Step 6, create a `/step6/objs6` folder in the working directory, copy the Seurat object to this folder, and run Step 6.

```
export SCRNABOX_PWD=/path/to/working/directory
cd ${SCRNABOX_PWD}
mkdir step6
cd step6
mkdir objs6
cp /path/to/Seurat/object ${SCRNABOX_PWD}/step6/objs6

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 6
```
- - - -

### From Step 7: Cluster annotation
If you have a Seurat object and want to initiate the pipeline at Step 7, create a `/step7/objs7` folder in the working directory, copy the Seurat object to this folder, and run Step 7.

```
export SCRNABOX_PWD=/path/to/working/directory
cd ${SCRNABOX_PWD}
mkdir step7
cd step7
mkdir objs7
cp /path/to/Seurat/object ${SCRNABOX_PWD}/step7/objs7

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7
```
- - - -

### From Step 8:  Differential gene expression contrasts
If you have a Seurat object and want to initiate the pipeline at Step 8, create a `/step8/objs8` folder in the working directory, copy the Seurat object to this folder, and run Step 8.

```
export SCRNABOX_PWD=/path/to/working/directory
cd ${SCRNABOX_PWD}
mkdir step8
cd step8
mkdir objs8
cp /path/to/Seurat/object ${SCRNABOX_PWD}/step8/objs8

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8
```
- - - -