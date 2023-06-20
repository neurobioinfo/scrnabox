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
## ScRNA  pipline  
## Contents

- [Introduction](#introduction)
  - [From step 2](#from-step-2)
  - [From step 3](#from-step-3)  
  - [From step 4](#from-step-4)
  - [From Step 5: integration](#from-step-5-integration)
  - [From Step 6: Clustering](#from-step-6-clustering)   

## Introduction 
This guide provides a concise introduction to analyzing data using the Scrnabox pipeline. The `scrnabox.slurm` pipeline is primarily designed to initiate Step 1, which involves running cellranger on fastq data. However, it can also be utilized with processed data, where some of the steps have already been completed. The following section explains how to use the pipeline for analyzing processed data. To begin, you need to set up the pipeline for analysis and determine the starting step based on your requirements.

### From Step 2
If the cellranger is already run the raw data, copy the results under `${SCRNABOX_PWD}/step1`, then follow the step2 to analyze data 
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 2
```

### From step 3
If the seurat objects are available and you want to run `QC and filter cell`; save the seurate objects under  `${SCRNABOX_PWD}/step2/objs` and run the following codes, 
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 3 \
--nFeature_RNA_L 1000 \
--nCount_RNA_U 65000 \
--mitochondria_percent_U 25
```

### From step 4
If the QC and filtering is already done, you can save under  `${SCRNABOX_PWD}/step3/objs` and run the following codes to get the hashtag labels by running the following code 
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 4 \
--msd T 
```

Add the old label and its new corresponding label in '${SCRNABOX_PWD}/job_output/parameters/step4_par.txt'. Run the following to run the demuplixing  
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 4 
```

### From Step 5: Integration 
If you want to integrate the surate objects, just save them in  `${SCRNABOX_PWD}/step4/objs` and follow Srep 5
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 5 
```

### From Step 6: Clustering 
In this step, you need to have just one seurat object which should be save in  `${SCRNABOX_PWD}/step5/objs`

```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 6 
```
