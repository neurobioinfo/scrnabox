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
`scrnabox.slurm` is designed to start from Step1, running cellranger on fastq, although it can can be used on the the processed data (a seurat objects that some of stepts dsccusssed in are already done). The following discussed how to use the pipeline to analyse the process data. First setup the pipeline to do the analysis. Once you setup the pipeline, you should from what step you want to start. 

### From Step 2
If the cellranger is already run the raw data, copy the results under `${SCRNABOX_PWD}/step1`, then follow the step2 to analyze data 
```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 2
```

### From step 3
If the seurat objects are available and you want to run `QC and filter cell`; save the seurate objects under  `${SCRNABOX_PWD}/step2/objs` and run the following codes, 
```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 3 \
--nFRNA 300 \
--nCRNA 6500 \
--pmt 25
```

### From step 4
If the QC and filtering is already done, you can save under  `${SCRNABOX_PWD}/step3/objs` and run the following codes to get the hashtag labels by running the following code 
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

### From Step 5: Integration 
If you want to integrate the surate objects, just save them in  `${SCRNABOX_PWD}/step4/objs` and follow Srep 5
```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 5 
```

### From Step 6: Clustering 
In this step, you need to have just one seurat object which should be save in  `${SCRNABOX_PWD}/step5/objs`

```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 6 
```
