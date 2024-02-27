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
## Intepreting Log Files from the scRNAbox pipeline
## Contents

- [Introduction](#introduction)
- [Accessing Log Files](#accessing-the-log-files)
- [Contents of the Log files](#contents-of-the-log-files)
    - [Execution parameters](#execution-parameters)
    - [R function outputs](#r-function-outputs)
    - [Intermediate checkpoints](#intermediate-checkpoints)

 - - - -

## Introduction 
This guide describes the information contained in the Log files produced by each step of the scRNAbox pipeline. The Log file can be used to further investigate the processes of each step and for debugging errors.
 
- - - -
## Accessing the Log Files
Log files for each step of the scRNAbox can be found in the job_info folder of the working directory.
```
working_directory
└──job_info
    └── logs
        ├── step_1.log
        ├── step_2.log
        ├── step_3.log
        ├── step_4.log
        ├── step_5.log
        ├── step_6.log
        ├── step_7_markergsea.log
        ├── step_7_knownmarkers.log
        ├── step_7_referenceannotation.log
        ├── step_8_addmeta.log
        └── step_8_rundge.log
```

The contents of the individual log files can be visualized using the cat command

```
cd ~/working_directory/job_info/logs
cat step_1.log
```
 - - - -
## Contents of the Log files
#### Execution parameters
At the beginning of every Log file, scRNABox documents the execution parameters defined by the user for the particular run. This is intended to increase the replicability of an analysis and to allow users to adequately document their methods. Shown below is an example of the documented parameters in the Log file for step 4 of the pipeline. 
 <p align="center">
 <img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/b7c74042-e76d-4856-bd44-2ee1f612f281">
 </p>

#### R function outputs
After listing the execution parameters, the Log files will display the outputs of the R functions comprising the particulsr step of the scRNAbox pipeline; these are the same outputs that would be visible in the R console if the code was run interactively. The R function outputs will provide further information to the user beyond the outputs returned by the scRNAbox pipeline. Shown below is an example of the R function outputs in the Log file for step 4 of the pipeline.
 <p align="left">
 <img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/57a32f72-e12e-4274-8668-e234cdb30cf5" width="500" height="500">
 </p>


#### Intermediate checkpoints
Finally, to help identify the source of errors encountered while using the scRNAbox pipeline, the Log files will document the initiation and completetion of intermediate steps. For instance, Step 4 of the Standard analysis track (Doublet detection), is divided into two intermediate steps: 1) Loading R libraries and parameters and 2) Doublet detection. Before loading the required R libraries for step 4 (intermediate step 1), the log file will display **"Loading libraries and parameters started"**. Upon successfully completing intermediate step 1, the log file will display **"Loading libraries and parameters has been achieved"**. Similarly for intermediate step 2, the log file will display **"Remove doublets started"** upon initiation and **"Remove doublets has been achieved"** upon completition. If, however, an intermediate step enconters an error, the successfull completion message will not be displayed.

The intermediate checkpoints were implemented to help users narrow down the source of their errors. If, for example, users encounter an error during the doublet detection intermediate step of step 4, we suggest revisiting the parameters corresponding to doublet detection and ensurig that there are no errors in the input. If the problem persists, we invite users to open an issue in the [GitHub repository](https://github.com/neurobioinfo/scrnabox) or contact the developers directly ([Help and Feedback](contributing.md)).