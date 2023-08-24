---
layout: post
title:  A Guide to Analyzing HTO with the Scrnabox Pipeline
description: A short introduction to  Hashtag oligonucleotide analyzing using scrnabox pipeline
date: 2023-06-16
author: Saeid Amiri
published: true
tags: scRNA FAQ
categories: 
comments: false
---
## Frequently asked questions

- [Possibility to consider multiple control](#possibility-to-consider-multiple-control)
- [How skip the integration step](#how-skip-the-integration-step)
- [Can we run without cellraanger](#can-we-run-without-cellraanger)
## Possibility to consider multiple control
User can redefine the genotype; in the step8, add `/job_info/parameters/step8_par.txt`, add `new_genotype='YES` and redefine the labels: 

```
new_genotype='YES'
old_antibody_label=c('B0251-TotalSeqB','B0252-TotalSeqB','B0253-TotalSeqB','B0254-TotalSeqB','B0255-TotalSeqB','B0256-TotalSeqB')
new_antibody_label=c('AIW002','SNCA-A53T','GBA-KO','Parkin-KO','PINK1-KO','SNCA-KO')
```

## How to skip the integration step
If you have a sole sample, there is no need for the integration step; once you run step 4, go to step 6 and add `par_skip_step5='YES'` to the step 6 parameter, `/job_info/parameters/step6_par.txt`.

## Can we run without cellraanger. 
if you have the matrix file, .......
