# scrnabox: A pipeline for scRNA 
This repository includes

## Contents
- [Analysis Flowchart]
- [scrnabox.svn]
- [scrnaboxR]
- [Rmarkdown]

## Analysis workflow
The analysis can be done via 

### Hashtag
The following figure shows the steps to analyze the hashtag scRNA
![hashtag](https://raw.githubusercontent.com/neurobioinfo/scrnabox/main/hashtag.png)

- Step 1, cellranger:  this step run cellranger  on the single cell data. 
- Step 2, Seurat object: this step run seurat on feature-barcode matrices obtained from step 1 to generate seurat object for each sample.  Seurat's object includes a lot of information 
SeuratObject@meta.data will return the data frame and relevant information on each cell. 
- Step 3, QC and filter:  Seurat object includes some quality measure that can be used to filter cell and genes against possible doublets, we often use the total UMI counts per cell (nCount_RNA), the total number of detected features per cell (nFeature_RNA), and  mitochondrial count (percent.mito). 
- Step 4, demuplixing: 
- Step 5, integration: this step integrates multiple single cell RNA-seq datasets.  Seurat uses the Comprehensive Integration of Single Cell Data (CCA) to perform integration; we identify anchors using the FindIntegrationAnchors function and pass them to the IntegrateData function to get a Seurat object.
- Step 6, Clustering: here, we run clustering (a k-nearest neighbour graph) on the intergrated PCA. 
- step 7,  Differetial gene expression (DEG):  DEG can be done in different ways, here we  run the function FindAllMarkers to compute a ranking for the highly differential genes in each cluster which determines the genes differentially expressed between each cluster and the rest of the cells. Then define contrast to run statstical tests to study the phenotype and genotypes. 
- step 8, Enrichment analysis: in this step, we obtain list of significant genes using .  


The step 1 - Step 7 can be done using `scrnabox.svn` and step 8 scrnaboxR. 

### Non-hashtaq
????????
-----------
## scrnabox.svn
`scrnabox.svn` is a pipeline developed to run step 1 to step 7 under HPC system, we are using the pipeline under [Beluga](https://docs.alliancecan.ca/wiki/B%C3%A9luga)

  - [Setup](#setup)
  - [Step 1: cellranger](#step-1-cellranger)
  - [Step 2: Seurat object](#step-2-seurat-object)  
  - [Step 3: QC and filter](#step-3-qc-and-filter)
  - [Step 4: demuplixing](#step-4-demuplixing)
  - [Step 5: integration](#step-5-integration)
  - [Step 6: Clustering](#step-6-clustering)   
  - [step 7: DEG contrast](#step-7-DEG-contrast)  

In order to run the pipeline, first create a folder to do the analysis and export the pipeline
```
mkdir -P  ~/scratch/des
export SCRNABOX_HOME=/lustre03/project/6070393/COMMON/Dark_Genome/samamiri/pipeline/scrnabox.svn
SCRNABOX_PWD=~/scratch/des
```

### Setup
Once its 'SCRNABOX_PWD' is defined, you need to create a folder entitled `samples_info` and write samples's `library.csv` and `features_ref.csv`. Then run the following code to setup pipelione, 
```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 0 
```
The pipeline creates `./job_output/scrnabox.config.ini` (include the configure arguments) and `./job_output/expected.done.files.txt` (recorder the done steps). 

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
This step run QC and save the results under `${SCRNABOX_PWD}/step3`. The followings do `nFeature_RNA > 300 & nFeature_RNA < 6500 & percent.mt < 25`. 
```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 3 \
--nFRNAl 300 \
--nFRNAu 6500 \
--pmt 25
```

### Step 4: Demuplixing 
In this step, you need to choose the right label, you can get the label by running the following code 
```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 4 \
--msd T 
```

Add the old label and its new corresponding label in '${SCRNABOX_PWD}/job_output/step4_par.txt'. Run the following to run the demuplixing  

```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 4 
```

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


### step 7: DGE contrast
In this step, one can run the contrast on clustered result, which can be done on genotype and genotype-cell and are referred as main and interact. First add the label to /job_output/step7_ clus_label.txt. 
#### genotype 
There is a file ${SCRNABOX_PWD}/job_output/step7_contrast_main.txt, with columns of cont_name,control,ex_control,all, you can write the genotype contrast here, then select `--main T` to run the genotype contrast. 
```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--main T
```

#### genotype-cell
To run interact between celltype and genptype, write your contrast in `${SCRNABOX_PWD}/job_output/step7_contrast_inte.txt`. To run Step 7 on interact contrast, run the following command. Select `--inte T` to run the main contrast. 
```
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--inte T
```

You can directly call the contrast to the pipeline, 
```
CONTINT=/lustre03/project/6070393/COMMON/Dark_Genome/samamiri/pipeline/scrnabox.svn_run/des/step7_contrast_inte.txt
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--inte T \
--cont ${CONTINT}
```

```
CONTMAIN=/lustre03/project/6070393/COMMON/Dark_Genome/samamiri/pipeline/scrnabox.svn_run/des/step7_contrast_main.txt
sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--main T \
--cont ${CONTMAIN}
```



## scrnaboxR
The R package can be downloaded using the following script. 
```
devtools::install_github("neurobioinfo/scrnabox/scrnaboxR")
```


## Rmarkdown
.....


## Contributing
This is an early version, any contribute or suggestion is appreciated, you can directly contact with [Saeid Amiri] or [Rhalena Thompson] 
## Changelog
Every release is documented on the [GitHub Releases page](https://github.com/neurobioinfo/scrnabox/releases).
## License
This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/neurobioinfo/scrnabox/blob/main/LICENSE) file for details
## Acknowledgement
## Todo

**[⬆ back to top](#contents)**
