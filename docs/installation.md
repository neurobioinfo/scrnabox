# Installation
To use the scRNAbox pipeline, the folowing must be installed on your High-Performance Computing (HPC) system:

- [scrnabox.slurm](#scrnaboxslurm-installation)
- [CellRanger](#cellranger-installation)
- [R and R packages](#r-library-preparation-and-r-package-installation)

 - - - -

### scrnabox.slurm installation

`scrnabox.slurm` is written in bash and can be used with any Slurm system. To download the latest version of `scrnabox.slurm` (v0.1.35) run the following command: 
```
wget https://github.com/neurobioinfo/scrnabox/releases/download/v0.1.35/scrnabox.slurm.zip
unzip scrnabox.slurm.zip
```

For a description of the options for running `scrnabox.slurm` run the following command:
```
bash /pathway/to/scrnabox.slurm/launch_scrnabox.sh -h 
```

If the `scrnabox.slurm` has been installed properly, the above command should return the folllowing:
```
        mandatory arguments:
                -d  (--dir)  = Working directory (where all the outputs will be printed) (give full path)
                --steps  =  Specify what steps, e.g., 2 to run just step 2, 2-4, run steps 2 through 4)

        optional arguments:
                -h  (--help)  = See helps regarding the pipeline options. 
                --method  = Choose what scRNA method you want to use; use HTO  and SCRNA for for hashtag nad Standard scRNA, respectively. 
                --nFeature_RNA_L  = Lower threshold of number of unique RNA transcripts for each cell, it filters nFeature_RNA > nFeature_RNA_L.  
                --nFeature_RNA_U  = Upper threshold of number of unique RNA transcripts for each cell, it filters --nFeature_RNA_U.  
                --nCount_RNA_L  = Lower threshold for nCount_RNA, it filters nCount_RNA > nCount_RNA_L   
                --nCount_RNA_U  = Upper threshold for  nCount_RNA, it filters nCount_RNA < nCount_RNA_U  
                --mitochondria_percent_L  = Lower threshold for the amount of mitochondrial transcript, it is in percent, mitochondria_percent > mitochondria_percent_L. 
                --mitochondria_percent_U  = Upper threshold for the amount of mitochondrial transcript, it is in percent, mitochondria_percent < mitochondria_percent_U. 
                --log10GenesPerUMI_U  = Upper threshold for the log number of genes per UMI for each cell, it is in percent,log10GenesPerUMI=log10(nFeature_RNA)/log10(nCount_RNA). mitochondria_percent < log10GenesPerUMI_U. 
                --log10GenesPerUMI_L  = Lower threshold for the log number of genes per UMI for each cell, log10GenesPerUMI=log10(nFeature_RNA)/log10(nCount_RNA). mitochondria_percent > log10GenesPerUMI_L.  
                --msd  = you can get the hashtag labels by running the following code 
                --marker  = Find marker. 
                --sinfo  = Do you need sample info? 
                --fta  = FindTransferAnchors 
                --enrich  = Annotation 
                --dgelist  = creates a DGEListobject from a table of counts obtained from seurate objects. 
                --genotype  = Run the genotype contrast. 
                --celltype  = Run the Genotype-cell contrast. 
                --cont  = You can directly call the contrast to the pipeline.  
                --seulist                = You can directly call the list of seurat objects to the pipeline. 
```
 - - - -

### CellRanger installation

For information regarding the installation of `CellRanger`, please visit the 10X Genomics [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation). If CellRanger is already installed on your HPC system, you may skip the CellRanger installation procedures.

 - - - -

### R library preparation and R package installation
Users must first install `R` onto their HPC system: 

```
# install R
module load r/4.2.1
```
Then, users should create a designated directory on their HPC system where the required R packages will be installed:

```
# make common R library
mkdir R_library
cd R_library

# open R
R 

# set common R library path
R_LIB_PATH="/pathway/to/R_library"
.libPaths(R_LIB_PATH)

# load packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(foreach)
library(doParallel)
library(Matrix)
library(DoubletFinder)
library(cowplot)
library(clustree)
library(xlsx)
library(enrichR)
library(stringi)
library(limma)
library(tidyverse)
library(edgeR)
library(vctrs)
library(RColorBrewer)
library(fossil)
library(openxlsx)
library(stringr)
library(ggpubr)
devtools::install_github(“neurobioinfo/scrnabox/scrnaboxR”)
```
 - - - -
Upon completing the installation procedures, users can proceed with the scRNAbox pipeline using either the [Standard scRNAseq Analysis Track](SCRNA.md) or [Cell Hashtag scRNAseq Analysis Track](HTO.md). 



