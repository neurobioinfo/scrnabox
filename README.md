# scRNAbox: a comprehensive single-cell RNA sequencing pipeline designed for high-performance computing systems  
The `scrnabox.slurm` is a single-cell RNA sequencing (scRNAseq) pipeline specifically designed for analyzing data under a High-Performance Computing (HPC) systems  using the Slurm Workload Manager ([slurm work load manager system](https://slurm.schedmd.com/)). The pipeline has been extensively utilized on the Digital Research Alliance of Canada's  [Beluga](https://docs.alliancecan.ca/wiki/B%C3%A9luga) cluster. ScRNAbox provides two distinct, yet high comparable Analysis Tracks: 
1. [Standard scRNAseq](#standard-scRNA-seq)
2. [Cell Hashtag scRNAseq](#cell-hashtags)

The Standard Analysis Track is designed for experiments where each sample is captured and sequenced separately, while the Cell Hashtag Analysis Track is designed for multiplexed experiments, whereby samples are tagged with sample-specific barcodes, pooled, and sequenced together. The Cell Hashtag Analysis Track is distinguished by an additional sample demultiplexing Step that assigns cells to their sample-of-origin via the sample-specific barcodes. 

<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/3a6df83e-e104-45d2-9b04-fe246642c6a8" height="300">

Please refer to the [documentation](https://neurobioinfo.github.io/scrnabox/site/) for comprehensive instructions regarding the application of each Analysis Track.   

## Contents
- [scRNAbox analysis workflow](#scRNAbox-analysis-workflow)
  - [Standard scRNAseq](#standard-scRNAseq)
  - [Cell Hashtag scRNAseq](#cell-hashtags)
- [Installation](#installing)
- [Tutorial](#tutorial)


---

## scRNAbox analysis workflow
The following figure illustrates the Anlytical Steps comprising each Analysis Track – Standard scRNAseq and Cell Hashtag scRNAseq – of the the scRNAbox pipeline. Prior to running each Analytical Step, users can adjust the execution parameters in the Step-specific parameters text file, which is automatically downloaded upon [Installation](#installing). Following each Analytical Step, intermediate Seurat objects are generated and, where applicable, results are reported as intuitive summary files, tables, or figures, deposited directly into the working directory. 

<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/3729254c-0ca1-4866-aa27-1bda6129e7ca" height="800">

Summaries of each Analytical Step comprising the [Standard scRNAseq](#standard-scRNA-seq) and [Cell Hashtag scRNAseq](#cell-hashtags) Analysis Tracks are provided below.

#### [Standard scRNAseq](https://github.com/neurobioinfo/scrnabox/tree/main/README_SCRNA.md)

The following steps describe how to analyze scRNA-seq data using the pipeline:<br />
- **Step 1: FASTQ pre-processing** - Feature-barcode expression matrices are generated from FASTQ files using the CellRanger _counts_ pipeline.<br />
- **Step 2.1: Ambient RNA removal** - The ambient RNA rate is estimated and the gene expression profiles are corrected for RNA contamination (optional) using SoupX (Young et al. 2020).<br />
- **Step 2.2: Create Seurat object** - CellRanger- or SoupX-generated feature-barcode expression matrices are transformed into Seurat objects. Genes expressed in less than a minimum number of cells and cells expressing less than a minimum number of genes can be filtered.<br />
- **Step 3: Quality control and filtering** - Low quality cells are filtered based on the user-defined thresholds for the number of genes detected per cell,	number of unique transcripts detected per cell, percentage of mitochondrial-encoded transcripts, and percentage of ribosomal-encoded transcripts. In addition, mitochondrial- and ribosomal-encoded genes can be filtered out.<br />
- **Step 4: Doublet removal** - Doublets are identified and removed from downstream analysis (optional) using the DoubletFinder tool (McGinnis et al. 2019).<br />
- **Step 5: Integration and linear dimensional reduction** - Individual Seurat objects are integrated to enable the joint analysis across sequencing runs using Seurat's integration algorithm (Stuart et al. 2019); if experiments are limited to a single sequencing run, the integration Step can be bypassed. Linear dimensional reduction is then performed on the resulting Seurat object to inform the optimal parameters for clustering in Step 6.<br />
- **Step 6: Clustering** - Clustering is performed to define groups of cells with similar expression profiles using the graph-based clustering approach implemented in the Seurat framework (Macosko et al. 2015).<br />
- **Step 7: Cluster annotation** - Cell populations, or clusters, with similar expression profiles are annotated to define cell types by three distinct methods:<br />
1. _Cluster marker gene set enrichment analysis (GSEA)_: Seurat's _FindAllMarkers_ function is used to identify differentially expressed marker genes (DEG) by the Wilcoxon rank-sum test (Macosko et al. 2015). DEGs in the positive direction (Log2 fold-change > 0.00) are then tested for enrichment across user-defined gene set libraries that define cell types using the EnrichR tool (Chen et al. 2013).<br />
2. _Reference-based annotation_: Seurat's _FindTransferAnchors_ and _TransferData_ functions are used to leverage cell-type annotations from a reference Seurat object and generate annotation predictions for the query dataset (Macosko et al. 2015). User's must define the location of their referene Seurat object in the parameters file of Step 7.<br />
3. _Cluster module score_: Seurat’s implementation (_AddModuleScore_) of Tirosh et al.’s algorithm is used to comparatively quantify the expression of gene sets across clusters at the single-cell level (Tirosh et al. 2016). Users must define their desired gene sets in the parameters file of Step 7.<br /> 
- **Step 8: Differential gene expression (DEG) contrasts** - There are multiple ways to perform differential gene expression analysis, but in this pipeline, we use the FindAllMarkers function to rank the highly differentially expressed genes in each cluster, which allows us to identify genes that are significantly differentially expressed between each cluster and the rest of the cells. From there, we can define contrasts to run statistical tests and investigate the phenotype and genotypes of each cluster.<br />
- step 9, Enrichment analysis: in this step, we obtain a list of significant genes using enrichment methods. The step 9 can be done using scrnaboxR, see [Practice](https://github.com/neurobioinfo/scrnabox/blob/main/tutorial/practice.md).

The Step 1 - Step 8 can be done using [scrnabox.slurm](https://github.com/neurobioinfo/scrnabox/tree/main/README_SCRNA.md) in the HPC system ([slurm work load manager system](https://slurm.schedmd.com/)).

#### [Cell Hashtags](https://github.com/neurobioinfo/scrnabox/tree/main/README_HTO.md)

The following steps explain how to analyze the Hashtag oligonucleotide (henceforth referred to as HTO)
- Step 1: cellranger - This step runs Cellranger on the scRNA-seq data to generate the feature-barcode matrices for each sample.
- Step 2: Seurat object - This step runs Seurat on the feature-barcode matrices obtained from step 1 to generate a Seurat object for each sample. The Seurat object contains a lot of information, and accessing SeuratObject@meta.data will provide a data frame with relevant information on each cell.
- Step 3: QC and filter - The Seurat object includes quality metrics that can be used to filter cells and genes against possible doublets. Metrics such as total UMI counts per cell (nCount_RNA), total number of detected features per cell (nFeature_RNA), and mitochondrial count (percent.mito) are often used.
- Step 4: Demuplixing- In this step, demultiplexing is performed to separate the reads in the sequencing run according to their sample of origin, based on their barcode information. The pipeline includes an option to remove doublets and negative cells during a later step after quality control and filtering.
- Step 5: Integration - This step integrates multiple scRNA-seq datasets using the Comprehensive Integration of Single Cell Data (CCA) method in Seurat, Tim, et al. (2019). The pipeline identifies anchors using the FindIntegrationAnchors function and passes them to the IntegrateData function to get a single Seurat object representing all the datasets.
- Step 6: Clustering- It involves clustering the data using a k-nearest neighbor graph based on the integrated PCA. This step produces UMAP and heatmaps of unlabelled clusters. However, it is up to the user to decide on the best cluster resolution outside the pipeline by examining the output and selecting the most appropriate annotation for the clusters.
- step 7: Cluster annotate - In this step, you can annotate the clusters with known cell types or use marker genes to predict the cell type of each cluster. 
- step 8: Differential gene expression (DEG)- There are multiple ways to perform differential gene expression analysis, but in this pipeline, we use the FindAllMarkers function to rank the highly differentially expressed genes in each cluster, which allows us to identify genes that are significantly differentially expressed between each cluster and the rest of the cells. From there, we can define contrasts to run statistical tests and investigate the phenotype and genotypes of each cluster.

The Step 1 - Step 8 can be done using [scrnabox.slurm](https://github.com/neurobioinfo/scrnabox/tree/main/README_HTO.md) in the HPC system ([slurm work load manager system](https://slurm.schedmd.com/)).
- step 9, Enrichment analysis: in this step, we obtain a list of significant genes using enrichment methods. The step 9 can be done using scrnaboxR, see [Practice](https://github.com/neurobioinfo/scrnabox/blob/main/tutorial/practice.md).

For a comprehensive decription of each Analytical Step please visit scRNAbox's [documentation](https://neurobioinfo.github.io/scrnabox/site/)

## Installing
The package is written in the bash, so it can be used with any slurm system. To download 
`scrnabox.slurm` run the below comments 
```
wget https://github.com/neurobioinfo/scrnabox/releases/download/v0.1.0/scrnabox.slurm.zip
unzip scrnabox.slurm.zip 
```

To obtain a brief guidance of the pipeline, execute the following code.
```
bash ./scrnabox.slurm/launch_scrnabox.sh -h 
```

`scrnabox.slurm` needs `R` and `cellranger`. For the R, you need to install 
`'Seurat','ggplot2', 'dplyr', 'foreach', 'doParallel', 'Matrix', 'DoubletFinder','cowplot','clustree'`. Then install `'scrnaboxR'`: 
```
devtools::install_github("neurobioinfo/scrnabox/scrnaboxR")
```
The `'scrnaboxR'` is an R package that provides a collection of functions for conducting enrichment analysis and other analyses associated with single-cell RNA sequencing (scRNA-seq) data. It serves as a companion to scrnabox, offering a range of tools and functionalities to enhance scRNA-seq data analysis. You need to add the R info in `scrnabox_config.ini`, you can define the path of R library in `R_LIB_PATH=`, version of R in `R_VERSION`, you can add the path of `cell ranger`in `MODULECELLRANGER` 


## Tutorial
You can find the details of how to use the pipeline in the 
the [documentation](https://neurobioinfo.github.io/scrnabox/site/).


---
#### Contributing
This is an early version, any contribute or suggestion is appreciated, you can directly contact with developers:  [Saeid Amiri](https://github.com/saeidamiri1), [Michael Fiorini](https://github.com/fiorini9), or [Rhalena Thomas](https://github.com/RhalenaThomas). 

#### Changelog
Every release is documented on the [GitHub Releases page](https://github.com/neurobioinfo/scrnabox/releases).

#### License
This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/neurobioinfo/scrnabox/blob/main/LICENSE) file for details.

#### Acknowledgement
The pipeline is done as part Dark Genome project, it is written by [Saeid Amiri](https://github.com/saeidamiri1) with associate of Rhalena Thomas, Sali Farhan, and Michael Fiorini at  Neuro Bioinformatics Core. Copyright belong MNI BIOINFO CORE (https://github.com/neurobioinfo). 

**[⬆ back to top](#contents)**

