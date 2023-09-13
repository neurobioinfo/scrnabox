# scRNAbox: a comprehensive single-cell RNA sequencing pipeline designed for high-performance computing systems  
The `scrnabox.slurm` is a single-cell RNA sequencing (scRNAseq) pipeline specifically designed for analyzing data under a High-Performance Computing (HPC) systems using the [Slurm Workload Manager](https://slurm.schedmd.com/). The pipeline has been extensively utilized on the Digital Research Alliance of Canada's  [Beluga](https://docs.alliancecan.ca/wiki/B%C3%A9luga) cluster. ScRNAbox provides two distinct, yet high comparable Analysis Tracks: 
1. [Standard scRNAseq](#standard-scRNA-seq)
2. [Cell Hashtag scRNAseq](#Cell-Hashtag-scRNAseq)

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

- **Step 1: FASTQ pre-processing** - Feature-barcode expression matrices are generated from FASTQ files using the CellRanger _counts_ pipeline.<br />
- **Step 2.1: Ambient RNA removal** - The ambient RNA rate is estimated and the gene expression profiles are corrected for RNA contamination (optional) using SoupX (Young et al. 2020).<br />
- **Step 2.2: Create Seurat object** - CellRanger- or SoupX-generated feature-barcode expression matrices are transformed into Seurat objects. Genes expressed in less than a minimum number of cells and cells expressing less than a minimum number of genes can be filtered.<br />
- **Step 3: Quality control and filtering** - Low quality cells are filtered based on the user-defined thresholds for the number of genes detected per cell,	number of unique transcripts detected per cell, percentage of mitochondrial-encoded transcripts, and percentage of ribosomal-encoded transcripts. In addition, mitochondrial- and ribosomal-encoded genes can be filtered out.<br />
- **Step 4: Doublet removal** - Doublets are identified and removed from downstream analysis (optional) using the DoubletFinder tool (McGinnis et al. 2019).<br />
- **Step 5: Integration and linear dimensional reduction** - Individual Seurat objects are integrated to enable the joint analysis across sequencing runs using Seurat's integration algorithm (Stuart et al. 2019); if experiments are limited to a single sequencing run, the integration Step can be bypassed. Linear dimensional reduction is then performed on the resulting Seurat object to inform the optimal parameters for clustering in Step 6.<br />
- **Step 6: Clustering** - Clustering is performed to define groups of cells with similar expression profiles using the graph-based clustering approach implemented in the Seurat framework (Macosko et al. 2015).<br />
- **Step 7: Cluster annotation** - Cell populations, or clusters, with similar expression profiles are annotated to define cell types by three distinct methods:<br />
    1. _Cluster marker gene set enrichment analysis (GSEA)_: Seurat's _FindAllMarkers_ function is used to identify differentially expressed marker genes (DEG) by the Wilcoxon rank-sum test (Macosko et al. 2015). DEGs in the positive direction     (Log2 fold-change > 0.00) are then tested for enrichment across user-defined gene set libraries that define cell types using the EnrichR tool (Chen et al. 2013).<br />
    2. _Reference-based annotation_: Seurat's _FindTransferAnchors_ and _TransferData_ functions are used to leverage cell-type annotations from a reference Seurat object and generate annotation predictions for the query dataset (Macosko et        al. 2015). User's must define the location of their referene Seurat object in the parameters file of Step 7.<br />
    3. _Module score_: Seurat’s implementation (_AddModuleScore_) of Tirosh et al.’s algorithm is used to comparatively quantify the expression of gene sets across clusters at the single-cell level (Tirosh et al. 2016). Users must          define their desired gene sets in the parameters file of Step 7.<br /> 
- **Step 8: Differential gene expression (DEG) contrasts** - There are multiple ways to perform differential gene expression analysis, but in this pipeline, we use the FindAllMarkers function to rank the highly differentially expressed genes in each cluster, which allows us to identify genes that are significantly differentially expressed between each cluster and the rest of the cells. From there, we can define contrasts to run statistical tests and investigate the phenotype and genotypes of each cluster.<br />

Steps 1-8 are performed using [scrnabox.slurm](https://github.com/neurobioinfo/scrnabox/tree/main/README_SCRNA.md) in the HPC system using the [Slurm Workload Manager](https://slurm.schedmd.com/).

- **Step 9: Differentially expressed genes enrichment analysis**: in this step, we obtain a list of significant genes using enrichment methods. This Step is performed using scrnaboxR. For more details see [Practice](https://github.com/neurobioinfo/scrnabox/blob/main/tutorial/practice.md).


#### [Cell Hashtag scRNAseq](https://github.com/neurobioinfo/scrnabox/tree/main/README_HTO.md)

- **Step 1: FASTQ pre-processing** - Feature-barcode expression matrices are generated from FASTQ files using the CellRanger _counts_ pipeline.<br />
- **Step 2.1: Ambient RNA removal** - The ambient RNA rate is estimated and the gene expression profiles are corrected for RNA contamination (optional) using SoupX (Young et al. 2020).<br />
- **Step 2.2: Create Seurat object** - CellRanger- or SoupX-generated feature-barcode expression matrices are transformed into Seurat objects. Genes expressed in less than a minimum number of cells and cells expressing less than a minimum number of genes can be filtered.<br />
- **Step 3: Quality control and filtering** - Low quality cells are filtered based on the user-defined thresholds for the number of genes detected per cell,	number of unique transcripts detected per cell, percentage of mitochondrial-encoded transcripts, and percentage of ribosomal-encoded transcripts. In addition, mitochondrial- and ribosomal-encoded genes can be filtered out.<br />
- **Step 4: Demultiplexing and Doublet removal** -  Seurat’s implementation (_MULTIseqDemux_) of the tag assignment algorithm outlined in Multi-seq to demultiplex pooled samples and identify doublets according to the expression matrices of the sample-specific barcodes (McGinnis et al 2019).<br />
- **Step 5: Integration and linear dimensional reduction** - Individual Seurat objects are integrated to enable the joint analysis across sequencing runs using Seurat's integration algorithm (Stuart et al. 2019); if experiments are limited to a single sequencing run, the integration Step can be bypassed. Linear dimensional reduction is then performed on the resulting Seurat object to inform the optimal parameters for clustering in Step 6.<br />
- **Step 6: Clustering** - Clustering is performed to define groups of cells with similar expression profiles using the graph-based clustering approach implemented in the Seurat framework (Macosko et al. 2015).<br />
- **Step 7: Cluster annotation** - Cell populations, or clusters, with similar expression profiles are annotated to define cell types by three distinct methods:<br />
    1. _Cluster marker gene set enrichment analysis (GSEA)_: Seurat's _FindAllMarkers_ function is used to identify differentially expressed marker genes (DEG) by the Wilcoxon rank-sum test (Macosko et al. 2015). DEGs in the positive direction     (Log2 fold-change > 0.00) are then tested for enrichment across user-defined gene set libraries that define cell types using the EnrichR tool (Chen et al. 2013).<br />
    2. _Reference-based annotation_: Seurat's _FindTransferAnchors_ and _TransferData_ functions are used to leverage cell-type annotations from a reference Seurat object and generate annotation predictions for the query dataset (Macosko et        al. 2015). User's must define the location of their referene Seurat object in the parameters file of Step 7.<br />
    3. _Module score_: Seurat’s implementation (_AddModuleScore_) of Tirosh et al.’s algorithm is used to comparatively quantify the expression of gene sets across clusters at the single-cell level (Tirosh et al. 2016). Users must          define their desired gene sets in the parameters file of Step 7.<br /> 
- **Step 8: Differential gene expression (DEG) contrasts** - There are multiple ways to perform differential gene expression analysis, but in this pipeline, we use the FindAllMarkers function to rank the highly differentially expressed genes in each cluster, which allows us to identify genes that are significantly differentially expressed between each cluster and the rest of the cells. From there, we can define contrasts to run statistical tests and investigate the phenotype and genotypes of each cluster.<br />

Steps 1-8 are performed using [scrnabox.slurm](https://github.com/neurobioinfo/scrnabox/tree/main/README_SCRNA.md) in the HPC system using the [Slurm Workload Manager](https://slurm.schedmd.com/).

- **Step 9: Differentially expressed genes enrichment analysis**: in this step, we obtain a list of significant genes using enrichment methods. This Step is performed using scrnaboxR. For more details see [Practice](https://github.com/neurobioinfo/scrnabox/blob/main/tutorial/practice.md).

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

