# ScRNAbox: Empowering Single-Cell RNA Sequencing on High Performance Computing Systems  


[![](https://img.shields.io/badge/Documentation-scrnabox-blue)](https://neurobioinfo.github.io/scrnabox/site/) 

-------------
## Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Pipeline steps](#pipeline-steps)
- [Tutorial](#tutorial)

---

## Introduction
ScRNAbox is a single-cell RNA sequencing (scRNAseq) pipeline specifically designed for analyzing data under a High-Performance Computing (HPC) systems using the [Slurm Workload Manager](https://slurm.schedmd.com/). The scRNAbox pipeline incorporates eight analytical steps into a comprehensive scRNAseq analysis that provides the foundation for further investigations. The eight analytical steps are outlined below. 

 <p align="center">
 <img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/3221e078-68d4-4fde-8a75-9d36691c1cf5" height="700">
 </p>

The scRNAbox pipeline provides two distinct, yet highly comparable analysis tracks:

1. **Standard analysis track**
2. **HTO analysis track**

The **standard analysis track** is designed for experiments where each sample is captured and sequenced separately, while the **HTO analysis track** is designed for multiplexed experiments where samples are tagged with sample-specific oligonucleotide tagged Hashtag antibodies (HTO), pooled, and sequenced together. The **HTO analysis track** is distinguished by an additional sample demultiplexing step that assigns cells to their sample-of-origin via the sample-specific HTOs. 

<p align="center">
<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/3a6df83e-e104-45d2-9b04-fe246642c6a8" height="300"> 
</p>
 
For a comprehenseive description of each step, please see the **Pipeline** section of the [scRNAbox documentation](https://neurobioinfo.github.io/scrnabox/site/) or see our [pre-print manuscript](). <br/>

For a tutorial that leverages the datasets used as the application cases in our pre-print manuscript, please see [scRNAbox analysis of the midbrain dataset](https://neurobioinfo.github.io/scrnabox/site/Dataset1/). 

---

## Installation
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
                -h  (--help)  = See helps regarding the pipeline arguments. 
                --method  = Select your preferred method: HTO and SCRNA for hashtag, and Standard scRNA, respectively. 
                --msd  = You can get the hashtag labels by running the following code 
                --markergsea  = Identify marker genes for each cluster and run marker gene set enrichment analysis (GSEA) using EnrichR libraries. 
                --knownmarkers  = Run module score and visualize the expression of known cell type marker genes. 
                --referenceannotation  = Run module score and visualize the expression of known cell type marker genes. 
                --annotate  = Run module score and visualize the expression of known cell type marker genes. 
                --addmeta  = Add metadata columns to the Seurat object 
                --rundge  = Perform differential gene expression contrasts 
                --seulist  = You can directly call the list of seurat objects to the pipeline.  
 
 ------------------- 
 For a comprehensive help, visit https://github.com/neurobioinfo/scrnabox for documentation.
```
 - - - -

### CellRanger installation

For information regarding the installation of `CellRanger`, please visit the 10X Genomics [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation). If CellRanger is already installed on your HPC system, you may skip the CellRanger installation procedures.

 - - - -

### R library preparation and R package installation
Users must first install `R` (v4.2 or later) onto their HPC system: 

```
# install R
module load r/4.2.1
```
Then, users must run the following installation code, which will create a directory where the R packages will be loaded and will install the required R packages:

```
# Folder for R packages 
R_PATH=~/path/to/R/library
mkdir -p $R_PATH

# Install package
Rscript ./scrnabox.slurm/soft/R/install_packages_scrnabox.R $R_PATH
```
Alternatively, users can install the packages manually. The packages required for each step of the scRNAbox pipeline are described at `./scrnabox.slurm/soft/R/R.library_hto.ini`
 - - - -

## Pipeline steps
The scRNAbox pipeline begins with 10X Genomics expression data from raw sequencing files and facilitates standard steps in scRNAseq processing through to differential gene expression between two different conditions. The pipeline is divided into 8 steps, which correspond to analytical tasks in the scRNAseq analysis workflow. Summaries of each analytical step are provided below.

**Step 1: FASTQ to gene expression matrix** <br />
In Step 1, gene expression matrices are generated from FASTQ files using the CellRanger counts pipeline. <br />

**Step 2: Create Seurat object and remove ambient RNA** <br />
In Step 2, the CellRanger outputs generated in Step 1 (expression matrix, features, and barcodes) are used to create a Seurat object for each sample. The ambient RNA quantity is estimated and there is an option to correct gene expression profiles for RNA contamination using SoupX ([Young et al. 2020](https://academic.oup.com/gigascience/article/9/12/giaa151/6049831). Then, CellRanger (if not removing ambient RNA) or SoupX (if removing ambient RNA) feature-barcode expression matrices are transformed into Seurat objects. Quality control measures are then computed to inform filtering in Step 3.<br />

**Step 3: Quality control and filtering** <br />
In Step 3, low quality cells are filtered based on the user-defined thresholds for: (i)the number of RNA unique transcripts (genes), (ii) the total number of RNA transcripts, (iii) the percentage of mitochondrial-encoded transcripts, and (iv) the percentage of ribosome gene transcripts. In addition, users can remove or regress a custom gene list from the dataset. <br />

**Step 4: Step 4: Doublet removal (standard track)** <br />
In Step 4 of the standard analysis track, doublets (barcodes produced by sequencing two or more cells) are identified and optionally removed from downstream analysis using the DoubletFinder tool ([McGinnis et al. 2019](https://www.cell.com/cell-systems/pdfExtended/S2405-4712(19)30073-0)).<br />

**Step 4: Demultiplexing and doublet detection (HTO track)** <br />
In Step 4 of the HTO track, Seurat’s implementation (MULTIseqDemux) of the tag assignment algorithm outlined in Multi-seq is used to demultiplex pooled samples and identify doublets according to the expression matrices of the sample-specific barcodes ([McGinnis et al. 2019](https://pubmed.ncbi.nlm.nih.gov/31209384/)).<br />

**Step 5: Integration**<br />
In Step 5, individual Seurat objects from each sample are combined to enable the joint analysis across samples. Users can either merge or integrate their Seurat objects ([Stuart et al. 2019](https://pubmed.ncbi.nlm.nih.gov/31178118/)). Scaling, and linear dimensional reduction is then performed to inform the optimal parameters for clustering in Step 6. <br />

**Step 6: Clustering**<br />
In Step 6, clustering is performed to define groups of cells with similar expression profiles using the Seurat implementation of the Louvain network detection with PCA dimensionality reduction as input ([Macosko et al. 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4481139/)).<br />

**Step 7: Cluster annotation** 
In Step 7, cluster annotation is performed to define the cell types comprising the clusters identified in Step 6. ScRNAbox provides three tools to identify cell types comprising the clusters:

_Tool 1: Cluster marker gene identification and gene set enrichment analysis (GSEA)_
Seurat's FindAllMarkers function is used to identify differentially expressed marker genes (DEG) by the Wilcoxon rank-sum test (Macosko et al. 2015). DEGs in the positive direction (Log2 fold-change > 0.00) are then tested for enrichment across user-defined gene set libraries that define cell types using the EnrichR tool (Chen et al. 2013). 

_Tool 2: Expression profiling of cell type markers and module scores_ 
Users can visualize the expression of individual genes and the aggregated expression of multiple genes. For each gene in a user-defined list, plots are produced to visualize its expression at the cluster or cell level. The aggregated expression of genes in a user-defined list are calculated using the Seurat AddModuleScore function (Tirosh et al. 2016. 

_Tool 3: Cell type predictions based on reference data_
Seurat's FindTransferAnchors and TransferData functions are used to leverage cell-type annotations from a reference Seurat object and generate annotation predictions for the query dataset (Macosko et al. 2015). 




Cell populations, or clusters, with similar expression profiles are annotated to define cell types by three distinct methods:<br />
    1. _Cluster marker gene set enrichment analysis (GSEA)_: Seurat's _FindAllMarkers_ function is used to identify differentially expressed marker genes (DEG) by the Wilcoxon rank-sum test (Macosko et al. 2015). DEGs in the positive direction     (Log2 fold-change > 0.00) are then tested for enrichment across user-defined gene set libraries that define cell types using the EnrichR tool (Chen et al. 2013).<br />
    2. _Reference-based annotation_: Seurat's _FindTransferAnchors_ and _TransferData_ functions are used to leverage cell-type annotations from a reference Seurat object and generate annotation predictions for the query dataset (Macosko et        al. 2015). User's must define the location of their referene Seurat object in the parameters file of Step 7.<br />
    3. _Module score_: Seurat’s implementation (_AddModuleScore_) of Tirosh et al.’s algorithm is used to comparatively quantify the expression of gene sets across clusters at the single-cell level (Tirosh et al. 2016). Users must          define their desired gene sets in the parameters file of Step 7.<br /> 
- **Step 8: Differential gene expression (DEG) contrasts** - There are multiple ways to perform differential gene expression analysis, but in this pipeline, we use the FindAllMarkers function to rank the highly differentially expressed genes in each cluster, which allows us to identify genes that are significantly differentially expressed between each cluster and the rest of the cells. From there, we can define contrasts to run statistical tests and investigate the phenotype and genotypes of each cluster.<br />

Steps 1-8 are performed using [scrnabox.slurm](https://github.com/neurobioinfo/scrnabox/tree/main/README_SCRNA.md) in the HPC system using the [Slurm Workload Manager](https://slurm.schedmd.com/).

- **Step 9: Differentially expressed genes enrichment analysis**: in this step, we obtain a list of significant genes using enrichment methods. This Step is performed using scrnaboxR. For more details see [Practice](https://github.com/neurobioinfo/scrnabox/blob/main/tutorial/temp/practice.md).


#### [Cell Hashtag scRNAseq](https://github.com/neurobioinfo/scrnabox/tree/main/README_HTO.md)

- **Step 1: FASTQ pre-processing** - Feature-barcode expression matrices are generated from FASTQ files using the CellRanger _counts_ pipeline.<br />
- **Step 2: Create Seurat object** - CellRanger- or SoupX-generated feature-barcode expression matrices are transformed into Seurat objects. Genes expressed in less than a minimum number of cells and cells expressing less than a minimum number of genes can be filtered.<br />
- **Step 3: Quality control and filtering** - Low quality cells are filtered based on the user-defined thresholds for the number of genes detected per cell,	number of unique transcripts detected per cell, percentage of mitochondrial-encoded transcripts, and percentage of ribosomal-encoded transcripts. In addition, mitochondrial- and ribosomal-encoded genes can be filtered out.<br />
- **Step 4: Demultiplexing and Doublet removal** -  Seurat’s implementation (_MULTIseqDemux_) of the tag assignment algorithm outlined in Multi-seq is used to demultiplex pooled samples and identify doublets according to the expression matrices of the sample-specific barcodes (McGinnis et al 2019).<br />
- **Step 5: Integration and linear dimensional reduction** - Individual Seurat objects are integrated to enable the joint analysis across sequencing runs using Seurat's integration algorithm (Stuart et al. 2019); if experiments are limited to a single sequencing run, the integration Step can be bypassed. Linear dimensional reduction is then performed on the resulting Seurat object to inform the optimal parameters for clustering in Step 6.<br />
- **Step 6: Clustering** - Clustering is performed to define groups of cells with similar expression profiles using the graph-based clustering approach implemented in the Seurat framework (Macosko et al. 2015).<br />
- **Step 7: Cluster annotation** - Cell populations, or clusters, with similar expression profiles are annotated to define cell types by three distinct methods:<br />
    1. _Cluster marker gene set enrichment analysis (GSEA)_: Seurat's _FindAllMarkers_ function is used to identify differentially expressed marker genes (DEG) by the Wilcoxon rank-sum test (Macosko et al. 2015). DEGs in the positive direction     (Log2 fold-change > 0.00) are then tested for enrichment across user-defined gene set libraries that define cell types using the EnrichR tool (Chen et al. 2013).<br />
    2. _Reference-based annotation_: Seurat's _FindTransferAnchors_ and _TransferData_ functions are used to leverage cell-type annotations from a reference Seurat object and generate annotation predictions for the query dataset (Macosko et        al. 2015). User's must define the location of their referene Seurat object in the parameters file of Step 7.<br />
    3. _Module score_: Seurat’s implementation (_AddModuleScore_) of Tirosh et al.’s algorithm is used to comparatively quantify the expression of gene sets across clusters at the single-cell level (Tirosh et al. 2016). Users must          define their desired gene sets in the parameters file of Step 7.<br /> 
- **Step 8: Differential gene expression (DEG) contrasts** - There are multiple ways to perform differential gene expression analysis, but in this pipeline, we use the FindAllMarkers function to rank the highly differentially expressed genes in each cluster, which allows us to identify genes that are significantly differentially expressed between each cluster and the rest of the cells. From there, we can define contrasts to run statistical tests and investigate the phenotype and genotypes of each cluster.<br />

Steps 1-8 are performed using [scrnabox.slurm](https://github.com/neurobioinfo/scrnabox/tree/main/README_SCRNA.md) in the HPC system using the [Slurm Workload Manager](https://slurm.schedmd.com/).

- **Step 9: Differentially expressed genes enrichment analysis**: in this step, we obtain a list of significant genes using enrichment methods. This Step is performed using scrnaboxR. For more details see [Practice](https://github.com/neurobioinfo/scrnabox/blob/main/tutorial/temp/practice.md).

For a comprehensive decription of each Analytical Step please visit scRNAbox's [documentation](https://neurobioinfo.github.io/scrnabox/site/).




## Tutorial
Comprehensive instructions for running both Analytical Tracks of the scRNAbox pipeline are provided [here](https://neurobioinfo.github.io/scrnabox/site/).

---
#### Contributing
This is an early version of scRNAbox and any contributions or suggestions are appreciated. To do so, you can directly contact the developers:  [Saeid Amiri](https://github.com/saeidamiri1), [Michael Fiorini](https://github.com/fiorini9), or [Rhalena Thomas](https://github.com/RhalenaThomas). 

#### Changelog
Every release is documented on the [GitHub Releases page](https://github.com/neurobioinfo/scrnabox/releases).

#### License
This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/neurobioinfo/scrnabox/blob/main/LICENSE) file for details.

#### Acknowledgement
The scRNAbox pipeline is a component of the Dark Genome project and has been developed by [Saeid Amiri](https://github.com/saeidamiri1), Michael Fiorini,  Rhalena Thomas, and Sali Farhan at  Neuro Bioinformatics Core. Copyright belongs tp MNI BIOINFO CORE (https://github.com/neurobioinfo). 

**[⬆ back to top](#contents)**

