# ScRNAbox: Empowering Single-Cell RNA Sequencing on High Performance Computing Systems  


[![](https://img.shields.io/badge/Documentation-scrnabox-blue)](https://neurobioinfo.github.io/scrnabox/site/)[![](https://img.shields.io/badge/Container-scrnabox-blue)](https://cloud.sylabs.io/library/saeidamiri1/mni/scrnabox.sif) [![biorxiv](https://img.shields.io/badge/biorxiv-scrnabox-blue)](https://www.biorxiv.org/content/biorxiv/early/2023/11/15/2023.11.13.566851.full.pdf) 

-------------
## Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Pipeline steps](#pipeline-steps)
- [Running scRNAbox](#running-scrnabox)
- [Tutorial](#tutorial)

---

## Introduction
ScRNAbox is a single-cell RNA sequencing (scRNAseq) pipeline specifically designed for analyzing data under a High-Performance Computing (HPC) systems using the [Slurm Workload Manager](https://slurm.schedmd.com/). The scRNAbox pipeline incorporates eight analytical steps into a comprehensive scRNAseq analysis that provides the foundation for further investigations. The eight analytical steps are outlined below. 

 <p align="center">
 <img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/bd671431-4a28-4a3a-967a-28f04580a5f8" height="700">
 </p>

The scRNAbox pipeline provides two distinct, yet highly comparable analysis tracks:

1. **Standard analysis track**
2. **HTO analysis track**

The **standard analysis track** is designed for experiments where each sample is captured and sequenced separately, while the **HTO analysis track** is designed for multiplexed experiments where samples are tagged with sample-specific oligonucleotide tagged Hashtag antibodies (HTO), pooled, and sequenced together. The **HTO analysis track** is distinguished by an additional sample demultiplexing step that assigns cells to their sample-of-origin via the sample-specific HTOs. 

<p align="center">
<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/029c737c-c696-4cd3-8a00-3c434cf4e731" height="300"> 
</p>
 
For a comprehenseive description of each step, please see the **Pipeline** section of the [scRNAbox documentation](https://neurobioinfo.github.io/scrnabox/site/) or see our [pre-print manuscript](https://www.biorxiv.org/content/10.1101/2023.11.13.566851v1). <br/>

For a tutorial that leverages the datasets used as the application cases in our pre-print manuscript, please see [scRNAbox analysis of the midbrain dataset](https://neurobioinfo.github.io/scrnabox/site/Dataset1/). 

---

## Installation
To use the scRNAbox pipeline, the folowing must be installed on your High-Performance Computing (HPC) system:

- [scrnabox.slurm](#scrnaboxslurm-installation)
- [CellRanger](#cellranger-installation)
- [R and R packages](#r-library-preparation-and-r-package-installation)

 - - - -

### scrnabox.slurm installation

`scrnabox.slurm` is written in bash and can be used with any Slurm system. To download the latest version of `scrnabox.slurm` (v0.1.52) run the following command: 
```
wget https://github.com/neurobioinfo/scrnabox/releases/download/v0.1.52/scrnabox.slurm.zip
unzip scrnabox.slurm.zip
```

For a description of the options for running `scrnabox.slurm` run the following command:
```
bash /pathway/to/scrnabox.slurm/launch_scrnabox.sh -h 
```

If the `scrnabox.slurm` has been installed properly, the above command should return the folllowing:
```
scrnabox pipeline version 0.1.52
------------------- 
        mandatory arguments:
                -d  (--dir)  = Working directory (where all the outputs will be printed) (give full path)
                --steps  =  Specify what steps, e.g., 2 to run just step 2, 2-4, run steps 2 through 4)

        optional arguments:
                -h  (--help)  = See helps regarding the pipeline arguments. 
                --method  = Select your preferred method: HTO and SCRNA for hashtag, and Standard scRNA, respectively. 
                --msd  = You can get the hashtag labels by running the following code (HTO Step 4). 
                --markergsea  = Identify marker genes for each cluster and run marker gene set enrichment analysis (GSEA) using EnrichR libraries (Step 7). 
                --knownmarkers  = Profile the individual or aggregated expression of known marker genes. 
                --referenceannotation  = Generate annotation predictions based on the annotations of a reference Seurat object (Step 7). 
                --annotate  = Add clustering annotations to Seurat object metadata (Step 7). 
                --addmeta  = Add metadata columns to the Seurat object (Step 8). 
                --rundge  = Perform differential gene expression contrasts (Step 8). 
                --seulist  = You can directly call the list of Seurat objects to the pipeline.  
 
                --rcheck  = You can identify which libraries are not installed.  
 
 ------------------- 
 For a comprehensive help, visit  https://neurobioinfo.github.io/scrnabox/site/ for documentation. 
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

Alternatively, users can install the packages manually. The R packages required for sRNAbox are shown [here](https://github.com/neurobioinfo/scrnabox/blob/main/scrnabox.slurm/soft/R/R.library.ini).
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

**Step 7: Cluster annotation** <br />
In Step 7, cluster annotation is performed to define the cell types comprising the clusters identified in Step 6. ScRNAbox provides **three tools** to identify cell types comprising the clusters:<br />

 _**Tool 1: Cluster marker gene identification and gene set enrichment analysis (GSEA)**_ <br />
 Seurat's FindAllMarkers function is used to identify differentially expressed marker genes (DEG) by the Wilcoxon rank-sum test ([Macosko et al. 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4481139/)). DEGs in the positive direction (Log2 fold-change > 0.00) are then tested for enrichment across user-defined gene set libraries that define cell types using the EnrichR tool ([Chen et al. 2013](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-128)). <br />

 _**Tool 2: Expression profiling of cell type markers and module scores**_ <br />
 Users can visualize the expression of individual genes and the aggregated expression of multiple genes. For each gene in a user-defined list, plots are produced to visualize its expression at the cluster or cell level. The aggregated expression of genes in a user-defined list are calculated using the Seurat AddModuleScore function ([Tirosh et al. 2016](https://pubmed.ncbi.nlm.nih.gov/27124452/)). <br /> 
 
 _**Tool 3: Cell type predictions based on reference data**_ <br />
 Seurat's FindTransferAnchors and TransferData functions are used to leverage cell-type annotations from a reference Seurat object and generate annotation predictions for the query dataset ([Butler et al. 2019](https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867419305598%3Fshowall%3Dtrue)). <br />


**Step 8: Differential gene expression (DGE) analysis** <br />
In Step 8, DGE analysis is computed to identify differentially expressed genes (DEG) between two conditions. ScRNAbox can compute DGE between conditions using all cell types or cell type groups. Furthermore, scRNAbox provides **two frameworks** for computing DGE: 

_**Framework 1: Cell-based DGE**_ <br />
Cells are used as replicates and DGE is computed using the Seurat FindMarkers ([Macosko et al. 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4481139/)). While FindMarkers supports several statistical frameworks to compute DGE, we set the default method in our implementation to MAST, which is tailored for scRNAseq data ([Finak et al. 2015](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0844-5)). <br />

_**Framework 2: Sample-based DGE**_ <br />
Samples are used as replicates by applying a pseudo-bulk analysis. The Seurat AggregateExpression function is used to compute the sum of RNA counts for each gene across all cells from a particular sample ([Cao et al. 2022](https://academic.oup.com/nar/article/50/21/e121/6709246?login=false)). The DESq2 statistical framework is then used to compute DGE between conditions using the aggregated counts. ([Love et al. 2014](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)). <br />

For a comprehensive decription of each step please visit scRNAbox's [documentation](https://neurobioinfo.github.io/scrnabox/site/) or see our [pre-print manuscript]().

 - - - -
 
## Running scRNAbox
To run the scRNAbox pipeline, begin by creating a dedicated folder for the analysis from the command line. Then, export the path to the working directory and the path to scrnabox.slurm:
```
mkdir working_directory
cd /pathway/to/working_directory

export SCRNABOX_HOME=/pathway/to/scrnabox.slurm
export SCRNABOX_PWD=/pathway/to/working_directory
```

Users can then run each step of the scRNAbox pipeline by adjusting the "--steps" flag in the following command:
```
cd /pathway/to/working_directory 

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 1 
```

For a comprehensive decription of how to run each step please visit scRNAbox's [documentation](https://neurobioinfo.github.io/scrnabox/site/).

 - - - -
## Tutorial ##

For a tutorial that leverages the datasets used as the application cases in our pre-print manuscript, please see [Standard analysis track: Midbrain dataset](https://neurobioinfo.github.io/scrnabox/site/Dataset1/) or [HTO analysis track: PBMC dataset](https://neurobioinfo.github.io/scrnabox/site/Dataset2/). 

---
#### Contributing
Any contributions or suggestions for improving the scRNAbox pipeline are welcomed and appreciated. You may directly contact [Saeid Amiri](https://github.com/saeidamiri1), [Michael Fiorini](https://github.com/fiorini9) or [Rhalena Thomas](https://github.com/RhalenaThomas). 

If you encounter any issues, please open an issue in the [GitHub repository](https://github.com/neurobioinfo/scrnabox).

Alternatively, you are welcomed to email the developers directly. <br />
For questions related to your HPC, please contact Saeid Amiri: [saeid.amiri@mcgill.ca]() <br />

For questions related to running the scRNAbox pipeline, please contact Michael Fiorini: [michael.fiorini@mail.mcgill.ca]() <br />

For questions related to scRNAseq analytical concepts and experimental design, please contact Rhalena Thomas: [rhalena.thomas@mcgill.ca]() <br /> 

#### Changelog
Every release is documented on the [GitHub Releases page](https://github.com/neurobioinfo/scrnabox/releases).

#### License
This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/neurobioinfo/scrnabox/blob/main/LICENSE) file for details.

#### Acknowledgement
The scRNAbox pipeline was produced for projects funded by Canadian Institute of Health Research and Michael J. Fox Foundation. The scRNAbox pipeline was produced as part Dark Genome Project. It is written by [Saeid Amiri](https://github.com/saeidamiri1), [Michael Fiorini](https://github.com/fiorini9), and [Rhalena Thomas](https://github.com/RhalenaThomas) with associate of Sali Farhan and Edward Fon at the Montreal Neurological Institute-Hospital. Copyright belongs [MNI BIOINFO CORE](https://github.com/neurobioinfo). 

**[⬆ back to top](#contents)**

