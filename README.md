# scrnabox: A pipeline for scRNA analysis under HPC  
This repository includes the tutorial for the scRNA-seq data using scrnabox pipeline under HPC ([slurm work load manager system](https://slurm.schedmd.com/)). 

## Contents
- [Workflow analysis](#analysis-workflow)
  - [Standard scRNA-seq]
  - [Cell hashtags]
- [scrnabox.slurm](#scrnaboxsvn)
- [scrnaboxR](#scrnaboxr)
- [Processed Data](#processed-data)
- [Tutorial](#tutorial)
- [References](#references)


## Analysis workflow
The following figure shows the steps to analyze the Standard and Cell Hashtags scRNA-seq data:

![Steps of Standard scRNA-seq ](https://raw.githubusercontent.com/neurobioinfo/scrnabox/main/figs/scrna.png)

![Steps of Cell Hashtags scRNA-seq](https://raw.githubusercontent.com/neurobioinfo/scrnabox/main/figs/hto.png)

At the moment, the standard and cell hashtags are implemented in the pipeline.  

#### [Standard scRNA-seq](https://github.com/neurobioinfo/scrnabox/tree/main/README_SCRNA.md)
The following steps explain how to analyze the, 
- Step 1, cellranger: this step runs cellranger on the scRNA-seq data.
- Step 2, Seurat object: this step runs seurat on feature-barcode matrices obtained from step 1 to generate seurat object for each sample. Seurat's object includes a lot of information; SeuratObject@meta.data will return the data frame and relevant information on each cell.
- Step 3, QC and filter: Seurat object includes some quality measure that can be used to filter cell and genes against possible doublets, we often use the total UMI counts per cell (nCount_RNA), the total number of detected features per cell (nFeature_RNA), and mitochondrial count (percent.mito).
- Step 4, This step can be used to remove the 'Doublet'; the defalut is to remove it. 
- Step 5, integration: this step integrates multiple single cell RNA-seq datasets. Seurat uses the Comprehensive Integration of Single Cell Data (CCA), see Stuart, Tim, et al. (2019) to perform integration; we identify anchors using the FindIntegrationAnchors function and pass them to the IntegrateData function to get a Seurat object.
- Step 6, Clustering: here, we run clustering (a k-nearest neighbour graph) on the integrated PCA. At the end of this step we have some UMAP and heatmaps of unlabelled clusters. We need to decide on the cluster annotation outside the pipeline. User needs to look at the output and decide on the best cluster resolution.
- step 7, Cluster annotate: you can annotate the cluster with cell types in this step.
- step 8, Differential gene expression (DEG): this step can be done in different ways, here we run the function FindAllMarkers to compute a ranking for the highly differential genes in each cluster which determines the genes differentially expressed between each cluster and the rest of the cells. Then define contrast to run statistical tests to study the phenotype and genotypes.

The Step 1 - Step 8 can be done using [scrnabox.slurm](https://github.com/neurobioinfo/scrnabox/tree/main/README_SCRNA.md) in the HPC system ([slurm work load manager system](https://slurm.schedmd.com/)).
- step 9, Enrichment analysis: in this step, we obtain a list of significant genes using enrichment methods. The step 9 can be done using scrnaboxR, see [Practice](https://github.com/neurobioinfo/scrnabox/blob/main/tutorial/practice.md).


#### [Cell Hashtags](https://github.com/neurobioinfo/scrnabox/tree/main/README_HTO.md)
The following steps explain how to analyze the  Hashtag oligonucleotide (henceforth referred to as HTO)
- Step 1, cellranger: this step runs cellranger on the single cell data.
- Step 2, Seurat object: this step runs seurat on feature-barcode matrices obtained from step 1 to generate seurat object for each sample. Seurat's object includes a lot of information; SeuratObject@meta.data will return the data frame and relevant information on each cell.
- Step 3, QC and filter: Seurat object includes some quality measure that can be used to filter cell and genes against possible doublets, we often use the total UMI counts per cell (nCount_RNA), the total number of detected features per cell (nFeature_RNA), and mitochondrial count (percent.mito).
- Step 4, demuplixing: the process uses the barcode information in order to know which sequences came from which samples after they had all been sequenced together. Pipeline has an option to remove the doublets and negative cells as this step.
- Step 5, integration: this step integrates multiple single cell RNA-seq datasets. Seurat uses the Comprehensive Integration of Single Cell Data (CCA), see Stuart, Tim, et al. (2019) to perform integration; we identify anchors using the FindIntegrationAnchors function and pass them to the IntegrateData function to get a Seurat object.
- Step 6, Clustering: here, we run clustering (a k-nearest neighbour graph) on the integrated PCA. At the end of this step we have some UMAP and heatmaps of unlabelled clusters. We need to decide on the cluster annotation outside the pipeline. User needs to look at the output and decide on the best cluster resolution.
- step 7, Cluster annotate: you can annotate the cluster with cell types in this step.
- step 8, Differential gene expression (DEG): this step can be done in different ways, here we run the function FindAllMarkers to compute a ranking for the highly differential genes in each cluster which determines the genes differentially expressed between each cluster and the rest of the cells. Then define contrast to run statistical tests to study the phenotype and genotypes.

The Step 1 - Step 8 can be done using [scrnabox.slurm](https://github.com/neurobioinfo/scrnabox/tree/main/README_HTO.md) in the HPC system ([slurm work load manager system](https://slurm.schedmd.com/)).
- step 9, Enrichment analysis: in this step, we obtain a list of significant genes using enrichment methods. The step 9 can be done using scrnaboxR, see [Practice](https://github.com/neurobioinfo/scrnabox/blob/main/tutorial/practice.md).


## scrnabox.slurm
`scrnabox.slurm` is a pipeline developed to run step 1 to step 8 under HPC system ([slurm work load manager system](https://slurm.schedmd.com/)), the pipeline has been using under [Beluga](https://docs.alliancecan.ca/wiki/B%C3%A9luga), the detail of how to use it is discussed in [scrnabox.slurm](https://github.com/neurobioinfo/scrnabox/tree/main/scrnabox.slurm). 

## scrnaboxR
The `scrnaboxR` is a R package of some functions to run the enrichment analysis and other analysis. This library can be installed using the following script. 
```
devtools::install_github("neurobioinfo/scrnabox/scrnaboxR")
```

## Processed Data
The pipeline can be used to work with the [processed data](https://github.com/neurobioinfo/scrnabox/blob/main/README_PROC.md), i.e., you can bring data from different project and skip the proposed steps in the analysis data. 

## Tutorial
To analyze DGE look at [Practice](https://github.com/neurobioinfo/scrnabox/blob/main/tutorial/practice.md). 

## References
- Stuart, T., Butler, A., Hoffman, P., Hafemeister, C., Papalexi, E., Mauck III, W. M., ... & Satija, R. (2019). Comprehensive integration of single-cell data. Cell, 177(7), 1888-1902.

## Contributing
This is an early version, any contribute or suggestion is appreciated, you can directly contact with [Saeid Amiri](https://github.com/saeidamiri1) or [Rhalena Thomas](https://github.com/RhalenaThomas). 

## Changelog
Every release is documented on the [GitHub Releases page](https://github.com/neurobioinfo/scrnabox/releases).
## License
This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/neurobioinfo/scrnabox/blob/main/LICENSE) file for details
## Acknowledgement
The pipeline is done as part Dark Genome project, it is written by [Saeid Amiri](https://github.com/saeidamiri1) with associate of Rhalena Thomas and  Roxanne Larivière. 

## Todo
**[⬆ back to top](#contents)**
