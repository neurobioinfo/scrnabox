# scrnabox: A pipeline for scRNA analysis under HPC  
This repository includes a tutorial for the scRNA-seq data analysis using scrnabox pipeline under an HPC environment ([slurm work load manager system](https://slurm.schedmd.com/)). 

## Contents
- [Workflow analysis](#workflow-analysis)
  - [Standard scRNA-seq]
  - [Cell hashtags]
- [scrnabox.slurm](#scrnaboxsvn)
- [scrnaboxR](#scrnaboxr)
- [Processed Data](#processed-data)
- [Tutorial](#tutorial)
- [References](#references)

<details id=0>
<summary>
  
## Workflow analysis
  
</summary>
  
The following figures illustrate the steps involved in analyzing scRNA-seq data using the Standard and Cell Hashtags with the scrnabox pipeline
<br />
<br />
<kbd>
![Steps of Standard scRNA-seq ](https://raw.githubusercontent.com/neurobioinfo/scrnabox/main/figs/scrna.png)
</kbd>
<kbd>
![Steps of Cell Hashtags scRNA-seq](https://raw.githubusercontent.com/neurobioinfo/scrnabox/main/figs/hto.png)
</kbd>

</details>
  
This pipeline currently includes implementation of the standard and cell hashtags.

<details id=0>
<summary>
  
#### [Standard scRNA-seq](https://github.com/neurobioinfo/scrnabox/tree/main/README_SCRNA.md)
</summary>
The following steps describe how to analyze scRNA-seq data using the pipeline:
- Step 1: cellranger - This step runs Cellranger on the scRNA-seq data to generate the feature-barcode matrices for each sample.
- Step 2: Seurat object - This step runs Seurat on the feature-barcode matrices obtained from step 1 to generate a Seurat object for each sample. The Seurat object contains a lot of information, and accessing SeuratObject@meta.data will provide a data frame with relevant information on each cell.
- Step 3: QC and filter - The Seurat object includes quality metrics that can be used to filter cells and genes against possible doublets. Metrics such as total UMI counts per cell (nCount_RNA), total number of detected features per cell (nFeature_RNA), and mitochondrial count (percent.mito) are often used.
- Step 4: Doublet removal - This step can be used to remove doublets from the data. By default, the pipeline removes doublets, but you can choose to keep them by changing the parameter from 'yes' to 'no'.
- Step 5: Integration - This step integrates multiple scRNA-seq datasets using the Comprehensive Integration of Single Cell Data (CCA) method in Seurat, Tim, et al. (2019). The pipeline identifies anchors using the FindIntegrationAnchors function and passes them to the IntegrateData function to get a single Seurat object representing all the datasets.


- Step 6: Clustering- It involves clustering the data using a k-nearest neighbor graph based on the integrated PCA. This step produces UMAP and heatmaps of unlabelled clusters. However, it is up to the user to decide on the best cluster resolution outside the pipeline by examining the output and selecting the most appropriate annotation for the clusters.
- step 7: Cluster annotate - In this step, you can annotate the clusters with known cell types or use marker genes to predict the cell type of each cluster. 
- step 8: Differential gene expression (DEG)- There are multiple ways to perform differential gene expression analysis, but in this pipeline, we use the FindAllMarkers function to rank the highly differentially expressed genes in each cluster, which allows us to identify genes that are significantly differentially expressed between each cluster and the rest of the cells. From there, we can define contrasts to run statistical tests and investigate the phenotype and genotypes of each cluster.

The Step 1 - Step 8 can be done using [scrnabox.slurm](https://github.com/neurobioinfo/scrnabox/tree/main/README_SCRNA.md) in the HPC system ([slurm work load manager system](https://slurm.schedmd.com/)).
- step 9, Enrichment analysis: in this step, we obtain a list of significant genes using enrichment methods. The step 9 can be done using scrnaboxR, see [Practice](https://github.com/neurobioinfo/scrnabox/blob/main/tutorial/practice.md).
</details>
<details id=1>
<summary>
  
#### [Cell Hashtags](https://github.com/neurobioinfo/scrnabox/tree/main/README_HTO.md)
  
</summary>
  
The following steps explain how to analyze the  Hashtag oligonucleotide (henceforth referred to as HTO)
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
</details>

## [scrnabox.slurm](https://github.com/neurobioinfo/scrnabox/tree/main/scrnabox.slurm)
`scrnabox.slurm` is a pipeline developed to analyse under HPC system ([slurm work load manager system](https://slurm.schedmd.com/)), the pipeline has been using under [Beluga](https://docs.alliancecan.ca/wiki/B%C3%A9luga), details on how to use it are discussed in the [scrnabox.slurm](https://github.com/neurobioinfo/scrnabox/tree/main/scrnabox.slurm). 

## [scrnaboxR](https://github.com/neurobioinfo/scrnabox/tree/main/scrnaboxR)
The `scrnaboxR` is an R package containg various functions for running enrichment analysis and other analyses related to single-cell RNA sequence data. . To install this package, you can use the following script:
```
devtools::install_github("neurobioinfo/scrnabox/scrnaboxR")
```

## [Processed Data](https://github.com/neurobioinfo/scrnabox/blob/main/README_PROC.md)
The pipeline can also be used with [processed data](https://github.com/neurobioinfo/scrnabox/blob/main/README_PROC.md) from different projects, allowing users.  

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

