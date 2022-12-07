# scrnabox: A pipeline for scRNA 
This repository includes

## Contents
- [Workflow analysis](#analysis-workflow)
- [scrnabox.svn](#scrnaboxsvn)
- [scrnaboxR](#scrnaboxr)
- [References](#references)


## Analysis workflow
#### Hashtag
The following figure shows the steps to analyze the hashtag scRNA
![hashtag](https://raw.githubusercontent.com/neurobioinfo/scrnabox/main/hashtag.png)

- Step 1, cellranger:  this step runs cellranger on the single cell data. 
- Step 2, Seurat object: this step runs seurat on feature-barcode matrices obtained from step 1 to generate seurat object for each sample.  Seurat's object includes a lot of information; SeuratObject@meta.data will return the data frame and relevant information on each cell. 
- Step 3, QC and filter: Seurat object includes some quality measure that can be used to filter cell and genes against possible doublets, we often use the total UMI counts per cell (nCount_RNA), the total number of detected features per cell (nFeature_RNA), and mitochondrial count (percent.mito). 
- Step 4, demuplixing: It process where use the barcode information in order to know which sequences came from which samples after they had all be sequenced together. Pipline has option to remove the doublets and negative cells as this step. 
- Step 5, integration: this step integrates multiple single cell RNA-seq datasets. Seurat uses the Comprehensive Integration of Single Cell Data (CCA), see Stuart, Tim, et al. (2019) to perform integration; we identify anchors using the FindIntegrationAnchors function and pass them to the IntegrateData function to get a Seurat object.
- Step 6, Clustering: here, we run clustering (a k-nearest neighbour graph) on the integrated PCA. At the end of this step we have some UMAP and heatmaps of unlabelled clusters.  We need to decide on the cluster annotation outside the pipeline. User needs to look at the output and decided on the best cluster resolution. 
- step 7, Cluster annotate: you can annotate the cluster with cell types in this step. 
- step 8, Differetial gene expression (DEG): DEG can be done in different ways, here we run the function FindAllMarkers to compute a ranking for the highly differential genes in each cluster which determines the genes differentially expressed between each cluster and the rest of the cells. Then define contrast to run statistical tests to study the phenotype and genotypes. 
- step 9, Enrichment analysis: in this step, we obtain list of significant genes using enrichment methods.  

The step 1 - Step 8 can be done using [scrnabox.svn](https://github.com/neurobioinfo/scrnabox/tree/main/scrnabox.svn), the step 9 should can done using scrnaboxR, see [Practice](https://github.com/neurobioinfo/scrnabox/blob/main/scrnaboxrmd/practice.md). 

#### Non-hashtaq
Ongoing 

-----------

## scrnabox.svn
`scrnabox.svn` is a pipeline developed to run step 1 to step 7 under HPC system, the pipeline has been using under [Beluga](https://docs.alliancecan.ca/wiki/B%C3%A9luga), the detail of how to use it is discussed in [scrnabox.svn](https://github.com/neurobioinfo/scrnabox/tree/main/scrnabox.svn)

## scrnaboxR
This R package includes some functions to the enrichment analysis. This library can be installed using the following script. 
```
devtools::install_github("neurobioinfo/scrnabox/scrnaboxR")
```

To analyze DGE look at [Practice](https://github.com/neurobioinfo/scrnabox/blob/main/scrnaboxrmd/practice.md)


## References
- Stuart, T., Butler, A., Hoffman, P., Hafemeister, C., Papalexi, E., Mauck III, W. M., ... & Satija, R. (2019). Comprehensive integration of single-cell data. Cell, 177(7), 1888-1902.

## Contributing
This is an early version, any contribute or suggestion is appreciated, you can directly contact with [Saeid Amiri](https://github.com/saeidamiri1) or [Rhalena Thomas](https://github.com/RhalenaThomas). 

## Citation
Amiri, S., Thomas, R., & Larivière, R. (2022). scrnabox: A pipeline for scRNA (Version 0.1.3) [Computer software]. https://github.com/neurobioinfo/scrnabox

## Changelog
Every release is documented on the [GitHub Releases page](https://github.com/neurobioinfo/scrnabox/releases).
## License
This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/neurobioinfo/scrnabox/blob/main/LICENSE) file for details
## Acknowledgement

The pipeline is done as part Dark Genome project, it is written by [Saeid Amiri](https://github.com/saeidamiri1) with associate of Rhalena Thomas and  Roxanne Larivière. 

## Todo
- Add nba to pipeline
**[⬆ back to top](#contents)**
