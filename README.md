# scrnabox: A pipeline for scRNA 
This repository includes

## Contents
- [Workflow analysis]
- [scrnabox.svn]
- [scrnaboxR]
- [Rmarkdown]

## Analysis workflow
The analysis can be done via 

#### Hashtag
The following figure shows the steps to analyze the hashtag scRNA
![hashtag](https://raw.githubusercontent.com/neurobioinfo/scrnabox/main/hashtag.png)

- Step 1, cellranger:  this step run cellranger  on the single cell data. 
- Step 2, Seurat object: this step run seurat on feature-barcode matrices obtained from step 1 to generate seurat object for each sample.  Seurat's object includes a lot of information 
SeuratObject@meta.data will return the data frame and relevant information on each cell. 
- Step 3, QC and filter:  Seurat object includes some quality measure that can be used to filter cell and genes against possible doublets, we often use the total UMI counts per cell (nCount_RNA), the total number of detected features per cell (nFeature_RNA), and  mitochondrial count (percent.mito). 
- Step 4, demuplixing: 
- Step 5, integration: this step integrates multiple single cell RNA-seq datasets.  Seurat uses the Comprehensive Integration of Single Cell Data (CCA), see Stuart, Tim, et al. (2019) to perform integration; we identify anchors using the FindIntegrationAnchors function and pass them to the IntegrateData function to get a Seurat object.
- Step 6, Clustering: here, we run clustering (a k-nearest neighbour graph) on the intergrated PCA. 
- step 7,  Differetial gene expression (DEG):  DEG can be done in different ways, here we  run the function FindAllMarkers to compute a ranking for the highly differential genes in each cluster which determines the genes differentially expressed between each cluster and the rest of the cells. Then define contrast to run statstical tests to study the phenotype and genotypes. 
- step 8, Enrichment analysis: in this step, we obtain list of significant genes using .  


The step 1 - Step 7 can be done using `scrnabox.svn` and step 8 scrnaboxR. 

#### Non-hashtaq

-----------
## scrnabox.svn
`scrnabox.svn` is a pipeline developed to run step 1 to step 7 under HPC system, the pipeline has been using under [Beluga](https://docs.alliancecan.ca/wiki/B%C3%A9luga), the detail of how use it in [scrnabox.svn](https://github.com/neurobioinfo/scrnabox/blob/main/scrnabox_svn.md)


## scrnaboxR
The R package includes some functions to ????. The library can be downloaded using the following script. 
```
devtools::install_github("neurobioinfo/scrnabox/scrnaboxR")
```

## Rmarkdown
To run Enrichment Analysis look at [EA](https://github.com/neurobioinfo/scrnabox/blob/main/scrnaboxrmd/rmd1.md)


## References
- Stuart, T., Butler, A., Hoffman, P., Hafemeister, C., Papalexi, E., Mauck III, W. M., ... & Satija, R. (2019). Comprehensive integration of single-cell data. Cell, 177(7), 1888-1902.

## Contributing
This is an early version, any contribute or suggestion is appreciated, you can directly contact with [Saeid Amiri] or [Rhalena Thompson] 
## Changelog
Every release is documented on the [GitHub Releases page](https://github.com/neurobioinfo/scrnabox/releases).
## License
This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/neurobioinfo/scrnabox/blob/main/LICENSE) file for details
## Acknowledgement
## Todo

**[⬆ back to top](#contents)**
