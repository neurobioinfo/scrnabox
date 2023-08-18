# Installation
The package is written in the bash, so it can be used with any slurm system. To download  `scrnabox.slurm`, run the below comments 
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
The `scrnaboxR` is an R package that provides a collection of functions for conducting enrichment analysis and other analyses associated with single-cell RNA sequencing (scRNA-seq) data. It serves as a companion to scrnabox, offering a range of tools and functionalities to enhance scRNA-seq data analysis. You need to add the R info in `scrnabox_config.ini`, you can define the path of R library in `R_LIB_PATH=`, version of R in `R_VERSION`, you can add the path of `cell ranger`in `MODULECELLRANGER` 