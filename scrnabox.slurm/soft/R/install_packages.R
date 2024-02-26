#!/usr/bin/env Rscript

# Path of R library ----------------------------------------------------------
## Create if it does not exist
if (!dir.exists("~/scratch/R/x86_64-pc-linux-gnu-library/4.2")){
  dir.create("~/scratch/R/x86_64-pc-linux-gnu-library/4.2", recursive = TRUE)
}

## Specify the path 
R_LIB_PATH="~/scratch/R/x86_64-pc-linux-gnu-library/4.2";.libPaths(R_LIB_PATH) 
## Where are my packages saved?
.libPaths()
# R package requirements ----------------------------------------------------
## Installl devtool, if you are using the unix, make sure you already installed: build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev
install.packages("devtools", dependencies = TRUE)
require(devtools)
install_version("Matrix", version = "1.6.0", repos = "http://cran.us.r-project.org")
install_version("igraph", version = "1.5.0.1", repos = "http://cran.us.r-project.org")
install_version("ggplot2", version = "3.4.2", repos = "http://cran.us.r-project.org")
install_version("dplyr", version = "1.1.2", repos = "http://cran.us.r-project.org")
install_version("foreach", version = "1.5.2", repos = "http://cran.us.r-project.org")
install_version("doParallel", version = "1.0.17", repos = "http://cran.us.r-project.org")
install_version("tidyverse", version ="2.0.0", repos = "http://cran.us.r-project.org", dependencies=FALSE)
install_version("Seurat", "4.3.0.1", repos = "http://cran.us.r-project.org", dependencies=FALSE)
install_version("SeuratObject", "4.1.3", repos = "http://cran.us.r-project.org", dependencies=FALSE)
install_github('https://github.com/ekernf01/DoubletFinder', force = T)
####################
install_version("ggpubr" , version ="0.6.0", repos = "http://cran.us.r-project.org",dependencies=FALSE)
install_version("SoupX" , version ="1.6.2", repos = "http://cran.us.r-project.org",dependencies=FALSE)
install_version("xlsx" , version ="0.6.5", repos = "http://cran.us.r-project.org",dependencies=FALSE)
install_version("data.table" , version ="1.14.8", repos = "http://cran.us.r-project.org",dependencies=FALSE)
install.packages("BiocManager",version ="1.30.21.1",repos='https://cran.r-project.org')

### packages needed for scrnabox pipeline
pkgs <- c('cowplot','xlsx','clustree','enrichR','stringi','limma', 'edgeR', 'DESeq2','org.Hs.eg.db', 'SoupX','MatrixGenerics','BiocGenerics','S4Vectors','IRanges','GenomeInfoDb','GenomicRanges', 'Biobase', 'SummarizedExperiment', 'SingleCellExperiment', 'DropletUtils', 'EnhancedVolcano','MAST')

## install Bioconductor 
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager",version ="1.30.21.1",repos='https://cran.r-project.org', dependencies=FALSE)
}

## install and check package loading
for (pkg in basename(pkgs)) {
  BiocManager::install(pkg, ask = FALSE, update = FALSE)
  if (! library(pkg, character.only = TRUE, logical.return = TRUE)) {
    write(paste0("Installation of package ",
                 pkg,
                 " exited with non-zero exit status"),
          stdout())
    quit(status = 1, save = "no")
  }
}

## Install scrnabox
devtools::install_github("neurobioinfo/scrnabox/scrnaboxR",dependencies=FALSE)
