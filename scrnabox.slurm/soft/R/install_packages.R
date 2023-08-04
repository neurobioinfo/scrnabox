#!/usr/bin/env Rscript

# Path of R library ----------------------------------------------------------
## Create if it does not exist
if (!dir.exists("~/scratch/R/x86_64-pc-linux-gnu-library/4.2")){
  dir.create("~/scratch/R/x86_64-pc-linux-gnu-library2/4.2", recursive = TRUE)
}

## Specify the path 
R_LIB_PATH="~/scratch/R/x86_64-pc-linux-gnu-library/4.2";.libPaths(R_LIB_PATH) 
## Where are my packages saved?
.libPaths()
# R package requirements ----------------------------------------------------
## Installl devtool, if you are using the unix, make sure you already installed: build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev
install.packages("devtools", dependencies = TRUE)
devtools::install_github("igraph/rigraph")
devtools::install_github("tidyverse/ggplot2")
devtools::install_github("tidyverse/dplyr")
devtools::install_github("tidyverse/tidyverse")
devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')
### packages needed for scrnabox pipeline
pkgs <- c('Seurat', 'foreach', 'doParallel', 'Matrix','cowplot','xlsx','clustree','enrichR','stringi','limma', 'edgeR', 'org.Hs.eg.db')
pkgs <- c('SoupX','MatrixGenerics','BiocGenerics','S4Vectors','IRanges','GenomeInfoDb','GenomicRanges')

## install Bioconductor 
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
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
devtools::install_github("neurobioinfo/scrnabox/scrnaboxR")
