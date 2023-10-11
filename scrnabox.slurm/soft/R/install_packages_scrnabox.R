#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
r_lib_path=args[1]
## Specify the path 
.libPaths(r_lib_path) 
## Where are my packages saved?
# .libPaths()
# R package requirements ----------------------------------------------------
#  if you are using the unix, make sure you already installed: build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev
## Installl devtool
install.packages("devtools", repos='https://cran.r-project.org', dependencies=TRUE)
if (! library("devtools", character.only = TRUE, logical.return = TRUE)) {
    write(paste0("Installation of package ", "devtools" ," exited with non-zero exit status"),stdout())
}


## install via Github
write("Installation of packages via Github",  stdout())

devtools::install_github("igraph/rigraph")
if (! library("igraph", character.only = TRUE, logical.return = TRUE)) {
    write(paste0("Installation of package ", "ggplot2" ," exited with non-zero exit status"),stdout())
}
devtools::install_github("tidyverse/tidyverse",force = TRUE)
if (! library("tidyverse", character.only = TRUE, logical.return = TRUE)) {
    write(paste0("Installation of package ", "tidyverse" ," exited with non-zero exit status"),stdout())
}

devtools::install_github('chris-mcginnis-ucsf/DoubletFinder',force = TRUE)
if (! library("DoubletFinder", character.only = TRUE, logical.return = TRUE)) {
    write(paste0("Installation of package ", "DoubletFinder" ," exited with non-zero exit status"),stdout())
}


## install via Cran
write("Installation of packages via cran",  stdout())

# Seurat need modules of gsl, hdf5, gcc
pkgs_seurat<-'Seurat'
install.packages(pkg, repos='https://cran.r-project.org', dependencies=TRUE)


if (! library(pkgs_seurat, character.only = TRUE, logical.return = TRUE)) {
  write(paste0("Installation of package ", pkgs_seurat ," exited with non-zero exit status, you should install it manually"),
        stdout())
  quit(status = 1, save = "no")
}



pkgs_cran<-c('tidyverse', 'foreach', 'doParallel', 'Seurat', 'caTools','gplots','ROCR','vctrs', 'fossil', 'openxlsx',
'stringr','ggpubr','data.table','ggrepel', 'scCustomize')

for (pkg in (pkgs_cran)) {
  install.packages(pkg, repos='https://cran.r-project.org', dependencies=TRUE)
  if (! library(pkg, character.only = TRUE, logical.return = TRUE)) {
    write(paste0("Installation of package ", pkg ," exited with non-zero exit status"),
          stdout())
    quit(status = 1, save = "no")
  }
}


## install Bioconductor 
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager",repos='https://cran.r-project.org', dependencies=TRUE)
}

pkgs_BiocManager <- c('Matrix','cowplot','xlsx','clustree','enrichR','stringi','limma', 'edgeR', 
'org.Hs.eg.db', 'SoupX','MatrixGenerics','BiocGenerics','S4Vectors','IRanges','GenomeInfoDb','GenomicRanges', 'Biobase', 
'SummarizedExperiment', 'SingleCellExperiment', 'DropletUtils', 'EnhancedVolcano')

## install via BiocManager
write("Installation of package via BiocManager",  stdout())
for (pkg in basename(pkgs_BiocManager)) {
  BiocManager::install(pkg, ask = FALSE, update = FALSE)
  if (! library(pkg, character.only = TRUE, logical.return = TRUE)) {
    write(paste0("Installation of package ",pkg," exited with non-zero exit status"),  stdout())
    quit(status = 1, save = "no")
  }
}

## Install scrnabox
devtools::install_github("neurobioinfo/scrnabox/scrnaboxR")

if (! library("scrnaboxR", character.only = TRUE, logical.return = TRUE)) {
    write(paste0("Installation of package ", "scrnaboxR" ," exited with non-zero exit status"),stdout())
}
