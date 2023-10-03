# set common R library path
R_LIB_PATH="/home/fiorini9/scratch/scrnabox_final/4.2"
.libPaths(R_LIB_PATH)

# load packages
install.packages('caTools') 
library(caTools) 

install.packages('gplots')
library(gplots) 

install.packages('ROCR')
library(ROCR) 

install.packages('Seurat') 
library(Seurat) 

install.packages('ggplot2') 
library(ggplot2) 

install.packages('dplyr') 
library(dplyr) 

install.packages('foreach') 
library(foreach) 

install.packages('doParallel') 
library(doParallel) 

install.packages('Matrix') 
library(Matrix) 

install.packages('remotes') 
library(remotes) 

remotes::install_github('chris-mcginnis-ucsf/DoubletFinder') 
library(DoubletFinder) 

install.packages('cowplot') 
library(cowplot) 

install.packages('clustree') 
library(clustree) 

install.packages('xlsx') 
library(xlsx) 

install.packages('enrichR') 
library(enrichR) 

install.packages('stringi') 
library(stringi) 

if (!requireNamespace("BiocManager", quietly = TRUE)) 
    install.packages("BiocManager")

BiocManager::install("limma") ## might be able to remove this
library(limma) 

install.packages('tidyverse') 
library(tidyverse) 

BiocManager::install("edgeR", force = TRUE)  ## might be able to remove this
library(edgeR) 

install.packages('vctrs') 
library(vctrs) 

install.packages('fossil') 
library(fossil) 

install.packages('openxlsx') 
library(openxlsx) 

install.packages('stringr') 
library(stringr) 

install.packages('ggpubr') 
library(ggpubr) 

install.packages('devtools') 
library(devtools) 

library(data.table)

install.packages('SoupX') 
library(SoupX) 

BiocManager::install("MatrixGenerics") 
library(MatrixGenerics) 

BiocManager::install("BiocGenerics") 
library(BiocGenerics) 

BiocManager::install("S4Vectors") 
library(S4Vectors) 

BiocManager::install("IRanges") 
library(IRanges) 

BiocManager::install("GenomeInfoDb") 
library(GenomeInfoDb) 

BiocManager::install("GenomicRanges") 
library(GenomicRanges) 

BiocManager::install("Biobase") 
library(Biobase) 

BiocManager::install("SummarizedExperiment") 
library(SummarizedExperiment) 

BiocManager::install("SingleCellExperiment") 
library(SingleCellExperiment) 

BiocManager::install("DropletUtils") 
library(DropletUtils)  
