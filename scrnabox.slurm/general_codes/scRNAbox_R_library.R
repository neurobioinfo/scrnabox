# set common R library path
R_LIB_PATH="/home/fiorini9/scratch/scrnabox_final/4.2"
.libPaths(R_LIB_PATH)

# load packages
install.packages('caTools') #in
library(caTools) 

install.packages('gplots') #in
library(gplots) 

install.packages('ROCR') #in
library(ROCR) 

install.packages('Seurat') #already in
library(Seurat) 

install.packages('ggplot2') #already in
library(ggplot2) 

install.packages('dplyr') #already in
library(dplyr) 

install.packages('foreach') #already in 
library(foreach) 

install.packages('doParallel') #already in 
library(doParallel) 

install.packages('Matrix') #already in 
library(Matrix) 

install.packages('remotes') #dont need
library(remotes) 

remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')  #already in 
library(DoubletFinder) 

install.packages('cowplot') #already in 
library(cowplot) 

install.packages('clustree') 
library(clustree) #already in 

install.packages('xlsx') #already in 
library(xlsx) 

install.packages('enrichR') #already in 
library(enrichR) 

install.packages('stringi') #already in 
library(stringi) 

if (!requireNamespace("BiocManager", quietly = TRUE)) #already in 
    install.packages("BiocManager")

BiocManager::install("limma") ## might be able to remove this
library(limma) 

install.packages('tidyverse') #in 
library(tidyverse) 

BiocManager::install("edgeR", force = TRUE)  ## might be able to remove this #in 
library(edgeR) 

install.packages('vctrs') #added
library(vctrs) 

install.packages('fossil') #added
library(fossil) 

install.packages('openxlsx') #added
library(openxlsx) 

install.packages('stringr') #added
library(stringr) 

install.packages('ggpubr') #added
library(ggpubr) 

install.packages('devtools') #already in 
library(devtools) 

library(data.table) #added

install.packages('SoupX') #added
library(SoupX) 

BiocManager::install("MatrixGenerics") #in
library(MatrixGenerics) 

BiocManager::install("BiocGenerics") #in
library(BiocGenerics) 

BiocManager::install("S4Vectors") #in
library(S4Vectors) 

BiocManager::install("IRanges") #in
library(IRanges) 

BiocManager::install("GenomeInfoDb") #in
library(GenomeInfoDb) 

BiocManager::install("GenomicRanges") #in
library(GenomicRanges) 

BiocManager::install("Biobase") #added
library(Biobase) 

BiocManager::install("SummarizedExperiment") #added
library(SummarizedExperiment) 

BiocManager::install("SingleCellExperiment") #added
library(SingleCellExperiment) 

BiocManager::install("DropletUtils") #added
library(DropletUtils)  
