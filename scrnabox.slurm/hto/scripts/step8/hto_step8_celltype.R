#!/usr/bin/env Rscript

####################
# step8 -- genotype-cell type contrast
####################

## load parameters
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
contrast=args[3]

## load library
.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr','stringi','limma','tidyverse')
lapply(packages, library, character.only = TRUE)

## load contrast matrix
if (contrast=="F") {
   dd<-read.csv(paste(output_dir,'/job_info/parameters/step8_contrast_celltype.txt',sep='/'), sep="")
} else {
   dd<-read.csv(paste(contrast), sep="")
}

## load Seurat object
seu.int.c<-readRDS(paste(output_dir,'/step8/objs8',"/seu_step8.rds", sep=""))

## load dge object
dge<-readRDS(paste(output_dir,'/step8/info8',"/dge.rds", sep=""))

## function to perform genotype-celltype DGE contrast
run_func_celltype<- function(design,des,maincell,pheno,selvar,gp2){
    xcontrol<-paste(pheno,maincell, sep='') 
    xrest<-paste('ALLNEED',maincell, sep='') 
    aaa<-colnames(design)
    aaa<-stri_replace_all_regex(aaa,pattern=c('d1', 'ct', ':', '-'),replacement=c(''),vectorize=FALSE)
    colnames(design)<-stri_replace_all_regex(aaa,pattern=gp2[which(gp2%in%selvar)], ,replacement=c('ALLNEED'),vectorize=FALSE)
    abc0<-c(paste(xrest, xcontrol, sep=' - '))
    abc0
    contrast.mat <- limma::makeContrasts(contrasts=abc0, levels = design)
    vm <- limma::voom(dge, design = design, plot = TRUE)
    fit <- limma::lmFit(vm, design = design)
    fit.contrasts <- limma::contrasts.fit(fit, contrast.mat)
    fit.contrasts <- limma::eBayes(fit.contrasts)
    limma::topTable(fit.contrasts, number = 10, sort.by = "P")
    limma_de <- limma::topTable(fit.contrasts, number = Inf, sort.by = "P")
    aa<-as.data.frame(limma_de)
   write.csv(aa,file=paste(output_dir,paste('/step8/cont_celltype',paste(des,'.csv',sep=''),sep='/'),sep=''))
}

## set up user-defined DGE contrast design
gp2<-unique(seu.int.c$MULTI_ID_Lables)
gp2<-stri_replace_all_regex(gp2,pattern=c(':', '-'),replacement=c(''),vectorize=FALSE)
design <- model.matrix(~ 0 + d1:ct , data = seu.int.c@meta.data)

## perform DGE contrast
for(i in 1:dim(dd)[1]){
    des<-str_split(dd[i,1],",")[[1]]
    maincell<-str_split(dd[i,2],",")[[1]]
    pheno<-str_split(dd[i,3],",")[[1]]
    pheno<-stri_replace_all_regex(pheno,pattern=c(':', '-'),replacement=c(''),vectorize=FALSE)
    selvar<-str_split(dd[i,4],",")[[1]]
    selvar<-stri_replace_all_regex(selvar,pattern=c(':', '-'),replacement=c(''),vectorize=FALSE)
    run_func_celltype(design,des,maincell,pheno,selvar,gp2)
}

## save session informationd
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step8/info8/sessionInfo_cont_celltype.txt', sep=""))