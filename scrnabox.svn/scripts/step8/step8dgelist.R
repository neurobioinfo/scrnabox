#!/usr/bin/env Rscript

####################
# step7dgelist
## Create seurat object 

args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]

.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr')
lapply(packages, library, character.only = TRUE)

library('stringi')
library('limma')
library('tidyverse')

seu.int.c<-readRDS(paste(output_dir,'/step6/objs','/seu_int_clu.rds', sep=''))


# if (contrast=="F") {
   # dd<-read.csv(paste(output_dir,'/job_output/step7_contrast_inte.txt',sep='/'), sep="")
# } else {
   # dd<-read.csv(paste(contrast), sep="")
# }


clus_label<-read.csv(paste(output_dir,'/job_output/step8_clus_label.txt',sep=''), header=FALSE, sep="")
cluster.ids<-str_split(clus_label[1,],",")[[1]]
# cluster.ids <- c("DAneurons","DAneurons","DAneurons","other","neurons","progenitors","neurons","neurons","neurons","ependymal","DAneurons","DAneurons","progenitors","other","Glia","progenitors","neuroblasts","neurons","neurons","neurons","progenitors","ependymal","DAneurons","DAneurons","neurons","DAneurons","neurons","progenitors","Glia")
names(cluster.ids) <- levels(seu.int.c)    # where seu.q is your seurat object
seu.int.c <- RenameIdents(seu.int.c, cluster.ids) 
seu.int.c$cell.types.pool <- Idents(seu.int.c)   # cell.types is a new slot in metadata you can call this whatever you want
de_genes <- Seurat::FindAllMarkers(seu.int.c,  min.pct = 0.25,only.pos = TRUE)
DefaultAssay(seu.int.c) <- "RNA"
counts <- Seurat::GetAssayData(seu.int.c, slot = "counts")
counts <- counts[rowSums(counts) != 0,]
# dim(counts)
dge <- edgeR::DGEList(counts = counts)
dge <- edgeR::calcNormFactors(dge)  

seu.int.c$ct<-seu.int.c$cell.types.pool
seu.int.c$d1<-seu.int.c$Lables

saveRDS(seu.int.c, paste(output_dir,'/step8/objs',"/seu.int.c.rds", sep=""))
saveRDS(de_genes, paste(output_dir,'/step8/objs',"/de_genes.rds", sep=""))
saveRDS(dge, paste(output_dir,'/step8/objs',"/dge.rds", sep=""))



# run_func_inte<- function(design,des,maincell,pheno,selvar,gp2){
#     xcontrol<-paste(pheno,maincell, sep='') 
#     xrest<-paste('ALLNEED',maincell, sep='') 
#     aaa<-colnames(design)
#     aaa<-stri_replace_all_regex(aaa,pattern=c('d1', 'ct', ':', '-'),replacement=c(''),vectorize=FALSE)
#     colnames(design)<-stri_replace_all_regex(aaa,pattern=gp2[which(gp2%in%selvar)], ,replacement=c('ALLNEED'),vectorize=FALSE)
#     abc0<-c(paste(xrest, xcontrol, sep=' - '))
#     abc0
#     contrast.mat <- limma::makeContrasts(contrasts=abc0, levels = design)
#     vm <- limma::voom(dge, design = design, plot = TRUE)
#     fit <- limma::lmFit(vm, design = design)
#     fit.contrasts <- limma::contrasts.fit(fit, contrast.mat)
#     fit.contrasts <- limma::eBayes(fit.contrasts)
#     limma::topTable(fit.contrasts, number = 10, sort.by = "P")
#     limma_de <- limma::topTable(fit.contrasts, number = Inf, sort.by = "P")
#     aa<-as.data.frame(limma_de)
#    write.csv(aa,file=paste(output_dir,paste('/step7/cont_inte',paste(des,'.csv',sep=''),sep='/'),sep=''))
# }

# gp2<-unique(seu.int.c$Lables)
# gp2<-stri_replace_all_regex(gp2,pattern=c(':', '-'),replacement=c(''),vectorize=FALSE)
# design <- model.matrix(~ 0 + d1:ct , data = seu.int.c@meta.data)

# for(i in 1:dim(dd)[1]){
#     des<-str_split(dd[i,1],",")[[1]]
#     maincell<-str_split(dd[i,2],",")[[1]]
#     # control
#     pheno<-str_split(dd[i,3],",")[[1]]
#     pheno<-stri_replace_all_regex(pheno,pattern=c(':', '-'),replacement=c(''),vectorize=FALSE)
#     # ex_control
#     selvar<-str_split(dd[i,4],",")[[1]]
#     selvar<-stri_replace_all_regex(selvar,pattern=c(':', '-'),replacement=c(''),vectorize=FALSE)
#     # selvar<-selvar[1]
#     # all 
#    #  gp2<-str_split(dd[i,5],",")[[1]]
#     run_func_inte(design,des,maincell,pheno,selvar,gp2)
# }


