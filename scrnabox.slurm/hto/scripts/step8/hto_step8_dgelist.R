#!/usr/bin/env Rscript

####################
# step 8 -- dge list
####################

## load parameters
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]

## load library
.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr','stringi','limma','tidyverse','edgeR')
lapply(packages, library, character.only = TRUE)

## load existing Seurat object
sample_name<-list.files(path = paste(output_dir, "/step6/objs6",sep=""),pattern = "*.rds")
seu.int.c<-readRDS(paste(output_dir,'/step6/objs6/',sample_name, sep=''))

## load parameters
source(paste(output_dir,'/job_info/parameters/step8_par.txt',sep=""))

## set new sample labels (optional)
if (tolower(par_new_genotype)=='yes') {
old.names<-par_old_antibody_label
new.names<-par_new_antibody_label
Idents(seu.int.c)  <- "MULTI_ID"
for (i in 1:length(old.names)){
    newIdent <- new.names[i]
    names(newIdent) <- old.names[i]
    seu.int.c <- RenameIdents(object = seu.int.c, newIdent)
}
seu.int.c[["MULTI_ID_Lables"]] <- Idents(seu.int.c)
}

## set user-defined cluster annotations
cluster.ids<-par_step8_clus_label

## set cluster level and identity
seu.int.c <- SetIdent(seu.int.c, value = par_select_cluster)
names(cluster.ids) <- levels(seu.int.c)    # where seu.q is your seurat object
seu.int.c <- RenameIdents(seu.int.c, cluster.ids) 
seu.int.c$cell.types.pool <- Idents(seu.int.c)   # cell.types is a new slot in metadata you can call this whatever you want

## compute differential gene expression on RNA assay
de_genes <- Seurat::FindAllMarkers(seu.int.c,  min.pct = 0.25,only.pos = TRUE)
DefaultAssay(seu.int.c) <- "RNA"
counts <- Seurat::GetAssayData(seu.int.c, slot = "counts")
counts <- counts[rowSums(counts) != 0,]
dge <- edgeR::DGEList(counts = counts)
dge <- edgeR::calcNormFactors(dge)  
seu.int.c$ct<-seu.int.c$cell.types.pool
seu.int.c$d1<-seu.int.c$MULTI_ID_Lables

## save Seurat RDS object
saveRDS(seu.int.c, paste(output_dir,'/step8/objs8',"/seu_step8.rds", sep=""))

## save DEG RDS objects
saveRDS(de_genes, paste(output_dir,'/step8/info8',"/de_genes.rds", sep=""))
saveRDS(dge, paste(output_dir,'/step8/info8',"/dge.rds", sep=""))

## save DEG information
write.csv(colnames(seu.int.c[[]]), file= paste(output_dir,'/step8/info8/meta_info_seu_step8.txt', sep=""))
write.csv(head(de_genes),          file= paste(output_dir,'/step8/info8/meta_info_de_genes.txt', sep=""))
write.csv(summary(dge),            file= paste(output_dir,'/step8/info8/meta_info_dge.txt', sep=""))

## save RNA expression matrix
if (tolower(par_save_RNA)=='yes') {
    mat <- GetAssayData(object = seu.int.c, assay = "RNA", slot = "data")
    writeMM(mat,file= paste(output_dir,'/step8/info8/seu',"_RNA.txt", sep=""))
}

## save metadata dataframe
if (tolower(par_save_metadata)=='yes') {
    write.csv(seu.int.c[[]], file = paste(output_dir,'/step8/info8/seu_MetaData.txt', sep=""), quote = TRUE, sep = ",")
}

##  save session information
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step8/info8/sessionInfo_dgelist.txt', sep=""))