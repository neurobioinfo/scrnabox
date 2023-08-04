#!/usr/bin/env Rscript

####################
# step7dgelist

args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]

.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr','stringi','limma','tidyverse','edgeR')
lapply(packages, library, character.only = TRUE)

# library('stringi')
# library('limma')
# library('tidyverse')

sample_name<-list.files(path = paste(output_dir, "/step6/objs6",sep=""),pattern = "*.rds")
seu.int.c<-readRDS(paste(output_dir,'/step6/objs6/',sample_name, sep=''))


source(paste(output_dir,'/job_info/parameters/step8_par.txt',sep=""))


# clus_label<-read.csv(paste(output_dir,'/job_info/parameters/step8_clus_label.txt',sep=''), header=FALSE, sep="")
# cluster.ids<-str_split(clus_label[1,],",")[[1]]

cluster.ids<-step8_clus_label

seu.int.c <- SetIdent(seu.int.c, value = select_cluster)
names(cluster.ids) <- levels(seu.int.c)    # where seu.q is your seurat object
seu.int.c <- RenameIdents(seu.int.c, cluster.ids) 
seu.int.c$cell.types.pool <- Idents(seu.int.c)   # cell.types is a new slot in metadata you can call this whatever you want
de_genes <- Seurat::FindAllMarkers(seu.int.c,  min.pct = 0.25,only.pos = TRUE)
DefaultAssay(seu.int.c) <- "RNA"
counts <- Seurat::GetAssayData(seu.int.c, slot = "counts")
counts <- counts[rowSums(counts) != 0,]

dge <- edgeR::DGEList(counts = counts)
dge <- edgeR::calcNormFactors(dge)  

seu.int.c$ct<-seu.int.c$cell.types.pool
seu.int.c$d1<-seu.int.c$MULTI_ID_Lables

saveRDS(seu.int.c, paste(output_dir,'/step8/objs8',"/seu_step8.rds", sep=""))
saveRDS(de_genes, paste(output_dir,'/step8/info8',"/de_genes.rds", sep=""))
saveRDS(dge, paste(output_dir,'/step8/info8',"/dge.rds", sep=""))

write.csv(colnames(seu.int.c[[]]), file= paste(output_dir,'/step8/info8/meta_info_seu_step8.txt', sep=""))
write.csv(head(de_genes),          file= paste(output_dir,'/step8/info8/meta_info_de_genes.txt', sep=""))
write.csv(summary(dge),            file= paste(output_dir,'/step8/info8/meta_info_dge.txt', sep=""))

if (tolower(Save_RNA)=='yes') {
    mat <- GetAssayData(object = seu.int.c, assay = "RNA", slot = "data")
    writeMM(mat,file= paste(output_dir,'/step8/info8/seu',"_RNA.txt", sep=""))
}

if (tolower(Save_metadata)=='yes') {
    write.csv(seu.int.c[[]], file = paste(output_dir,'/step8/info8/seu_MetaData.txt', sep=""), quote = TRUE, sep = ",")
}

writeLines(capture.output(sessionInfo()), paste(output_dir,'/step8/info8/sessionInfo_dgelist.txt', sep=""))
