#!/usr/bin/env Rscript

####################
# step8 -- dgelist
####################

## load parameters
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]

## load library
.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr','stringi','limma','tidyverse','edgeR')
lapply(packages, library, character.only = TRUE)

## load existing Seurat objects
sample_name<-list.files(path = paste(output_dir, "/step6/objs6",sep=""),pattern = "*.rds")
seu.int.c<-readRDS(paste(output_dir,'/step6/objs6/',sample_name, sep=''))

## load parameters
source(paste(output_dir,'/job_info/parameters/step8_par.txt',sep=""))

## reset sample labels (optional)
if (tolower(par_new_genotype)=='yes') {
    old.names<-par_old_sample_label
    new.names<-par_new_sample_label
    Idents(seu.int.c)  <- "Sample_ID"
    for (i in 1:length(old.names)){
        newIdent <- new.names[i]
        names(newIdent) <- old.names[i]
        seu.int.c <- RenameIdents(object = seu.int.c, newIdent)
    }
    seu.int.c[["New_Sample_ID"]] <- Idents(seu.int.c)
}

## rename old sample labels metadata column name if users do not want to change labels -- this is necessary for DGE contrats
if (tolower(par_new_genotype)=='no') {
    Idents(seu.int.c)  <- "Sample_ID"
    seu.int.c[["New_Sample_ID"]] <- Idents(seu.int.c)
}

## set cluster annotations obtained from cluster anntation modules (step 7)
cluster.ids<-par_step8_clus_label

## set cluster resolution and rename cluster identities
seu.int.c <- SetIdent(seu.int.c, value = par_select_cluster)
names(cluster.ids) <- levels(seu.int.c)    
seu.int.c <- RenameIdents(seu.int.c, cluster.ids) 
seu.int.c$cell.types.pool <- Idents(seu.int.c) 

## print UMAP with final cluster annotation
dir.create(paste(output_dir,'/step8/figs8',sep=''))
DimPlot(seu.int.c, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(file = paste(output_dir,'/step8/figs8','/final_cluster_annotation.pdf', sep=''))

## identify cluster-specific marker genes
de_genes <- Seurat::FindAllMarkers(seu.int.c,  min.pct = 0.25,only.pos = TRUE)
DefaultAssay(seu.int.c) <- "RNA"
counts <- Seurat::GetAssayData(seu.int.c, slot = "counts")
counts <- counts[rowSums(counts) != 0,]

## compute differentially expressed genes (DEG)
dge <- edgeR::DGEList(counts = counts)
dge <- edgeR::calcNormFactors(dge)  

## set cell type
seu.int.c$ct<-seu.int.c$cell.types.pool

## set sample ID
seu.int.c$d1<-seu.int.c$New_Sample_ID

## save Seurat and DEG RDS objects
saveRDS(seu.int.c, paste(output_dir,'/step8/objs8',"/seu_step8.rds", sep=""))
saveRDS(de_genes, paste(output_dir,'/step8/info8',"/de_genes.rds", sep=""))
saveRDS(dge, paste(output_dir,'/step8/info8',"/dge.rds", sep=""))

## save DEG info
write.csv(colnames(seu.int.c[[]]), file= paste(output_dir,'/step8/info8/meta_info_seu_step8.txt', sep=""))
write.csv(head(de_genes),          file= paste(output_dir,'/step8/info8/meta_info_de_genes.txt', sep=""))
write.csv(summary(dge),            file= paste(output_dir,'/step8/info8/meta_info_dge.txt', sep=""))

## save rna expression matrix
if (tolower(par_save_RNA)=='yes') {
    mat <- GetAssayData(object = seu.int.c, assay = "RNA", slot = "data")
    #  write.csv(mat, paste(output_dir,'/step2/info2/',sample_name[i],"_RNA.csv", sep=""))
    writeMM(mat,file= paste(output_dir,'/step8/info8/seu',"_RNA.txt", sep=""))
}

## save metadata dataframe
if (tolower(par_save_metadata)=='yes') {
    write.csv(seu.int.c[[]], file = paste(output_dir,'/step8/info8/seu_MetaData.txt', sep=""), quote = TRUE, sep = ",")
}

## write session info
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step8/info8/sessionInfo_dgelist.txt', sep=""))
file.remove("Rplots.pdf")