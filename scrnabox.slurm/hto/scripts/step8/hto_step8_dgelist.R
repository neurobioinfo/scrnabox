#!/usr/bin/env Rscript

####################
# step8 -- DGE contrast
####################

## load parameters
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]

## load library
.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr','stringi','tidyverse','Matrix', 'ggrepel')
lapply(packages, library, character.only = TRUE)

## load parameters
source(paste(output_dir,'/job_info/parameters/step8_par.txt',sep=""))

############################################################################
## add metadata
############################################################################

if (tolower(par_run_add_metadata)=='yes') {
## load step7 seurat object
sample_name<-list.files(path = paste(output_dir, "/step7/objs7",sep=""),pattern = "*.rds")
seu_int<-readRDS(paste(output_dir,'/step7/objs7/',sample_name, sep=''))

## load meatadata
meta_data <- read.delim(par_metadata, header = T, sep = ",") 
new_meta <- colnames(meta_data)

## create existing metdata dataframe
metadata_df <- data.frame(seu_int@meta.data)

## merge existing metdata and new metdata
merge_metadata_df <- merge(metadata_df,meta_data, by = "orig.ident" )
df<- merge_metadata_df[,new_meta]
nrow(metadata_df) == nrow(df)

## remove orig.ident column
df <- subset(df, select = -c(orig.ident))
nrow(metadata_df) == nrow(df)

## set colnames into a list
colnames <- colnames(df)
nrow(df)

## add metadata
for (i in 1:ncol(df)){
seu_int <- AddMetaData(seu_int, metadata = df[,i], col.name = colnames[i])
}

## print RDS object
saveRDS(seu_int, paste(output_dir,'/step8/objs8',"/seu_step8.rds", sep=""))

## save rna expression matrix
if (tolower(par_save_RNA)=='yes') {
    mat <- GetAssayData(object = seu_int, assay = "RNA", slot = "data")
    #  write.csv(mat, paste(output_dir,'/step2/info2/',sample_name[i],"_RNA.csv", sep=""))
    writeMM(mat,file= paste(output_dir,'/step8/info8/seu',"_RNA.txt", sep=""))
}

## save metadata dataframe
if (tolower(par_save_metadata)=='yes') {
    write.csv(seu_int[[]], file = paste(output_dir,'/step8/info8/seu_MetaData.txt', sep=""), quote = TRUE, sep = ",")
}

## write session info
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step8/info8/sessionInfo_add_metadata.txt', sep=""))

if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}
}


############################################################################
## sample-sample contrast (wilcoxon)
############################################################################

if (tolower(par_run_sample_sample_wilcoxon)=='yes') {

################### import the right Seurat object ###################
## load name of existing Seurat objects
sample_name<-list.files(path = paste(output_dir, "/step7/objs7",sep=""),pattern = "*.rds")

if(file.exists(paste(output_dir,'/step8/objs8/','seu_step8.rds', sep = ""))){
    seu_int<-readRDS(paste(output_dir,'/step8/objs8/','seu_step8.rds', sep=''))
}else{
    seu_int<-readRDS(paste(output_dir,'/step7/objs7/',sample_name, sep=''))
}
################### ############################## ###################

## add info directory
OUT_DIR_info <- paste(output_dir,"/step8/info8",sep='') 
OUT_dir_info_sample <- paste(OUT_DIR_info,"/sample_sample_contrasts/",sep='') 
dir.create(OUT_dir_info_sample)

## add figs directory
OUT_DIR_figs <- paste(output_dir,"/step8/figs8",sep='') 
dir.create(OUT_DIR_figs)
OUT_DIR_figs_sample <- paste(OUT_DIR_figs,"/sample_sample_contrasts/",sep='') 
dir.create(OUT_DIR_figs_sample)

## set default assay to RNA
DefaultAssay(seu_int) <- "RNA"

### sample-sample dge wilcoxon analysis ####
dd<-read.csv(paste(output_dir,'/job_info/parameters/step8_contrast_genotype.txt',sep='/'), sep="")

for(i in 1:nrow(dd)){  
Idents(seu_int) <- dd[i,2]    
DGE <- FindMarkers(seu_int, ident.1 = dd[i,3], ident.2 = dd[i,4],  logfc.threshold = 0)
#write dge
write.csv(DGE, file = paste(OUT_dir_info_sample, dd[i,1],'_DEG.csv', sep=""), quote = FALSE, sep = ",")

# volcano plot
DGE_frame <- data.frame(DGE)
DGE_frame$gene <- rownames(DGE_frame) 
DGE_frame$col <- "grey66"
DGE_frame$col[DGE_frame$avg_log2FC >= 0.25 & DGE_frame$p_val < 0.05] <- "indianred3"
DGE_frame$col[DGE_frame$avg_log2FC <= -0.25 & DGE_frame$p_val < 0.05] <- "dodgerblue"
#plot
ggplot(DGE_frame, aes(x = avg_log2FC, y =-log10(p_val), col = col)) + 
theme_classic() + 
geom_hline(yintercept = -log10(0.05), linetype="dashed", colour = "grey") +
geom_vline(xintercept = 0.25, linetype="dashed", colour = "grey") +
geom_vline(xintercept = -0.25, linetype="dashed", colour = "grey") +
geom_point() +  
scale_colour_identity() +
xlab('Log2(FC)') +
ylab('-Log10(p-value)') 
ggsave(file = paste(OUT_DIR_figs_sample,dd[i,1],'_volcano_plot.pdf', sep=''))
}
###################

## print RDS object
saveRDS(seu_int, paste(output_dir,'/step8/objs8',"/seu_step8.rds", sep=""))

## save rna expression matrix
if (tolower(par_save_RNA)=='yes') {
    mat <- GetAssayData(object = seu_int, assay = "RNA", slot = "data")
    #  write.csv(mat, paste(output_dir,'/step2/info2/',sample_name[i],"_RNA.csv", sep=""))
    writeMM(mat,file= paste(output_dir,'/step8/info8/seu',"_RNA.txt", sep=""))
}

## save metadata dataframe
if (tolower(par_save_metadata)=='yes') {
    write.csv(seu_int[[]], file = paste(output_dir,'/step8/info8/seu_MetaData.txt', sep=""), quote = TRUE, sep = ",")
}

## write session info
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step8/info8/sessionInfo_sample_sample_cont.txt', sep=""))
if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}
}

############################################################################
## sample-cell contrast (wilcoxon)
############################################################################

if (tolower(par_run_sample_cell_wilcoxon)=='yes') {

################### import the right Seurat object ###################
## load name of existing Seurat objects
sample_name<-list.files(path = paste(output_dir, "/step7/objs7",sep=""),pattern = "*.rds")

if(file.exists(paste(output_dir,'/step8/objs8/','seu_step8.rds', sep = ""))){
    seu_int<-readRDS(paste(output_dir,'/step8/objs8/','seu_step8.rds', sep=''))
}else{
    seu_int<-readRDS(paste(output_dir,'/step7/objs7/',sample_name, sep=''))
}
################### ############################## ###################

## add info directory
OUT_DIR_info <- paste(output_dir,"/step8/info8",sep='') 
OUT_dir_info_sample <- paste(OUT_DIR_info,"/sample_cell_contrasts/",sep='') 
dir.create(OUT_dir_info_sample)

## add figs directory
OUT_DIR_figs <- paste(output_dir,"/step8/figs8",sep='')
dir.create(OUT_DIR_figs)
OUT_DIR_figs_sample <- paste(OUT_DIR_figs,"/sample_cell_contrasts/",sep='') 
dir.create(OUT_DIR_figs_sample)

## set default assay to RNA
DefaultAssay(seu_int) <- "RNA"

################### sample-cell dge analysis ###################
dd<-read.csv(paste(output_dir,'/job_info/parameters/step8_contrast_celltype.txt',sep='/'), sep="")

for(i in 1:nrow(dd)){
Idents(seu_int) <- dd[i,2]    
celltype.sub.seu <- subset(seu_int, idents = dd[i,3])
Idents(celltype.sub.seu) <- dd[i,4]    
DGE <- FindMarkers(celltype.sub.seu, ident.1 = dd[i,5], ident.2 = dd[i,6],  logfc.threshold = 0)
write.csv(DGE, file = paste(OUT_dir_info_sample,"/", dd[i,1],'_DEG.csv', sep=""), quote = FALSE, sep = ",")

# volcano plot
DGE_frame <- data.frame(DGE)
DGE_frame$gene <- rownames(DGE_frame) 
DGE_frame$col <- "grey66"
DGE_frame$col[DGE_frame$avg_log2FC >= 0.25 & DGE_frame$p_val < 0.05] <- "indianred3"
DGE_frame$col[DGE_frame$avg_log2FC <= -0.25 & DGE_frame$p_val < 0.05] <- "dodgerblue"
#plot
ggplot(DGE_frame, aes(x = avg_log2FC, y =-log10(p_val), col = col)) + 
theme_classic() + 
geom_hline(yintercept = -log10(0.05), linetype="dashed", colour = "grey") +
geom_vline(xintercept = 0.25, linetype="dashed", colour = "grey") +
geom_vline(xintercept = -0.25, linetype="dashed", colour = "grey") +
geom_point() +  
scale_colour_identity() +
xlab('Log2(FC)') +
ylab('-Log10(p-value)')
ggsave(file = paste(OUT_DIR_figs_sample,dd[i,1],'_volcano_plot.pdf', sep=''))

}
################### ########################## ###################

## print RDS object
saveRDS(seu_int, paste(output_dir,'/step8/objs8',"/seu_step8.rds", sep=""))

## save rna expression matrix
if (tolower(par_save_RNA)=='yes') {
    mat <- GetAssayData(object = seu_int, assay = "RNA", slot = "data")
    #  write.csv(mat, paste(output_dir,'/step2/info2/',sample_name[i],"_RNA.csv", sep=""))
    writeMM(mat,file= paste(output_dir,'/step8/info8/seu',"_RNA.txt", sep=""))
}

## save metadata dataframe
if (tolower(par_save_metadata)=='yes') {
    write.csv(seu_int[[]], file = paste(output_dir,'/step8/info8/seu_MetaData.txt', sep=""), quote = TRUE, sep = ",")
}

## write session info
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step8/info8/sessionInfo_sample_cell_cont.txt', sep=""))
if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}
}

