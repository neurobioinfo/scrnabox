##########################################################################################
## Figure 3A
##########################################################################################
## load parameters
output_dir = "/path/to/output/directory"
r_lib_path = "/path/to/R/library"

## load library
.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr', 'foreach', 'doParallel','Matrix', 'ggpubr')
lapply(packages, library, character.only = TRUE)

## load parameters text file
source(paste(output_dir,'/job_info/parameters/step3_par.txt',sep=""))

## if user has existing Seurat object 
if (exists("par_seurat_object")) {                                                  
sample_name<-list.files(path = par_seurat_object)
sample_nameb<-gsub(".rds","",sample_name)
if(length(sample_name)<1) {
print("You do not have any existing Seurat object")
}
} else {
sample_name<-list.files(path = paste(output_dir, "/step2/objs2",sep=""))
sample_nameb<-gsub(".rds","",sample_name)
if(length(sample_name)<1) {
print("You do not have any object from step 2 ")
}
}

## create a list of available seurat objects in the 
for (i in 1:length(sample_name)) {
if (!grepl(".rds",tolower(sample_name[i]), fixed = TRUE)){
print(c(sample_name[i],"is not R rds"))
}
} 

## set seed for replicability
set.seed(1234)

## detect number of cores
numCores <- detectCores()
cl <- makeCluster(numCores-1)
registerDoParallel(cl) 

###### QC and filtering if users have a preporcessed Seurat object

i = 1  
set.seed(1234)
if (exists("par_seurat_object")) { 
seu<-readRDS(paste(par_seurat_object, "/", sample_name[i], sep=""))
} else  {
seu<-readRDS(paste(output_dir,'/step2/objs2/',sample_name[i], sep=""))
}
print(sample_name[i])

## calculate percent MT and Ribo if users are starting from Step 3 as they may not have this claculated in their Seurat object.
if (exists("par_seurat_object")) { 
## calculate percent MT 
seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")
## calculate percent ribo  
seu[["percent.ribo"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^RP[SL]")    #NEW CODE
}

## fiter according to user-defined thresholds.
if (exists("par_nFeature_RNA_L")) {
cat("nFeature_RNA Lower: \n")
print(par_nFeature_RNA_L)
seu <- subset(seu, subset = nFeature_RNA > par_nFeature_RNA_L) 
}
if (exists("par_nFeature_RNA_U")){
cat("nFeature_RNA Upper: \n")
print(par_nFeature_RNA_U)
seu <- subset(seu, subset = nFeature_RNA < par_nFeature_RNA_U) 
}
if (exists("par_nCount_RNA_L")){
cat("nCount_RNA Lower: \n")
print(par_nCount_RNA_L)
seu <- subset(seu, subset = nCount_RNA > par_nCount_RNA_L) 
}
if (exists("par_nCount_RNA_U")){
cat("nCount_RNA Upper: \n")
print(par_nCount_RNA_U)        
seu <- subset(seu, subset = nCount_RNA < par_nCount_RNA_U) 
}
if (exists("par_mitochondria_percent_L")) {
cat("Mitochondria_percent Lower: \n")
print(par_mitochondria_percent_L)        
seu <- subset(seu, subset = percent.mt > par_mitochondria_percent_L) 
}
if (exists("par_mitochondria_percent_U")){
cat("Mitochondria_percent Upper: \n")
print(par_mitochondria_percent_U)                
seu <- subset(seu, subset = percent.mt < par_mitochondria_percent_U) 
}
if (exists("par_ribosomal_percent_L")) {
cat("Ribosomal_percent Lower: \n")
print(par_ribosomal_percent_L)                        
seu <- subset(seu, subset = percent.ribo > par_ribosomal_percent_L) 
}
if (exists("par_ribosomal_percent_U")) {
cat("Ribosomal_percent Upper: \n")
print(par_ribosomal_percent_U)                                
seu <- subset(seu, subset = percent.ribo < par_ribosomal_percent_U) 
}

## optional: filter out mitochondrial genes
if (tolower(par_remove_mitochondrial_genes)=='yes') {
MT_genes <- grep( "^MT-", rownames(seu), value = T)
counts <- GetAssayData(seu, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% MT_genes)),]
seu[["RNA"]] <- subset(seu[["RNA"]], features = rownames(counts))     
}

## optional: filter out ribosomal genes
if (tolower(par_remove_ribosomal_genes)=='yes') {
Ribo_genes <- grep( "^RP[SL]", rownames(seu), value = T)
counts <- GetAssayData(seu, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% Ribo_genes)),]
seu[["RNA"]] <- subset(seu[["RNA"]], features = rownames(counts))       
}

## remove custom list of genes
if (exists("par_remove_genes")) {
counts <- GetAssayData(seu, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% par_remove_genes)),]
seu[["RNA"]] <- subset(seu[["RNA"]], features = rownames(counts))    
}

## normalize after filtering
seu <- Seurat::NormalizeData(seu,normalization.method = par_normalization.method,scale.factor =par_scale.factor)

## find variable features after filtering
seu<- FindVariableFeatures(seu, selection.method = par_selection.method, nfeatures = par_nfeatures)
topsel <- head(Seurat::VariableFeatures(seu), par_top)
write.csv(topsel, file = paste(output_dir,'/step3/info3/most_variable_genes_',sample_nameb[i],'.txt', sep=""), quote = TRUE, sep = ",")

## print variable features plot
vf_plot <- Seurat::VariableFeaturePlot(seu)
Seurat::LabelPoints(plot = vf_plot,points = topsel, repel = TRUE)
ggsave(paste(output_dir,'/step3/figs3/VariableFeaturePlot_',sample_nameb[i],".pdf",sep=""))

##Do not regress out cc genes
if (tolower(par_regress_cell_cycle_genes)=='no') {
seu<- ScaleData(seu, verbose = FALSE)
}

## Regress out cc genes
if (tolower(par_regress_cell_cycle_genes)=='yes') {
seu<- ScaleData(seu, vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)
}

## perform linear dimensional reduction
seu <- RunPCA(seu, npcs = par_npcs_pca, verbose = FALSE)
## print PCA 
DimPlot(seu, reduction = "pca", raster = FALSE)
ggsave(paste(output_dir,'/step3/figs3/',"dimplot_pca_",sample_nameb[i],".pdf",sep=""))
## print elbow plot
ElbowPlot(seu, ndims = par_npcs_pca)
ggsave(paste(output_dir,'/step3/figs3/',"elbowplot_",sample_nameb[i],".pdf",sep=""))

## save each individual Seurat object as RDS
saveRDS(seu, paste(output_dir,'/step3/objs3/',sample_nameb[i],".rds", sep=""))

## nefeature
n_feature <- Seurat::VlnPlot(seu, group.by= "orig.ident", features = "nFeature_RNA", pt.size = 0.001,raster = FALSE) + NoLegend()  + 
ggtitle("") + 
xlab("Control1") + 
ylab("Unique RNA transcripts") +
theme(axis.text.x = element_blank(),
axis.text.y = element_text(size = 12),
axis.title = element_text(size = 12, face = "bold")) +
scale_y_continuous(limits = c(0,15000), breaks = c(0,2500,5000,7500,10000,12500,15000)) +
geom_hline(yintercept = 1000, col = "dodgerblue2", linetype = "dashed", linewidth =1)
ggsave(paste(output_dir,'/manuscript_figures/','filt_feat_' ,sample_name[i],".pdf", sep=""), width = 3, height =4)


## ncount
n_count <- Seurat::VlnPlot(seu, group.by= "orig.ident", features = "nCount_RNA", pt.size = 0.01,raster = FALSE) + NoLegend() + 
ggtitle("") + 
xlab("Control1") + 
ylab("Total RNA transcripts") +
theme(axis.text.x = element_blank(),
axis.text.y = element_text(size = 12),
axis.title = element_text(size = 12, face = "bold")) +
scale_y_continuous(limits = c(0,250000),breaks = c(0,25000, 50000, 75000, 100000,125000, 150000, 175000, 200000, 225000, 250000)) +
geom_hline(yintercept = 1500, col = "dodgerblue2", linetype = "dashed", linewidth =1)
ggsave(paste(output_dir,'/manuscript_figures/','filt_count_' ,sample_name[i],".pdf", sep=""), width = 3, height =4)


## percent mt
mito <- Seurat::VlnPlot(seu, group.by= "orig.ident", features = "percent.mt", pt.size = 0.01,raster = FALSE) + NoLegend() + 
ggtitle("") + 
xlab("Control1") + 
ylab("Percent mitochondrial RNA") +
theme(axis.text.x = element_blank(),
axis.text.y = element_text(size = 12),
axis.title = element_text(size = 12, face = "bold")) +
scale_y_continuous(limits = c(0,100),breaks = c(0,5,10,11,100) ) +
geom_hline(yintercept = 10, col = "dodgerblue2", linetype = "dashed", linewidth =1) +
scale_y_break(c(11, 100),space = 0)
ggsave(paste(output_dir,'/manuscript_figures/','filt_mito_' ,sample_name[i],".pdf", sep=""), width = 3, height =4)

## percent ribo
ribo <- Seurat::VlnPlot(seu, group.by= "orig.ident", features = "percent.ribo", pt.size = 0.01,raster = FALSE) + NoLegend() + 
ggtitle("") + 
xlab("Control1") + 
ylab("Percent ribosomal RNA") +
theme(axis.text.x = element_blank(),
axis.text.y = element_text(size = 12),
axis.title = element_text(size = 12, face = "bold")) +
scale_y_continuous(limits = c(0,100),breaks = c(0,5,10,11,100) ) +
geom_hline(yintercept = 10, col = "dodgerblue2", linetype = "dashed", linewidth =1) +
scale_y_break(c(11, 100),space = 0)
ggsave(paste(output_dir,'/manuscript_figures/','filt_ribo_' ,sample_name[i],".pdf", sep=""), width = 3, height =4)

##########################################################################################
## Figure 3B-C
##########################################################################################
## load library
packages<-c('Seurat','ggplot2', 'dplyr', 'foreach', 'doParallel','DoubletFinder', 'Matrix', 'ggpubr', 'data.table')
lapply(packages, library, character.only = TRUE)

## load existing Seurat objects
source(paste(output_dir,'/job_info/parameters/step4_par.txt',sep=""))
sample_name<-list.files(path = paste(output_dir, "/step3/objs3",sep=""),pattern = "*.rds")
sample_nameb<-gsub(".rds","",sample_name)
if(length(sample_name)<1) {
print("You do not have any object from step 2 ")
}

for (i in 1:length(sample_name)) {
if (!grepl(".rds",tolower(sample_name[i]), fixed = TRUE)){
print(c(sample_name[i],"is not R rds"))
}
}   

## identify available cores for parallel processing
numCores <- detectCores()
cl <- makeCluster(numCores-1)
registerDoParallel(cl) 

## make a dataframe for expected doublet rates for each sample -- defined by the user in the parameters
doublet_rate_df <- data.frame(par_sample_names, par_expected_doublet_rate)
doublet_rate_df

## create empty list to be populated with existing Seurat objects
seu_list<-list()

##### if users want to remove doublets from downstream analyses

print('The following are deleted')
i_s = 1
seu<-readRDS(paste(output_dir,'/step3/objs3/',sample_name[i_s], sep=""))
print(sample_name[i_s])
sample_select <- sub(".rds.*", "", sample_name[i_s]) 
sample_select
seu[["Sample_ID"]] <- sample_select

## RunUMAP
seu <- RunUMAP(seu, dims = 1:par_RunUMAP_dims, n.neighbors =par_RunUMAP_n.neighbors)

## parameter sweep 
sweep.res.list_pbmc <- paramSweep_v3(seu, PCs = 1:par_PCs, sct = par_sct)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)

## select the pK that corresponds to max bcmvn to optimize doublet detection
pK <- bcmvn_pbmc %>% filter(BCmetric == max(BCmetric)) %>%select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate 
annotations <- seu@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           

## Seurat object-specific expected doublet rate
nExp_poi <- round(doublet_rate_df$par_expected_doublet_rate[doublet_rate_df$par_sample_names == sample_select]*nrow(seu@meta.data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## run doubletFinder 
seu <- doubletFinder_v3(seu, PCs = 1:par_PCs, pN = par_pN,pK = pK, nExp = nExp_poi.adj,reuse.pANN = FALSE, sct = par_sct)

## set droplet classification 
DF.classifications=colnames(seu[[]])[which(grepl('DF.classifications', colnames(seu[[]]), fixed=TRUE))]
DF.classifications_u=eval(parse(text = paste('unique(seu$',DF.classifications,')', sep='')))        

## print UMAP with doublet/singlet classification
xx <- DimPlot(seu, reduction = 'umap', group.by = DF.classifications, pt.size = 1) + 
ggtitle("") + 
xlab("UMAP1") + 
ylab("UMAP") +
theme(axis.text.x = element_text(size = 16),
axis.text.y = element_text(size = 16),
axis.title = element_text(size = 16, face = "bold"),
legend.text = element_text(size = 16),
legend.position = "none") +
scale_colour_manual(values = c("#F8766D", "lightgrey"))

ggsave(paste(output_dir,'/manuscript_figures/',"New_DF.classifications.pdf",sep=""))

##print doublet detection summary plot
#meta df
meta_df <- data.frame(seu@meta.data)
#string detect column
## pANN
df_pANN <- meta_df[,colnames(meta_df) %like% c("pANN")]
## classification
df_class <- meta_df[,colnames(meta_df) %like% c("DF.class")]
## bind
df_bind <- cbind(df_pANN, df_class)
df_bind <- data.frame(df_bind)
df_bind$df_pANN <- as.numeric(df_bind$df_pANN )

## plot pANN violin plot
pANN <- ggplot(df_bind, aes(x = df_class, y = df_pANN, fill = df_class)) + 
geom_violin() + 
geom_jitter(size=0.001) +
theme_classic() + 
xlab("Identity") + 
ylab("pANN") +
theme(axis.text.x = element_text(size = 16, colour = "black"),
axis.text.y = element_text(size = 16, colour = "black"),
axis.title = element_text(size = 16, face = "bold"),
legend.text = element_text(size = 16),
legend.title = element_blank(),
legend.position = "none") +
scale_fill_manual(values = c("#F8766D", "lightgrey"))

## plot number of doublets and singlets
count<- ggplot(df_bind, aes(x = df_class, fill = df_class)) + 
geom_bar() +
theme_classic() + 
xlab("Identity") + 
ylab("Number of droplets") +
geom_text(stat='count', aes(label=after_stat(count)), vjust=-1, size = 5) + 
theme(axis.text = element_text(size = 12)) +
theme(axis.text.x = element_text(size = 16, colour = "black"),
axis.text.y = element_text(size = 16, colour = "black"),
axis.title = element_text(size = 16, face = "bold"),
legend.text = element_text(size = 16),
legend.title = element_blank(), 
legend.position = "none") +
scale_fill_manual(values = c("#F8766D", "lightgrey"))

##print plot
ggarrange(pANN,count, ncol = 2)
ggsave(paste(output_dir,'/manuscript_figures/summary',sample_nameb[i_s],"_doublet_summary.pdf",sep=""))

##########################################################################################
## Figure 3D-E
##########################################################################################
### print UMAp -- integrated
library(RColorBrewer)
brewer.pal(11, 'Paired')

## load integrated UMAP
seu_int<-readRDS('/lustre03/project/6070393/COMMON/Dark_Genome/Saeid_test_temp3/smajic3/step5/objs5/seu_step5.rds')

##plot
DimPlot(seu_int, reduction = "umap", group.by="Sample_ID", raster = FALSE ) +
theme(axis.text.x = element_text(size = 16),
axis.text.y = element_text(size = 16),
axis.title = element_text(size = 16, face = "bold"),
legend.text = element_text(size = 16),
legend.title = element_blank()) +
xlab("UMAP1") + ylab("UMAP2") + 
scale_color_manual(values = rev(c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F","#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99"))) +
ggtitle("")
ggsave(paste(output_dir,'/manuscript_figures/integrateUMAP.pdf',sep=""))

## load merged UMAP
seu_int<-readRDS('/lustre03/project/6070393/COMMON/Dark_Genome/Saeid_test_temp3/smajic3/step5/objs5/seu_step5.rds')

##plot
DimPlot(seu_int, reduction = "umap", group.by="Sample_ID", raster = FALSE ) +
theme(axis.text.x = element_text(size = 16),
axis.text.y = element_text(size = 16),
axis.title = element_text(size = 16, face = "bold"),
legend.text = element_text(size = 16),
legend.title = element_blank()) +
xlab("UMAP1") + ylab("UMAP2") + 
scale_color_manual(values = rev(c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F","#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99"))) +
ggtitle("")
ggsave(paste(output_dir,'/manuscript_figures/mergeUMAP.pdf',sep=""))


