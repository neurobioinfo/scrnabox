#!/usr/bin/env Rscript

##########################################
# step3: Quality control and filtering
##########################################

## load parameters
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
NFRNAL=as.numeric(args[3])
NFRNAU=as.numeric(args[4])
NCRNAL=as.numeric(args[5])
NCRNAU=as.numeric(args[6])
PMTL=as.numeric(args[7])
PMTU=as.numeric(args[8])
GENEUMIL=as.numeric(args[9])
GENEUMIU=as.numeric(args[10])
#PRIBOL=as.numeric(args[11]) ######## NEEDS SAEIDS ATTENTION: Saeid, can you please modify the underlying code to make the two commented lines possible? This is to set the thresholds for percent ribosomal genes (it is the same aa mitochondria but for ribosome). For now I am setting them manually to NA
#PRIBOU=as.numeric(args[12])
PRIBOL <- 0
PRIBOU<- 10

## load library
.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr','Matrix', 'foreach', 'doParallel')
lapply(packages, library, character.only = TRUE)

## load parameters text file
source(paste(output_dir,'/job_info/parameters/step3_par.txt',sep=""))

## if user has existing Seurat object -- process Seurat objects and create list of distinct objects
if (exists("par_seurat_object")) {                                                  
sample_name<-list.files(path = par_seurat_object)
sample_nameb<-gsub(".rds","",sample_name)
if(length(sample_name)<1) {
   print("You do not have any existing Seurat object")
}

for (i in 1:length(sample_name)) {
  if (!grepl(".rds",tolower(sample_name[i]), fixed = TRUE)){
     print(c(sample_name[i],"is not R rds"))
  }
} 
} 

## if using Seurat object produced in Step 2 -- process Seurat objects and create list of distinct objects
if (!exists("par_seurat_object")) { 
sample_name<-list.files(path = paste(output_dir, "/step2/objs2",sep=""))
sample_nameb<-gsub(".rds","",sample_name)
if(length(sample_name)<1) {
   print("You do not have any object from step 2 ")
}

for (i in 1:length(sample_name)) {
  if (!grepl(".rds",tolower(sample_name[i]), fixed = TRUE)){
     print(c(sample_name[i],"is not R rds"))
  }
}   
}

## Set QC parameters if defined by the user
if (exists("par_nFeature_RNA_L")) NFRNAL=par_nFeature_RNA_L
if (exists("par_nFeature_RNA_U")) NFRNAU=par_nFeature_RNA_U

if (exists("par_nCount_RNA_L")) NCRNAL=par_nCount_RNA_L
if (exists("par_nCount_RNA_U")) NCRNAU=par_nCount_RNA_U

if (exists("par_mitochondria_percent_L")) PMTL=par_mitochondria_percent_L
if (exists("par_mitochondria_percent_U")) PMTU=par_mitochondria_percent_U

if (exists("par_log10GenesPerUMI_L")) GENEUMIL=par_log10GenesPerUMI_L
if (exists("par_log10GenesPerUMI_U")) GENEUMIU=par_log10GenesPerUMI_U

if (exists("par_ribosomal_percent_L")) PRIBOL=par_ribosomal_percent_L
if (exists("par_ribosomal_percent_U")) PRIBOU=par_ribosomal_percent_U

## print QC parameters
print(output_dir)
print(NFRNAL)
print(NFRNAU)
print(NCRNAL)
print(NCRNAU)
print(PMTL)
print(PMTU)
print(GENEUMIL)
print(GENEUMIU)
print(PRIBOL)
print(PRIBOU)

## set seed for replicability
set.seed(1234)

## detect number of cores
numCores <- detectCores()
cl <- makeCluster(numCores-1)
registerDoParallel(cl) 

###### QC and filtering if users have a preporcessed Seurat object
if (exists("par_seurat_object")) { 
foreach (i=1:length(sample_name)) %do% {    
    set.seed(1234)
    seu<-readRDS(paste(par_seurat_object, "/", sample_name[i], sep=""))

    ## calculate percentage of mitochondrial transcripts
    seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")

    ## calculate percentage of ribosomal transcripts
    seu[["percent.ribo"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^RP[SL]")    #NEW CODE
    
    ## if any QC parameters are not defined, take min/max of the dataset
    if (is.na(NFRNAL)) {
        NFRNAL=min(seu[["nFeature_RNA"]])
    }
    if (is.na(NFRNAU)) {
        NFRNAU=max(seu[["nFeature_RNA"]])
    }
    if (is.na(NCRNAL)) {
        NCRNAL=min(seu[["nCount_RNA"]])
    }
    if (is.na(NCRNAU)) {
        NCRNAU=max(seu[["nCount_RNA"]])
    }
    if (is.na(PMTL)) {
        PMTL=min(seu[["percent.mt"]])
    }
    if (is.na(PMTU)) {
        PMTU=max(seu[["percent.mt"]])
    }
    seu$log10GenesPerUMI <- log10(seu$nFeature_RNA) / log10(seu$nCount_RNA)
    if (is.na(GENEUMIL)) {
        GENEUMIL=min(seu[["log10GenesPerUMI"]])
    }
    if (is.na(GENEUMIU)) {
        GENEUMIU=max(seu[["log10GenesPerUMI"]])
    }
    if (is.na(PRIBOU)) {
        PRIBOU=max(seu[["percent.ribo"]]) #new code
    }
    if (is.na(PRIBOL)) {
        PRIBOL=min(seu[["percent.ribo"]])
    }

    ## subset dataset to only retain cells that pass the user-define QC metrics
    seu <- subset(seu, subset = nFeature_RNA > NFRNAL & nFeature_RNA < NFRNAU & nCount_RNA > NCRNAL & nCount_RNA < NCRNAU & PMTL < percent.mt & percent.mt < PMTU & log10GenesPerUMI > GENEUMIL & log10GenesPerUMI < GENEUMIU & percent.ribo < PRIBOU & PRIBOL < percent.ribo) 
    
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

    ## Normalize and scale individual Seurat object prior to cell-cycle scoring
     seu <- Seurat::NormalizeData(seu,normalization.method = par_normalization.method,scale.factor =par_scale.factor)
            ## Find variable features and print figure
            seu<- FindVariableFeatures(seu, selection.method = par_selection.method, nfeatures = par_nfeatures)
            topsel <- head(Seurat::VariableFeatures(seu), par_top)
            write.csv(topsel, file = paste(output_dir,'/step3/info3/most_variable_genes_',sample_nameb[i],'.txt', sep=""), quote = TRUE, sep = ",")
            vf_plot <- Seurat::VariableFeaturePlot(seu)
            Seurat::LabelPoints(plot = vf_plot,points = topsel, repel = TRUE)
            ggsave(paste(output_dir,'/step3/figs3/VariableFeaturePlot',sample_nameb[i],".png",sep=""))
            ## scale data
            seu<- ScaleData(seu, verbose = FALSE)
        
            ## perform linear dimensional reduction on individual Seurat objects
            seu <- RunPCA(seu, npcs = par_npcs_pca, verbose = FALSE)
            ## print PCA 
            DimPlot(seu, reduction = "pca")
            ggsave(paste(output_dir,'/step3/figs3/',"dimplot_pca",sample_nameb[i],".png",sep=""))
            ## print elbow plot
            ElbowPlot(seu)
            ggsave(paste(output_dir,'/step3/figs3/',"elbowplot",sample_nameb[i],".png",sep=""))
            ## print heatmap
            #Seurat::DimHeatmap(seu, dims = 1:par_dims, cells = par_cells, balanced = TRUE)            ## Saeid this does not produce a figure for some reason. I could not figure out whye. If we really want this figure then we can return and look into it.
            #ggsave(paste(output_dir,'/step3/figs3/',"dimheatplot.",sample_nameb[i],".png",sep=""))
            ## print UMAP
            seu <- RunUMAP(seu, dims = 1:par_dims_umap, n.neighbors =par_n.neighbors)
            Seurat::DimPlot(seu, reduction = "umap")
            ggsave(paste(output_dir,'/step3/figs3/',"dimplot_umap",sample_nameb[i],".png",sep=""))
        
    ## perform cell cycle scoring on individual Seurat objects
    seu <- CellCycleScoring(object = seu, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
    ## print violin plot for cell cycle score
    Seurat::VlnPlot(seu, features = c("S.Score", "G2M.Score"), group.by = "orig.ident",
    ncol = 4, pt.size = 0.1)
    ggsave(paste(output_dir,'/step3/figs3/cellcycle_',sample_nameb[i],".png", sep=""))

    ## save each individual Seurat object as RDS
    saveRDS(seu, paste(output_dir,'/step3/objs3/',sample_nameb[i],".rds", sep=""))
    Seurat::VlnPlot(seu, group.by= "orig.ident", features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo"), pt.size = 0.1,ncol = 4) + NoLegend() #new code
    ggsave(paste(output_dir,'/step3/figs3/QC_vioplot_',sample_nameb[i],".png", sep=""))
    write.csv(colnames(seu[[]]), file= paste(output_dir,'/step3/info3/meta_info_',sample_nameb[i],".txt", sep=""))
    
    ## save RNA expression matrix for each individual Seurat object
    if (tolower(par_save_RNA)=='yes') {
       mat <- GetAssayData(object = seu, assay = "RNA", slot = "data")
      writeMM(mat,file= paste(output_dir,'/step3/info3/',sample_nameb[i],"_RNA.txt", sep=""))
    }
    ## save metadata dataframe for each individual Seurat object
    if (tolower(par_save_metadata)=='yes') {
      write.csv(seu[[]], file = paste(output_dir,'/step3/info3/MetaData',sample_nameb[i],'.txt', sep=""), quote = TRUE, sep = ",")
    }
    ## save summary information for each individual Seurat object
    sink(paste(output_dir,'/step3/info3/summary_',sample_nameb[i],".txt", sep=""))
    cat("Summary of nCount_RNA: \n")
    print(summary(seu$nCount_RNA))
    cat("Summary of nFeature_RNA: \n")
    print(summary(seu$nFeature_RNA))
    cat("Summary of pt_mito: \n")
    print(summary(seu$percent.mt))
    cat("Summary of pt_ribo: \n") 
    print(summary(seu$percent.ribo)) 
    cat("The number of GEM/barcodes: \n")
    print(dim(seu))
    sink()

}

## save session information
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step3/info3/sessionInfo.txt', sep=""))
if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}
}


##### QC and filtering if using the Seurat objects produced in step 2
if (!exists("par_seurat_object")) {
foreach (i=1:length(sample_name)) %do% {    
    set.seed(1234)
    seu<-readRDS(paste(output_dir,'/step2/objs2/',sample_name[i], sep=""))

    ## calculate percentage of mitochondrial transcripts
    seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")

    ## calculate percentage of ribosomal transcripts
    seu[["percent.ribo"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^RP[SL]")    #NEW CODE
    
    ## if any QC parameters are not defined, take min/max of the dataset
    if (is.na(NFRNAL)) {
        NFRNAL=min(seu[["nFeature_RNA"]])
    }
    if (is.na(NFRNAU)) {
        NFRNAU=max(seu[["nFeature_RNA"]])
    }
    if (is.na(NCRNAL)) {
        NCRNAL=min(seu[["nCount_RNA"]])
    }
    if (is.na(NCRNAU)) {
        NCRNAU=max(seu[["nCount_RNA"]])
    }
    if (is.na(PMTL)) {
        PMTL=min(seu[["percent.mt"]])
    }
    if (is.na(PMTU)) {
        PMTU=max(seu[["percent.mt"]])
    }
    seu$log10GenesPerUMI <- log10(seu$nFeature_RNA) / log10(seu$nCount_RNA)
    if (is.na(GENEUMIL)) {
        GENEUMIL=min(seu[["log10GenesPerUMI"]])
    }
    if (is.na(GENEUMIU)) {
        GENEUMIU=max(seu[["log10GenesPerUMI"]])
    }
    if (is.na(PRIBOU)) {
        PRIBOU=max(seu[["percent.ribo"]]) #new code
    }
    if (is.na(PRIBOL)) {
        PRIBOL=min(seu[["percent.ribo"]])
    }

    ## subset dataset to only retain cells that pass the user-define QC metrics
    seu <- subset(seu, subset = nFeature_RNA > NFRNAL & nFeature_RNA < NFRNAU & nCount_RNA > NCRNAL & nCount_RNA < NCRNAU & PMTL < percent.mt & percent.mt < PMTU & log10GenesPerUMI > GENEUMIL & log10GenesPerUMI < GENEUMIU & percent.ribo < PRIBOU & PRIBOL < percent.ribo) 
    
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

    ## Normalize and scale individual Seurat object prior to cell-cycle scoring
     seu <- Seurat::NormalizeData(seu,normalization.method = par_normalization.method,scale.factor =par_scale.factor)
            ## Find variable features and print figure
            seu<- FindVariableFeatures(seu, selection.method = par_selection.method, nfeatures = par_nfeatures)
            topsel <- head(Seurat::VariableFeatures(seu), par_top)
            write.csv(topsel, file = paste(output_dir,'/step3/info3/most_variable_genes_',sample_nameb[i],'.txt', sep=""), quote = TRUE, sep = ",")
            vf_plot <- Seurat::VariableFeaturePlot(seu)
            Seurat::LabelPoints(plot = vf_plot,points = topsel, repel = TRUE)
            ggsave(paste(output_dir,'/step3/figs3/VariableFeaturePlot',sample_nameb[i],".png",sep=""))
            ## scale data
            seu<- ScaleData(seu, verbose = FALSE)
        
            ## perform linear dimensional reduction on individual Seurat objects
            seu <- RunPCA(seu, npcs = par_npcs_pca, verbose = FALSE)
            ## print PCA 
            DimPlot(seu, reduction = "pca")
            ggsave(paste(output_dir,'/step3/figs3/',"dimplot_pca",sample_nameb[i],".png",sep=""))
            ## print elbow plot
            ElbowPlot(seu)
            ggsave(paste(output_dir,'/step3/figs3/',"elbowplot",sample_nameb[i],".png",sep=""))
            ## print heatmap
            #Seurat::DimHeatmap(seu, dims = 1:par_dims, cells = par_cells, balanced = TRUE)            ## Saeid this does not produce a figure for some reason. I could not figure out whye. If we really want this figure then we can return and look into it.
            #ggsave(paste(output_dir,'/step3/figs3/',"dimheatplot.",sample_nameb[i],".png",sep=""))
            ## print UMAP
            seu <- RunUMAP(seu, dims = 1:par_dims_umap, n.neighbors =par_n.neighbors)
            Seurat::DimPlot(seu, reduction = "umap")
            ggsave(paste(output_dir,'/step3/figs3/',"dimplot_umap",sample_nameb[i],".png",sep=""))
        
    ## perform cell cycle scoring on individual Seurat objects
    seu <- CellCycleScoring(object = seu, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
    ## print violin plot for cell cycle score
    Seurat::VlnPlot(seu, features = c("S.Score", "G2M.Score"), group.by = "orig.ident",
    ncol = 4, pt.size = 0.1)
    ggsave(paste(output_dir,'/step3/figs3/cellcycle_',sample_nameb[i],".png", sep=""))

    ## save each individual Seurat object as RDS
    saveRDS(seu, paste(output_dir,'/step3/objs3/',sample_nameb[i],".rds", sep=""))
    Seurat::VlnPlot(seu, group.by= "orig.ident", features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo"), pt.size = 0.1,ncol = 4) + NoLegend() #new code
    ggsave(paste(output_dir,'/step3/figs3/QC_vioplot_',sample_nameb[i],".png", sep=""))
    write.csv(colnames(seu[[]]), file= paste(output_dir,'/step3/info3/meta_info_',sample_nameb[i],".txt", sep=""))
    
    ## save RNA expression matrix for each individual Seurat object
    if (tolower(par_save_RNA)=='yes') {
       mat <- GetAssayData(object = seu, assay = "RNA", slot = "data")
      writeMM(mat,file= paste(output_dir,'/step3/info3/',sample_nameb[i],"_RNA.txt", sep=""))
    }
    ## save metadata dataframe for each individual Seurat object
    if (tolower(par_save_metadata)=='yes') {
      write.csv(seu[[]], file = paste(output_dir,'/step3/info3/MetaData',sample_nameb[i],'.txt', sep=""), quote = TRUE, sep = ",")
    }
    ## save summary information for each individual Seurat object
    sink(paste(output_dir,'/step3/info3/summary_',sample_nameb[i],".txt", sep=""))
    cat("Summary of nCount_RNA: \n")
    print(summary(seu$nCount_RNA))
    cat("Summary of nFeature_RNA: \n")
    print(summary(seu$nFeature_RNA))
    cat("Summary of pt_mito: \n")
    print(summary(seu$percent.mt))
    cat("Summary of pt_ribo: \n") 
    print(summary(seu$percent.ribo)) 
    cat("The number of GEM/barcodes: \n")
    print(dim(seu))
    sink()

}
}

## save session information
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step3/info3/sessionInfo.txt', sep=""))
if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}