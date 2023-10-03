#!/usr/bin/env Rscript

##########################################
# v1.38
# step3: Quality control and filtering
##########################################

## load parameters
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]

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
foreach (i=1:length(sample_name)) %do% {    
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
    ggsave(paste(output_dir,'/step3/figs3/VariableFeaturePlot',sample_nameb[i],".png",sep=""))

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
            ggsave(paste(output_dir,'/step3/figs3/',"dimplot_pca",sample_nameb[i],".png",sep=""))
            ## print elbow plot
            ElbowPlot(seu, ndims = par_npcs_pca)
            ggsave(paste(output_dir,'/step3/figs3/',"elbowplot",sample_nameb[i],".png",sep=""))
    
    ## save each individual Seurat object as RDS
    saveRDS(seu, paste(output_dir,'/step3/objs3/',sample_nameb[i],".rds", sep=""))
    
    ## print QC violin plot
    Seurat::VlnPlot(seu, group.by= "orig.ident", features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo","S.Score", "G2M.Score"), pt.size = 0.001,ncol = 3, raster = FALSE) + NoLegend() #new code
    ggsave(paste(output_dir,'/step3/figs3/filtered_QC_vioplot_',sample_nameb[i],".png", sep=""))
    
    ## write meta infor available in the Seurat metdata
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
    cat("The number of features/genes and number of GEM/barcodes: \n")
    print(dim(seu))
    sink()

}

## save session information
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step3/info3/sessionInfo.txt', sep=""))
if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}