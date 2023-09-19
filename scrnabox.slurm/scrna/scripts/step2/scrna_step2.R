#!/usr/bin/env Rscript

####################
# step2 -- ambient rna removal and create Seurat object
####################

## load parameters
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]

## load library
.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr', 'foreach', 'doParallel','Matrix')
lapply(packages, library, character.only = TRUE)

## load parameters
source(paste(output_dir,'/job_info/parameters/step2_par.txt',sep=""))

## detect number of available cores
numCores <- detectCores()
cl <- makeCluster(numCores-1)
registerDoParallel(cl) 

###### use exisisting feature-barcode matrices without running step 1 (Cell ranger)
if (exists("par_count_matrices")) {
    sample_name <- list.dirs(par_count_matrices,full.names = FALSE,recursive = FALSE) 
    foreach (i=1:length(sample_name)) %do% {    
        ## load barcodes
        pre <- list.files(paste0(par_count_matrices,"/",sample_name[i]),"barcodes.tsv.gz", recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
        if (rlang::is_empty(pre)){
          message("ERROR: barcodes.tsv.gz does dont exit")
        }
        ## load features
        pre <- list.files(paste0(par_count_matrices,"/",sample_name[i]),"features.tsv.gz", recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
        if (rlang::is_empty(pre)){
          message("ERROR: features.tsv.gz does dont exit")
        }
        ## load counts
        pre <- list.files(paste0(par_count_matrices,"/",sample_name[i]),"matrix.mtx.gz", recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
        if (rlang::is_empty(pre)){
          message("ERROR: matrix.mtx.gz does dont exit")
        }
        datadirs=dirname(pre)
        names(datadirs)=sample_name[i]
        
        ## create Seurat object and filter according to user-defined parameters
        sparse_matrix <- Seurat::Read10X(data.dir = datadirs)
        seurat_object <- Seurat::CreateSeuratObject(counts = sparse_matrix, min.cells=par_min.cells_L, min.features=par_min.features_L)
        nam <- paste("seurat_object", sample_name[i], sep = ".")
        assign(nam, seurat_object)
        saveRDS(get(nam),paste(output_dir,'/step2/objs2/',sample_name[i],".rds", sep=""),compress=TRUE)
        seu<-get(nam)
        
        ## calculate percent mitochondria
        seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")
        seu <- subset(seu, subset = percent.mt < 100)
        print(i)
        ## calculate percent ribosomal 
        seu[["percent.ribo"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^RP[SL]") 
        seu <- subset(seu, subset = percent.ribo < 100)
        print(i)
        ## print violin plot for QC metrics
        Seurat::VlnPlot(seu, group.by= "orig.ident", features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo"), pt.size = 0.1,ncol = 3) + NoLegend()
        ggsave(paste(output_dir,'/step2/figs2/vioplot_',sample_name[i],".png", sep=""))
        ## print summary information
        sink(paste(output_dir,'/step2/info2/summary_',sample_name[i],".txt", sep=""))
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
        write.csv(colnames(seu[[]]), file= paste(output_dir,'/step2/info2/meta_info',sample_name[i],".txt", sep=""))
        
        ## save RNA expression matrix
        if (tolower(par_save_RNA)=='yes') {
            mat <- GetAssayData(object = seu, assay = "RNA", slot = "data")
        writeMM(mat,file= paste(output_dir,'/step2/info2/',sample_name[i],"_RNA.txt", sep=""))
        }

        ## save metadata dataframe
        if (tolower(par_save_metadata)=='yes') {
        write.csv(seu[[]], file = paste(output_dir,'/step2/info2/seu_MetaData',sample_name[i],'.txt', sep=""), quote = TRUE, sep = ",")
        }
        }
}

###### use feature-barcode matrices produced from step 1 (cell ranger) and remove ambient RNA
if (tolower(par_ambient_RNA)=="yes" & !exists("par_count_matrices")) {
  list<-dir(path = paste(output_dir, "/step1",sep=""),full.names = TRUE)
  sample_name<-dir(path = paste(output_dir, "/step1",sep=""))
  library(SoupX)
  library(MatrixGenerics)
  library(BiocGenerics)
  library(S4Vectors)
  library(IRanges)
  library(GenomeInfoDb)
  library(GenomicRanges)
  library(Biobase)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(DropletUtils)  
  
  ## ambient RNA removal using SoupX
  foreach (i=1:length(sample_name)) %do% { 
      sc = load10X(paste0(list[i],"/ouput_folder","/outs"))
      sc = autoEstCont(sc) 
      out = adjustCounts(sc)
      dir.create(file.path(output_dir, 'step2_ambient'), showWarnings = FALSE)
      dir0 <- paste0(output_dir, '/step2_ambient/',sample_name[i])
      if (file.exists(dir0)) {
        unlink(dir0,recursive = TRUE)
      }

      ## save ambient RNA-corrected feature barcode atrices
      DropletUtils:::write10xCounts(paste0(output_dir, '/step2_ambient/',sample_name[i]), out)
      saveRDS(sc, paste(output_dir,'/step2/info2/',sample_name[i],'_ambient_rna_summary.rds', sep=''))
      
      ## create Seurat object with feature-barcode matrices correct for ambient RNA expression and filter according to user-defined parameters
      datadirs <- file.path(paste0(output_dir, '/step2_ambient/',sample_name[i]))
      sparse_matrix <- Seurat::Read10X(data.dir = datadirs)
      seurat_object <- Seurat::CreateSeuratObject(counts = sparse_matrix, min.cells=par_min.cells_L, min.features=par_min.features_L)
      nam <- paste("seurat_object", sample_name[i], sep = ".")
      assign(nam, seurat_object)
      saveRDS(get(nam),paste(output_dir,'/step2/objs2/',sample_name[i],".rds", sep=""),compress=TRUE)
      seu<-get(nam)
      
      ## calculate percent mitochondrial
      seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")
      seu <- subset(seu, subset = percent.mt < 100)
      print(i)

      ## calculate percent ribosomal 
      seu[["percent.ribo"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^RP[SL]") 
      seu <- subset(seu, subset = percent.ribo < 100)
      print(i)
      
      ## print violin plot for QC metrics
      Seurat::VlnPlot(seu, group.by= "orig.ident", features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo"), pt.size = 0.1,ncol = 3) + NoLegend()
      ggsave(paste(output_dir,'/step2/figs2/vioplot_',sample_name[i],".png", sep=""))

       ## print summary information
        sink(paste(output_dir,'/step2/info2/summary_',sample_name[i],".txt", sep=""))
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
        write.csv(colnames(seu[[]]), file= paste(output_dir,'/step2/info2/meta_info',sample_name[i],".txt", sep=""))
      
      ## save RNA exprression matrix
      if (tolower(par_save_RNA)=='yes') {
          mat <- GetAssayData(object = seu, assay = "RNA", slot = "data")
          writeMM(mat,file= paste(output_dir,'/step2/info2/',sample_name[i],"_RNA.txt", sep=""))
          }

      ## save metadata dataframe
      if (tolower(par_save_metadata)=='yes') {
          write.csv(seu[[]], file = paste(output_dir,'/step2/info2/MetaData',sample_name[i],'.txt', sep=""), quote = TRUE, sep = ",")
          }
      }
}
  
###### use feature-barcode matrices produced from step 1 (cell ranger) and do not remove  ambient RNA
if (tolower(par_ambient_RNA)=="no" & !exists("par_count_matrices")) {   
  list<-dir(path = paste(output_dir, "/step1",sep=""),full.names = TRUE)
  sample_name<-dir(path = paste(output_dir, "/step1",sep=""))
  foreach (i=1:length(sample_name)) %do% {    
      ## create Seurat object for each sample
      datadirs <- file.path(list[i],   "ouput_folder","outs","raw_feature_bc_matrix")
      names(datadirs)=sample_name[i]
      sparse_matrix <- Seurat::Read10X(data.dir = datadirs)
      seurat_object <- Seurat::CreateSeuratObject(counts = sparse_matrix, min.cells=par_min.cells_L, min.features=par_min.features_L)
      ## save Seurat object for each sample
      nam <- paste("seurat_object", sample_name[i], sep = ".")
      assign(nam, seurat_object)
      saveRDS(get(nam),paste(output_dir,'/step2/objs2/',sample_name[i],'.rds', sep=''),compress=TRUE)
      seu<-get(nam)

      ## calculate percent mitochondrial
      seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")
      seu <- subset(seu, subset = percent.mt < 100)
      print(i)
      
      ## calculate percent ribosomal 
      seu[["percent.ribo"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^RP[SL]") 
      seu <- subset(seu, subset = percent.ribo < 100)
      print(i)
      
      ## print violin plot for QC metrics
      Seurat::VlnPlot(seu, group.by= "orig.ident", features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo"), pt.size = 0.1,ncol = 3) + NoLegend()
      ggsave(paste(output_dir,'/step2/figs2/vioplot_',sample_name[i],".png", sep=""))
      
       ## print summary information
        sink(paste(output_dir,'/step2/info2/summary_',sample_name[i],".txt", sep=""))
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
        write.csv(colnames(seu[[]]), file= paste(output_dir,'/step2/info2/meta_info',sample_name[i],".txt", sep=""))

      ## save RNA expression matrix
      if (tolower(par_save_RNA)=='yes') {
          mat <- GetAssayData(object = seu, assay = "RNA", slot = "data")
        writeMM(mat,file= paste(output_dir,'/step2/info2/',sample_name[i],"_RNA.txt", sep=""))
        }
      
      ## save metadata dataframe
      if (tolower(par_save_metadata)=='yes') {
        write.csv(seu[[]], file = paste(output_dir,'/step2/info2/MetaData',sample_name[i],'.txt', sep=""), quote = TRUE, sep = ",")
        }
      }
}

## write session information
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step2/info2/sessionInfo.txt', sep=""))
if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}