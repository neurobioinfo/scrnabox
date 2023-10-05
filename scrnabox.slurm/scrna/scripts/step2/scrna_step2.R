#!/usr/bin/env Rscript

####################
# v1.38
# step2 -- ambient rna removal and create Seurat object
####################

## load parameters
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]

## load library
.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr', 'foreach', 'doParallel','Matrix', 'ggpubr')
lapply(packages, library, character.only = TRUE)

## load parameters
source(paste(output_dir,'/job_info/parameters/step2_par.txt',sep=""))

## detect number of available cores
numCores <- detectCores()
cl <- makeCluster(numCores-1)
registerDoParallel(cl) 

########################################################################################################
## Remove ambient RNA
########################################################################################################
if (tolower(par_ambient_RNA)=="yes") {
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

      ## save ambient RNA-corrected feature barcode matrices
      DropletUtils:::write10xCounts(paste0(output_dir, '/step2_ambient/',sample_name[i]), out)
      saveRDS(sc, paste(output_dir,'/step2/info2/',sample_name[i],'_ambient_rna_summary.rds', sep=''))

      ##print txt file file with estimated ambient RNA
      est_amb <- data.frame(sc$fit$rhoEst)
      write.csv(est_amb, file = paste(output_dir,'/step2/info2/estimated_ambient_RNA_',sample_name[i],'.txt', sep=""), quote = TRUE, sep = ",")

      ## create Seurat object with feature-barcode matrices correct for ambient RNA expression and filter according to user-defined parameters
      datadirs <- file.path(paste0(output_dir, '/step2_ambient/',sample_name[i]))
      sparse_matrix <- Seurat::Read10X(data.dir = datadirs)
      seurat_object <- Seurat::CreateSeuratObject(counts = sparse_matrix, min.cells=par_min.cells_L,min.features= 1) # only keep cells expressing atleast one gene
      
      ## name Seurat object
      nam <- paste("seurat_object", sample_name[i], sep = ".")
      assign(nam, seurat_object)
      seu<-get(nam)
      
      ## calculate percent mitochondrial
      seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")
      seu <- subset(seu, subset = percent.mt < 100)
      print(i)

      ## calculate percent ribosomal 
      seu[["percent.ribo"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^RP[SL]") 
      seu <- subset(seu, subset = percent.ribo < 100)
      print(i)

      ## Normalize and scale individual Seurat object prior to cell-cycle scoring
      seu <- Seurat::NormalizeData(seu,normalization.method = par_normalization.method,scale.factor =par_scale.factor)
            ## Find variable features and print figure
            seu<- FindVariableFeatures(seu, selection.method = par_selection.method, nfeatures = par_nfeatures)
            ## scale data
            seu<- ScaleData(seu, features = c(VariableFeatures(seu), cc.genes$g2m.genes, cc.genes$s.genes),verbose = FALSE)  
            ## perform linear dimensional reduction on individual Seurat objects
            seu <- RunPCA(seu, verbose =FALSE)
            
      ## perform cell cycle scoring on individual Seurat objects
      seu <- CellCycleScoring(object = seu, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
      
      #save seurat object       
      saveRDS(seu,paste(output_dir,'/step2/objs2/',sample_name[i],'.rds', sep=''),compress=TRUE)

      ## PCA cell cycle score
      seu <- RunPCA(seu, features = c(cc.genes$g2m.genes, cc.genes$s.genes),raster = FALSE)
      DimPlot(seu, group.by = "Phase", pt.size =2, raster =FALSE)
      ggsave(paste(output_dir,'/step2/figs2/cell_cyle_dim_plot_',sample_name[i],".pdf", sep=""))
      
      ## print violin plot for QC metrics
      Seurat::VlnPlot(seu, group.by= "orig.ident", features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo","S.Score", "G2M.Score"), pt.size = 0.01,ncol = 3,raster = FALSE) + NoLegend()
      ggsave(paste(output_dir,'/step2/figs2/vioplot_',sample_name[i],".pdf", sep=""))

      ## zoomed in violinplot
      ## nefeature
      n_feature <- Seurat::VlnPlot(seu, group.by= "orig.ident", features = "nFeature_RNA", pt.size = 0.001,ncol = 3,raster = FALSE) + NoLegend() + scale_y_continuous(limits = function(x){
      c(min(x), (mean(x)))
      })
      ## ncount
      n_count <- Seurat::VlnPlot(seu, group.by= "orig.ident", features = "nCount_RNA", pt.size = 0.01,ncol = 3,raster = FALSE) + NoLegend() + scale_y_continuous(limits = function(x){
            c(min(x), (mean(x)))
            })
      ## percent mt
      mito <- Seurat::VlnPlot(seu, group.by= "orig.ident", features = "percent.mt", pt.size = 0.01,ncol = 3,raster = FALSE) + NoLegend() + scale_y_continuous(limits = function(x){
            c(min(x), (mean(x)))
            })
            
      ## percent ribo
      ribo <- Seurat::VlnPlot(seu, group.by= "orig.ident", features = "percent.ribo", pt.size = 0.01,ncol = 3,raster = FALSE) + NoLegend() + scale_y_continuous(limits = function(x){
            c(min(x), (mean(x)))
            })

      ## S.Score
      s_score <- Seurat::VlnPlot(seu, group.by= "orig.ident", features = "S.Score", pt.size = 0.01,ncol = 3,raster = FALSE) + NoLegend() + scale_y_continuous(limits = function(x){
            c(min(x), (mean(x)))
            })
      
      ## G2M.Score
      G2M_score <- Seurat::VlnPlot(seu, group.by= "orig.ident", features = "G2M.Score", pt.size = 0.01,ncol = 3,raster = FALSE) + NoLegend() + scale_y_continuous(limits = function(x){
            c(min(x), (mean(x)))
            })
            
      ## merge all plots 
      ggarrange(n_feature, n_count,mito,ribo, s_score, G2M_score,  ncol = 3, nrow = 2)
      ggsave(paste(output_dir,'/step2/figs2/zoomed_in_vioplot_',sample_name[i],".pdf", sep=""))

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
        cat("The number of features/genes and number of GEM/barcodes: \n")
        print(dim(seu))
        sink()
        write.csv(colnames(seu[[]]), file= paste(output_dir,'/step2/info2/meta_info_',sample_name[i],".txt", sep=""))
      
      ## save RNA exprression matrix
      if (tolower(par_save_RNA)=='yes') {
          mat <- GetAssayData(object = seu, assay = "RNA", slot = "data")
          writeMM(mat,file= paste(output_dir,'/step2/info2/',sample_name[i],"_RNA.txt", sep=""))
          }

      ## save metadata dataframe
      if (tolower(par_save_metadata)=='yes') {
          write.csv(seu[[]], file = paste(output_dir,'/step2/info2/MetaData_',sample_name[i],'.txt', sep=""), quote = TRUE, sep = ",")
          }
      }
}

########################################################################################################
## Do not remove  ambient RNA
########################################################################################################
if (tolower(par_ambient_RNA)=="no") {  

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

      ## save ambient RNA-corrected feature barcode matrices
      saveRDS(sc, paste(output_dir,'/step2/info2/',sample_name[i],'_ambient_rna_summary.rds', sep=''))

      ##print txt file file with estimated ambient RNA
      est_amb <- data.frame(sc$fit$rhoEst)
      write.csv(est_amb, file = paste(output_dir,'/step2/info2/estimated_ambient_RNA_',sample_name[i],'.txt', sep=""), quote = TRUE, sep = ",")

  }

#########  #########  #########  #########  #########  #########  #########  #########  #########  #########  #########  #########  #########  #########  #########  #########  #########  
  list<-dir(path = paste(output_dir, "/step1",sep=""),full.names = TRUE)
  sample_name<-dir(path = paste(output_dir, "/step1",sep=""))
  foreach (i=1:length(sample_name)) %do% {    
      ## create Seurat object for each sample
      datadirs <- file.path(list[i],   "ouput_folder","outs","filtered_feature_bc_matrix")
      names(datadirs)=sample_name[i]
      sparse_matrix <- Seurat::Read10X(data.dir = datadirs)
      seurat_object <- Seurat::CreateSeuratObject(counts = sparse_matrix, min.cells=par_min.cells_L,min.features= 1) # only keep cells expressing atleast one gene
      ## name seurat object
      nam <- paste("seurat_object", sample_name[i], sep = ".")
      assign(nam, seurat_object)
      seu<-get(nam)

      ## calculate percent mitochondrial
      seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")
      seu <- subset(seu, subset = percent.mt < 100)
      print(i)
      
      ## calculate percent ribosomal 
      seu[["percent.ribo"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^RP[SL]") 
      seu <- subset(seu, subset = percent.ribo < 100)
      print(i)

            ## Normalize and scale individual Seurat object prior to cell-cycle scoring
            seu <- Seurat::NormalizeData(seu,normalization.method = par_normalization.method,scale.factor =par_scale.factor)
            ## Find variable features and print figure
            seu<- FindVariableFeatures(seu, selection.method = par_selection.method, nfeatures = par_nfeatures)
            ## scale data
            #seu<- ScaleData(seu, features = rownames(seu),verbose = FALSE)
            seu<- ScaleData(seu, features = c(VariableFeatures(seu), cc.genes$g2m.genes, cc.genes$s.genes),verbose = FALSE)   #### we have to see if this works?
            ## perform linear dimensional reduction on individual Seurat objects
            seu <- RunPCA(seu, verbose =FALSE)
            
      ## perform cell cycle scoring on individual Seurat objects
      seu <- CellCycleScoring(object = seu, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
            
      #save seurat object       
      saveRDS(seu,paste(output_dir,'/step2/objs2/',sample_name[i],'.rds', sep=''),compress=TRUE)

      ## visualize cell cycle score PCA
      seu <- RunPCA(seu, features = c(cc.genes$g2m.genes, cc.genes$s.genes),raster = FALSE)
      DimPlot(seu, group.by = "Phase", pt.size =2, raster =FALSE)
      ggsave(paste(output_dir,'/step2/figs2/cell_cyle_dim_plot_',sample_name[i],".pdf", sep=""))

      ## print violin plot for QC metrics
      Seurat::VlnPlot(seu, group.by= "orig.ident", features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo","S.Score", "G2M.Score"), pt.size = 0.01,ncol = 3,raster = FALSE) + NoLegend()
      ggsave(paste(output_dir,'/step2/figs2/vioplot_',sample_name[i],".pdf", sep=""))

      ## zoomed in violinplot
      ## nefeature
      n_feature <- Seurat::VlnPlot(seu, group.by= "orig.ident", features = "nFeature_RNA", pt.size = 0.001,ncol = 3,raster = FALSE) + NoLegend() + scale_y_continuous(limits = function(x){
      c(min(x), (mean(x)))
      })
      ## ncount
      n_count <- Seurat::VlnPlot(seu, group.by= "orig.ident", features = "nCount_RNA", pt.size = 0.01,ncol = 3,raster = FALSE) + NoLegend() + scale_y_continuous(limits = function(x){
            c(min(x), (mean(x)))
            })
      ## percent mt
      mito <- Seurat::VlnPlot(seu, group.by= "orig.ident", features = "percent.mt", pt.size = 0.01,ncol = 3,raster = FALSE) + NoLegend() + scale_y_continuous(limits = function(x){
            c(min(x), (mean(x)))
            })
            
      ## percent ribo
      ribo <- Seurat::VlnPlot(seu, group.by= "orig.ident", features = "percent.ribo", pt.size = 0.01,ncol = 3,raster = FALSE) + NoLegend() + scale_y_continuous(limits = function(x){
            c(min(x), (mean(x)))
            })

      ## S.Score
      s_score <- Seurat::VlnPlot(seu, group.by= "orig.ident", features = "S.Score", pt.size = 0.01,ncol = 3,raster = FALSE) + NoLegend() + scale_y_continuous(limits = function(x){
            c(min(x), (mean(x)))
            })
      
      ## G2M.Score
      G2M_score <- Seurat::VlnPlot(seu, group.by= "orig.ident", features = "G2M.Score", pt.size = 0.01,ncol = 3,raster = FALSE) + NoLegend() + scale_y_continuous(limits = function(x){
            c(min(x), (mean(x)))
            })
            
      ## merge all plots 
      ggarrange(n_feature, n_count,mito,ribo, s_score, G2M_score,  ncol = 3, nrow = 2)
      ggsave(paste(output_dir,'/step2/figs2/zoomed_in_vioplot_',sample_name[i],".pdf", sep=""))
            
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
        cat("The number of features/genes and number of GEM/barcodes: \n")
        print(dim(seu))
        sink()
        write.csv(colnames(seu[[]]), file= paste(output_dir,'/step2/info2/meta_info_',sample_name[i],".txt", sep=""))

      ## save RNA expression matrix
      if (tolower(par_save_RNA)=='yes') {
          mat <- GetAssayData(object = seu, assay = "RNA", slot = "data")
        writeMM(mat,file= paste(output_dir,'/step2/info2/',sample_name[i],"_RNA.txt", sep=""))
        }
      
      ## save metadata dataframe
      if (tolower(par_save_metadata)=='yes') {
        write.csv(seu[[]], file = paste(output_dir,'/step2/info2/MetaData_',sample_name[i],'.txt', sep=""), quote = TRUE, sep = ",")
        }
      }
}

## write session information
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step2/info2/sessionInfo.txt', sep=""))

if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}
