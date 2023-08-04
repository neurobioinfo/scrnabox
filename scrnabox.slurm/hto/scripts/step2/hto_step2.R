#!/usr/bin/env Rscript

####################
# step2 
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]


.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr', 'foreach', 'doParallel','Matrix')
lapply(packages, library, character.only = TRUE)
source(paste(output_dir,'/job_info/parameters/step2_par.txt',sep=""))

list<-dir(path = paste(output_dir, "/step1",sep=""),full.names = TRUE)
sample_name<-dir(path = paste(output_dir, "/step1",sep=""))


numCores <- detectCores()
cl <- makeCluster(numCores-1)
registerDoParallel(cl) 
ambient_RNA="NO" #### TESTAGAIN
if (exists("count_matrices")) {
   ##########################
       sample_name <- list.dirs(count_matrices,full.names = FALSE,recursive = FALSE) 
       foreach (i=1:length(sample_name)) %do% {    
            pre <- list.files(paste0(count_matrices,"/",sample_name[i]),"barcodes.tsv.gz", recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
            if (rlang::is_empty(pre)){
                message("ERROR: barcodes.tsv.gz does dont exit")
            }
            pre <- list.files(paste0(count_matrices,"/",sample_name[i]),"features.tsv.gz", recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
            if (rlang::is_empty(pre)){
                message("ERROR: features.tsv.gz does dont exit")
            }
            pre <- list.files(paste0(count_matrices,"/",sample_name[i]),"matrix.mtx.gz", recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
            if (rlang::is_empty(pre)){
                message("ERROR: matrix.mtx.gz does dont exit")
            }
            datadirs=dirname(pre)
        #   datadirs <- file.path(list[i],   "ouput_folder","outs","raw_feature_bc_matrix")
            names(datadirs)=sample_name[i]
            sparse_matrix <- Seurat::Read10X(data.dir = datadirs)
            seurat_object <- Seurat::CreateSeuratObject(counts = sparse_matrix$`Gene Expression`, min.cells=min.cells_L, min.features=min.features_L)
            seurat_object[['HTO']] = Seurat::CreateAssayObject(counts = sparse_matrix$`Antibody Capture`, min.cells=min.cells_L, min.features=min.features_L)
            # seurat_object <- Seurat::CreateSeuratObject(counts = sparse_matrix$`Gene Expression`)
            # seurat_object[['HTO']] = Seurat::CreateAssayObject(counts = sparse_matrix$`Antibody Capture`)
            nam <- paste("seurat_object", sample_name[i], sep = ".")
            assign(nam, seurat_object)
            saveRDS(get(nam),paste(output_dir,'/step2/objs2',"/seu",i,".rds", sep=""),compress=TRUE)
            seu<-get(nam)
            seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")
            seu <- subset(seu, subset = percent.mt < 100)
            print(i)
            # png(file = paste(output_dir,'/step2/figs/vioplot',sample_name$V1[i],".png", sep=""))
            Seurat::VlnPlot(seu, group.by= "orig.ident", features = c("nFeature_RNA","nCount_RNA","percent.mt"), pt.size = 0.1,ncol = 3) + NoLegend()
            # dev.off()
            ggsave(paste(output_dir,'/step2/figs2/vioplot_seu',i,".png", sep=""))
            write.csv(colnames(seu[[]]), file= paste(output_dir,'/step2/info2/meta_info_seu',i,".txt", sep=""))
            #
            sink(paste(output_dir,'/step2/info2/summary_seu',i,".txt", sep=""))
            cat("Summary of nCount_RNA: \n")
            print(summary(seu$nCount_RNA))
            cat("Summary of nFeature_RNA: \n")
            print(summary(seu$nFeature_RNA))
            cat("Summary of pt_mito: \n")
            print(summary(seu$percent.mt))
            cat("The number of GEM/barcodes: \n")
            print(dim(seu))
            sink()
            if (tolower(Save_RNA)=='yes') {
                mat <- GetAssayData(object = seu, assay = "RNA", slot = "data")
            #  write.csv(mat, paste(output_dir,'/step2/info2/',sample_name[i],"_RNA.csv", sep=""))
                writeMM(mat,file= paste(output_dir,'/step2/info2/seu',i,"_RNA.txt", sep=""))
            }
            if (tolower(Save_metadata)=='yes') {
                write.csv(seu[[]], file = paste(output_dir,'/step2/info2/seu_MetaData',i,'.txt', sep=""), quote = TRUE, sep = ",")
            }
        }
  # writeLines(capture.output(sessionInfo()), paste(output_dir,'/step2/info2/sessionInfo.txt', sep=""))
}

if (tolower(ambient_RNA)=="yes" & !exists("count_matrices")) {
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
    foreach (i=1:length(sample_name)) %do% { 
        sc = load10X(paste0(list[i],"/ouput_folder","/outs"))
        sc = autoEstCont(sc) 
        out = adjustCounts(sc)
        dir.create(file.path(output_dir, 'step1_ambient'), showWarnings = FALSE)
        dir0 <- paste0(output_dir, '/step1_ambient/',sample_name[i])
        if (file.exists(dir0)) {
          unlink(dir0,recursive = TRUE)
        }
        DropletUtils:::write10xCounts(paste0(output_dir, '/step1_ambient/',sample_name[i]), out)
        saveRDS(sc, paste(output_dir,'/step2/info2/',sample_name[i],'_ambient_rna_summary.rds', sep=''))
        #######################
        datadirs <- file.path(paste0(output_dir, '/step1_ambient/',sample_name[i]))
        sparse_matrix <- Seurat::Read10X(data.dir = datadirs)
        seurat_object <- Seurat::CreateSeuratObject(counts = sparse_matrix$`Gene Expression`, min.cells=min.cells_L, min.features=min.features_L)
        seurat_object[['HTO']] = Seurat::CreateAssayObject(counts = sparse_matrix$`Antibody Capture`, min.cells=min.cells_L, min.features=min.features_L)
        # seurat_object <- Seurat::CreateSeuratObject(counts = sparse_matrix$`Gene Expression`)
        # seurat_object[['HTO']] = Seurat::CreateAssayObject(counts = sparse_matrix$`Antibody Capture`)
        nam <- paste("seurat_object", sample_name[i], sep = ".")
        assign(nam, seurat_object)
        saveRDS(get(nam),paste(output_dir,'/step2/objs2',"/seu",i,".rds", sep=""),compress=TRUE)
        seu<-get(nam)
        seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")
        seu <- subset(seu, subset = percent.mt < 100)
        print(i)
        # png(file = paste(output_dir,'/step2/figs/vioplot',sample_name$V1[i],".png", sep=""))
        Seurat::VlnPlot(seu, group.by= "orig.ident", features = c("nFeature_RNA","nCount_RNA","percent.mt"), pt.size = 0.1,ncol = 3) + NoLegend()
        # dev.off()
        ggsave(paste(output_dir,'/step2/figs2/vioplot_seu',i,".png", sep=""))
        write.csv(colnames(seu[[]]), file= paste(output_dir,'/step2/info2/meta_info_seu',i,".txt", sep=""))
        #
        sink(paste(output_dir,'/step2/info2/summary_seu',i,".txt", sep=""))
        cat("Summary of nCount_RNA: \n")
        print(summary(seu$nCount_RNA))
        cat("Summary of nFeature_RNA: \n")
        print(summary(seu$nFeature_RNA))
        cat("Summary of pt_mito: \n")
        print(summary(seu$percent.mt))
        cat("The number of GEM/barcodes: \n")
        print(dim(seu))
        sink()
        #
        if (tolower(Save_RNA)=='yes') {
            mat <- GetAssayData(object = seu, assay = "RNA", slot = "data")
          #  write.csv(mat, paste(output_dir,'/step2/info2/',sample_name[i],"_RNA.csv", sep=""))
          writeMM(mat,file= paste(output_dir,'/step2/info2/seu',i,"_RNA.txt", sep=""))
          }
        if (tolower(Save_metadata)=='yes') {
          write.csv(seu[[]], file = paste(output_dir,'/step2/info2/seu_MetaData',i,'.txt', sep=""), quote = TRUE, sep = ",")
          }
        }
      # writeLines(capture.output(sessionInfo()), paste(output_dir,'/step2/info2/sessionInfo.txt', sep=""))
    } 
if (tolower(ambient_RNA)=="no" & !exists("count_matrices")) {
   ##########################
        foreach (i=1:length(sample_name)) %do% {    
          datadirs <- file.path(list[i],   "ouput_folder","outs","raw_feature_bc_matrix")
          names(datadirs)=sample_name[i]
          sparse_matrix <- Seurat::Read10X(data.dir = datadirs)
          seurat_object <- Seurat::CreateSeuratObject(counts = sparse_matrix$`Gene Expression`, min.cells=min.cells_L, min.features=min.features_L)
          seurat_object[['HTO']] = Seurat::CreateAssayObject(counts = sparse_matrix$`Antibody Capture`, min.cells=min.cells_L, min.features=min.features_L)
          # seurat_object <- Seurat::CreateSeuratObject(counts = sparse_matrix$`Gene Expression`)
          # seurat_object[['HTO']] = Seurat::CreateAssayObject(counts = sparse_matrix$`Antibody Capture`)
          nam <- paste("seurat_object", sample_name[i], sep = ".")
          assign(nam, seurat_object)
          saveRDS(get(nam),paste(output_dir,'/step2/objs2',"/seu",i,".rds", sep=""),compress=TRUE)
          seu<-get(nam)
          seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")
          seu <- subset(seu, subset = percent.mt < 100)
          print(i)
          # png(file = paste(output_dir,'/step2/figs/vioplot',sample_name$V1[i],".png", sep=""))
          Seurat::VlnPlot(seu, group.by= "orig.ident", features = c("nFeature_RNA","nCount_RNA","percent.mt"), pt.size = 0.1,ncol = 3) + NoLegend()
          # dev.off()
          ggsave(paste(output_dir,'/step2/figs2/vioplot_seu',i,".png", sep=""))
          write.csv(colnames(seu[[]]), file= paste(output_dir,'/step2/info2/meta_info_seu',i,".txt", sep=""))
          #
          sink(paste(output_dir,'/step2/info2/summary_seu',i,".txt", sep=""))
          cat("Summary of nCount_RNA: \n")
          print(summary(seu$nCount_RNA))
          cat("Summary of nFeature_RNA: \n")
          print(summary(seu$nFeature_RNA))
          cat("Summary of pt_mito: \n")
          print(summary(seu$percent.mt))
          cat("The number of GEM/barcodes: \n")
          print(dim(seu))
          sink()
          if (tolower(Save_RNA)=='yes') {
              mat <- GetAssayData(object = seu, assay = "RNA", slot = "data")
            #  write.csv(mat, paste(output_dir,'/step2/info2/',sample_name[i],"_RNA.csv", sep=""))
            writeMM(mat,file= paste(output_dir,'/step2/info2/seu',i,"_RNA.txt", sep=""))
            }
          if (tolower(Save_metadata)=='yes') {
            write.csv(seu[[]], file = paste(output_dir,'/step2/info2/seu_MetaData',i,'.txt', sep=""), quote = TRUE, sep = ",")
            }
          }
  # writeLines(capture.output(sessionInfo()), paste(output_dir,'/step2/info2/sessionInfo.txt', sep=""))
}


writeLines(capture.output(sessionInfo()), paste(output_dir,'/step2/info2/sessionInfo.txt', sep=""))
if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}