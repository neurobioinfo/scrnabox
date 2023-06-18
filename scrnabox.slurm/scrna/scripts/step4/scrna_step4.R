#!/usr/bin/env Rscript

####################
# step5

args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
.libPaths(r_lib_path)

packages<-c('Seurat','ggplot2', 'dplyr')
lapply(packages, library, character.only = TRUE)
library(DoubletFinder)
# output_dir="/home/samamiri/scratch/Darkgenome_run/pipeline/scrnabox.svn_run/des"
# list<-read.csv(paste(output_dir, "/job_output/configs/sample_dir.list",sep=""),header=FALSE)
# sample_name<-read.csv(paste(output_dir, "/job_output/configs/sample.list",sep=""),header=FALSE)
# sample_name<-list.files(path = paste(output_dir, "/step4/objs",sep=""),pattern = "*.rds")
source(paste(output_dir,'/job_info/parameters/step4_par.txt',sep=""))
sample_name<-list.files(path = paste(output_dir, "/step3/objs3",sep=""),pattern = "*.rds")

if(length(sample_name)<1) {
   print("You do not have any object from step 2 ")
}

for (i in 1:length(sample_name)) {
  if (!grepl(".rds",tolower(sample_name[i]), fixed = TRUE)){
     print(c(sample_name[i],"is not R rds"))
  }
}   


library(foreach)
library(doParallel)
numCores <- detectCores()
cl <- makeCluster(numCores-1)
registerDoParallel(cl) 

seu_list<-list()

if (dropDN=='yes') {
    print('The following are deleted')
    print(label_dropDN)
    seu_list<-foreach (i_s=1:length(sample_name)) %do% {  
        seu<-readRDS(paste(output_dir,'/step3/objs3/',sample_name[i_s], sep=""))
        pbmc.seurat.filtered <- NormalizeData(object = seu)
        pbmc.seurat.filtered <- FindVariableFeatures(object = pbmc.seurat.filtered)
        pbmc.seurat.filtered <- ScaleData(object = pbmc.seurat.filtered)
        pbmc.seurat.filtered <- RunPCA(object = pbmc.seurat.filtered)
        pbmc.seurat.filtered <- FindNeighbors(object = pbmc.seurat.filtered, dims = 1:20)
        pbmc.seurat.filtered <- FindClusters(object = pbmc.seurat.filtered)
        pbmc.seurat.filtered <- RunUMAP(object = pbmc.seurat.filtered, dims = 1:20)
        sweep.res.list_pbmc <- paramSweep_v3(pbmc.seurat.filtered, PCs = 1:20, sct = FALSE)
        sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
        bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
        # select the pK that corresponds to max bcmvn to optimize doublet detection
        pK <- bcmvn_pbmc %>% filter(BCmetric == max(BCmetric)) %>%select(pK) 
        pK <- as.numeric(as.character(pK[[1]]))
        ## Homotypic Doublet Proportion Estimate ------------------------------------------------------------------------------------------------
        annotations <- pbmc.seurat.filtered@meta.data$seurat_clusters
        homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
        nExp_poi <- round(0.076*nrow(pbmc.seurat.filtered@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
        nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
        # run doubletFinder 
        pbmc.seurat.filtered <- doubletFinder_v3(pbmc.seurat.filtered, PCs = 1:20, pN = 0.25,pK = pK, nExp = nExp_poi.adj,reuse.pANN = FALSE, sct = FALSE)
        # table(pbmc.seurat.filtered$DF.classifications_0.25_0.2_251)
        # colnames(pbmc.seurat.filtered[[]])
        # 'DF.classifications'%in%colnames(pbmc.seurat.filtered[[]])
        DF.classifications=colnames(pbmc.seurat.filtered[[]])[which(grepl('DF.classifications', colnames(pbmc.seurat.filtered[[]]), fixed=TRUE))]
        # unique(pbmc.seurat.filtered$DF.classifications_0.25_0.2_254)
        # unique(pbmc.seurat.filtered$DF.classifications_0.25_0.2_254)
        DF.classifications_u=eval(parse(text = paste('unique(pbmc.seurat.filtered$',DF.classifications,')', sep='')))
        # DF.classifications_u=unique(c(pbmc.seurat.filtered[[DF.classifications]]))        
        DimPlot(pbmc.seurat.filtered, reduction = 'umap', group.by = DF.classifications)
        ggsave(paste(output_dir,'/step4/figs4/seu',i_s,"DF.classifications.png",sep=""))
        aa<-as.character(DF.classifications_u)[!as.character(DF.classifications_u) %in% "Doublet"]
        Idents(pbmc.seurat.filtered) <- DF.classifications
        pbmc.seurat.filtered=subset(pbmc.seurat.filtered,idents=aa)
        saveRDS(pbmc.seurat.filtered, paste(output_dir,'/step4/objs4/seu',i_s,'.rds', sep=""))
        write.csv(colnames(pbmc.seurat.filtered[[]]), file= paste(output_dir,'/step4/info4/meta_info_seu',i_s,".txt", sep=""))
        if (tolower(Save_RNA)=='yes') {
        mat <- GetAssayData(object = seu1, assay = "RNA", slot = "data")
        #  write.csv(mat, paste(output_dir,'/step2/info2/',sample_name[i],"_RNA.csv", sep=""))
        writeMM(mat,file= paste(output_dir,'/step4/info4/seu',i_s,"_RNA.txt", sep=""))
        }
        if (tolower(Save_metadata)=='yes') {
        write.csv(seu1[[]], file = paste(output_dir,'/step4/info4/seu_MetaData',i_s,'.txt', sep=""), quote = TRUE, sep = ",")
        }
    }
}
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step4/info4/sessionInfo.txt', sep=""))

