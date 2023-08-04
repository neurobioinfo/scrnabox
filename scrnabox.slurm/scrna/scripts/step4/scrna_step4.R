#!/usr/bin/env Rscript

####################
# step5

args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
.libPaths(r_lib_path)

packages<-c('Seurat','ggplot2', 'dplyr', 'foreach', 'doParallel','DoubletFinder')
lapply(packages, library, character.only = TRUE)
# library(DoubletFinder)
# output_dir="/home/samamiri/scratch/Darkgenome_run/pipeline/scrnabox.svn_run/des"
# list<-read.csv(paste(output_dir, "/job_output/configs/sample_dir.list",sep=""),header=FALSE)
# sample_name<-read.csv(paste(output_dir, "/job_output/configs/sample.list",sep=""),header=FALSE)
# sample_name<-list.files(path = paste(output_dir, "/step4/objs",sep=""),pattern = "*.rds")
source(paste(output_dir,'/job_info/parameters/step4_par.txt',sep=""))
sample_name<-list.files(path = paste(output_dir, "/step3/objs3",sep=""),pattern = "*.rds")
sample_nameb<-gsub(".rds","",sample_name)
# sample_nameb<-gsub(".RDS","",sample_name)
if(length(sample_name)<1) {
   print("You do not have any object from step 2 ")
}

for (i in 1:length(sample_name)) {
  if (!grepl(".rds",tolower(sample_name[i]), fixed = TRUE)){
     print(c(sample_name[i],"is not R rds"))
  }
}   

old.names<-old_label
new.names<-new_label

# library(foreach)
# library(doParallel)
numCores <- detectCores()
cl <- makeCluster(numCores-1)
registerDoParallel(cl) 

seu_list<-list()




if (tolower(dropDN)=='yes') {
    print('The following are deleted')
    print(label_dropDN)
    seu_list<-foreach (i_s=1:length(sample_name)) %do% {  
        seu<-readRDS(paste(output_dir,'/step3/objs3/',sample_name[i_s], sep=""))
        if (tolower(normlazation_and_scalaing)=='yes'){
            seu <- Seurat::NormalizeData(seu,normalization.method = par_normalization.method,scale.factor =par_scale.factor)
            seu<- FindVariableFeatures(seu, selection.method = par_selection.method, nfeatures = par_nfeatures)
            topsel <- head(Seurat::VariableFeatures(seu), par_top)
            write.csv(topsel, file = paste(output_dir,'/step4/info4/most_variable_genes_',sample_nameb[i],'.txt', sep=""), quote = TRUE, sep = ",")
            vf_plot <- Seurat::VariableFeaturePlot(seu)
            Seurat::LabelPoints(plot = vf_plot,points = topsel, repel = TRUE)
            ggsave(paste(output_dir,'/step4/figs4/VariableFeaturePlot',sample_nameb[i_s],".png",sep=""))
            seu<- ScaleData(seu, verbose = FALSE)
        }
        ### dimensionality reduction
        if (tolower(dimensionality_reduction)=='yes'){
            seu <- RunPCA(seu, npcs = par_npcs_pca, verbose = FALSE)
            DimPlot(seu, reduction = "pca")
            ggsave(paste(output_dir,'/step4/info4/',"dimplot_pca",sample_nameb[i_s],".png",sep=""))
            ElbowPlot(seu)
            ggsave(paste(output_dir,'/step4/info4/',"elbowplot",sample_nameb[i_s],".png",sep=""))
            Seurat::DimHeatmap(seu, dims = 1:par_dims, cells = par_cells, balanced = TRUE)
            ggsave(paste(output_dir,'/step4/info4/',"dimheatplot.",sample_nameb[i_s],".png",sep=""))
            # FeaturePlot(seu, reduction = par_reduction, features = par_features)
            # ggsave(paste(output_dir,'/step4/info4/',"featureplot",sample_name[i_s],".png",sep=""))
            seu <- RunUMAP(seu, dims = 1:par_dims_umap, n.neighbors =par_n.neighbors)
            Seurat::DimPlot(seu, reduction = "umap")
            ggsave(paste(output_dir,'/step4/figs4/',"dimplot_umap",sample_nameb[i_s],".png",sep=""))
        }
        sweep.res.list_pbmc <- paramSweep_v3(seu, PCs = 1:par_PCs, sct = par_sct)
        sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
        bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
        # select the pK that corresponds to max bcmvn to optimize doublet detection
        pK <- bcmvn_pbmc %>% filter(BCmetric == max(BCmetric)) %>%select(pK) 
        pK <- as.numeric(as.character(pK[[1]]))
        ## Homotypic Doublet Proportion Estimate ------------------------------------------------------------------------------------------------
        annotations <- seu@meta.data$seurat_clusters
        homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
        nExp_poi <- round(par_rate_nExp*nrow(seu@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
        nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
        # run doubletFinder 
        seu <- doubletFinder_v3(seu, PCs = 1:par_PCs, pN = par_pN,pK = pK, nExp = nExp_poi.adj,reuse.pANN = FALSE, sct = par_sct)
        # table(seu$DF.classifications_0.25_0.2_251)
        # colnames(seu[[]])
        # 'DF.classifications'%in%colnames(seu[[]])
        DF.classifications=colnames(seu[[]])[which(grepl('DF.classifications', colnames(seu[[]]), fixed=TRUE))]
        # unique(seu$DF.classifications_0.25_0.2_254)
        # unique(seu$DF.classifications_0.25_0.2_254)
        DF.classifications_u=eval(parse(text = paste('unique(seu$',DF.classifications,')', sep='')))
        # DF.classifications_u=unique(c(seu[[DF.classifications]]))        
        DimPlot(seu, reduction = 'umap', group.by = DF.classifications)
        ggsave(paste(output_dir,'/step4/figs4/',sample_nameb[i_s],"DF.classifications.png",sep=""))
        aa<-as.character(DF.classifications_u)[!as.character(DF.classifications_u) %in% "Doublet"]
        Idents(seu) <- DF.classifications
        seu=subset(seu,idents=aa)
################
        seu[["MULTI_ID"]] <- seu$orig.ident
        multi.names <- unique(seu@meta.data$MULTI_ID)
        Idents(seu)  <- "MULTI_ID"
        levels(seu@meta.data$MULTI_ID)
        for (i in 1:length(old.names)){
            newIdent <- new.names[i]
            names(newIdent) <- old.names[i]
              try( seu <- RenameIdents(object = seu, newIdent), silent = TRUE)  
        }
        seu[["MULTI_ID_Lables"]] <- Idents(seu)
###############       
        write.csv(table(seu$MULTI_ID), paste(output_dir,'/step4/info4/',sample_nameb[i_s],"MULTI_ID.csv",sep="")) 
        Idents(seu) <- "MULTI_ID"
        RidgePlot(seu, assay = "RNA", features = rownames(seu[["RNA"]]), ncol =par_RidgePlot_ncol) 
        ggsave(paste(output_dir,'/step4/figs4/',sample_nameb[i_s],"RidgeplotRNA.png",sep=""))
        # Idents(seu) <- "HTO_classification.global"
        VlnPlot(seu, features = "nCount_RNA", pt.size = 0.01, log = TRUE)
        ggsave(paste(output_dir,'/step4/figs4/seu',i_s,"VlnPlot_nCounts_RNA.png",sep=""))
        DoHeatmap(seu, features = rownames(seu[["RNA"]]))
        ggsave(paste(output_dir,'/step4/figs4/seu',i_s,"HeatmapRNA.png",sep=""))
        # DotPlot(seu, group.by = "MULTI_ID", features = rownames(seu[["HTO"]])) + theme(axis.text.x = element_text(angle = 90)) # Done 
        # DotPlot(seu, features = rownames(seu[["RNA"]]), group.by = "MULTI_ID")  # sug by MF
        # ggsave(paste(output_dir,'/step4/figs4/seu',i_s,"DotPlotRNA.png",sep=""))
        # multi.names <- unique(seu@meta.data$MULTI_ID)
        # Idents(seu)  <- "MULTI_ID"
        # levels(seu@meta.data$MULTI_ID)

        saveRDS(seu, paste(output_dir,'/step4/objs4/',sample_nameb[i_s],'.rds', sep=""))
        write.csv(colnames(seu[[]]), file= paste(output_dir,'/step4/info4/meta_info_',sample_nameb[i_s],".txt", sep=""))
        if (tolower(Save_RNA)=='yes') {
        mat <- GetAssayData(object = seu, assay = "RNA", slot = "data")
        #  write.csv(mat, paste(output_dir,'/step2/info2/',sample_name[i],"_RNA.csv", sep=""))
        writeMM(mat,file= paste(output_dir,'/step4/info4/',sample_nameb[i_s],"_RNA.txt", sep=""))
        }
        if (tolower(Save_metadata)=='yes') {
        write.csv(seu[[]], file = paste(output_dir,'/step4/info4/seu_MetaData',sample_nameb[i_s],'.txt', sep=""), quote = TRUE, sep = ",")
        }
    }
}


if (tolower(dropDN)=='no') {
    print('The droplet are not deleted')
    seu_list<-foreach (i_s=1:length(sample_name)) %do% {  
        seu<-readRDS(paste(output_dir,'/step3/objs3/',sample_name[i_s], sep=""))
        # seu <- seu
        if (tolower(normlazation_and_scalaing)=='yes'){
            seu <- Seurat::NormalizeData(seu,normalization.method = par_normalization.method,scale.factor =par_scale.factor)
            seu<- FindVariableFeatures(seu, selection.method = par_selection.method, nfeatures = par_nfeatures)
            topsel <- head(Seurat::VariableFeatures(seu), par_top)
            write.csv(topsel, file = paste(output_dir,'/step4/info4/most_variable_genes_',sample_nameb[i],'.txt', sep=""), quote = TRUE, sep = ",")
            vf_plot <- Seurat::VariableFeaturePlot(seu)
            Seurat::LabelPoints(plot = vf_plot,points = topsel, repel = TRUE)
            ggsave(paste(output_dir,'/step4/figs4/VariableFeaturePlot',sample_nameb[i_s],".png",sep=""))
            seu<- ScaleData(seu, verbose = FALSE)
        }
        ### dimensionality reduction
        if (tolower(dimensionality_reduction)=='yes'){
            seu <- RunPCA(seu, npcs = par_npcs_pca, verbose = FALSE)
            DimPlot(seu, reduction = "pca")
            ggsave(paste(output_dir,'/step4/info4/',"dimplot_pca",sample_nameb[i_s],".png",sep=""))
            ElbowPlot(seu)
            ggsave(paste(output_dir,'/step4/info4/',"elbowplot",sample_nameb[i_s],".png",sep=""))
            # Seurat::DimHeatmap(seu, dims = 1:par_dims, cells = par_cells, balanced = TRUE)
            # ggsave(paste(output_dir,'/step4/info4/',"dimheatplot.",sample_nameb[i_s],".png",sep=""))
            # FeaturePlot(seu, reduction = par_reduction, features = par_features)
            # ggsave(paste(output_dir,'/step4/info4/',"featureplot",sample_name[i_s],".png",sep=""))
            seu <- RunUMAP(seu, dims = 1:par_dims_umap, n.neighbors =par_n.neighbors)
            Seurat::DimPlot(seu, reduction = "umap")
            ggsave(paste(output_dir,'/step4/figs4/',"dimplot_umap",sample_nameb[i_s],".png",sep=""))
        }
        ################
        seu[["MULTI_ID"]] <- seu$orig.ident
        multi.names <- unique(seu@meta.data$MULTI_ID)
        Idents(seu)  <- "MULTI_ID"
        levels(seu@meta.data$MULTI_ID)
        for (i in 1:length(old.names)){
            newIdent <- new.names[i]
            names(newIdent) <- old.names[i]
                try( seu <- RenameIdents(object = seu, newIdent), silent = TRUE)  
        }
        seu[["MULTI_ID_Lables"]] <- Idents(seu)
        ###############    
        write.csv(table(seu$MULTI_ID), paste(output_dir,'/step4/info4/',sample_nameb[i_s],"MULTI_ID.csv",sep="")) 
        Idents(seu) <- "MULTI_ID"
        RidgePlot(seu, assay = "RNA", features = rownames(seu[["RNA"]]), ncol =par_RidgePlot_ncol) 
        ggsave(paste(output_dir,'/step4/figs4/',sample_nameb[i_s],"RidgeplotRNA.png",sep=""))
        # Idents(seu) <- "HTO_classification.global"
        VlnPlot(seu, features = "nCount_RNA", pt.size = 0.01, log = TRUE)
        ggsave(paste(output_dir,'/step4/figs4/seu',i_s,"VlnPlot_nCounts_RNA.png",sep=""))
        DoHeatmap(seu, features = rownames(seu[["RNA"]]))
        ggsave(paste(output_dir,'/step4/figs4/seu',i_s,"HeatmapRNA.png",sep=""))
        # DotPlot(seu, group.by = "MULTI_ID", features = rownames(seu[["HTO"]])) + theme(axis.text.x = element_text(angle = 90)) # Done 
        # DotPlot(seu, features = rownames(seu[["RNA"]]), group.by = "MULTI_ID")  # sug by MF
        # ggsave(paste(output_dir,'/step4/figs4/seu',i_s,"DotPlotRNA.png",sep=""))
        # multi.names <- unique(seu@meta.data$MULTI_ID)
        # Idents(seu)  <- "MULTI_ID"
        # levels(seu@meta.data$MULTI_ID)    
        saveRDS(seu, paste(output_dir,'/step4/objs4/',sample_nameb[i_s],'.rds', sep=""))
        write.csv(colnames(seu[[]]), file= paste(output_dir,'/step4/info4/meta_info_',sample_nameb[i_s],".txt", sep=""))
        if (tolower(Save_RNA)=='yes') {
            mat <- GetAssayData(object = seu, assay = "RNA", slot = "data")
            #  write.csv(mat, paste(output_dir,'/step2/info2/',sample_name[i],"_RNA.csv", sep=""))
            writeMM(mat,file= paste(output_dir,'/step4/info4/',sample_nameb[i_s],"_RNA.txt", sep=""))
        }
        if (tolower(Save_metadata)=='yes') {
            write.csv(seu[[]], file = paste(output_dir,'/step4/info4/seu_MetaData',sample_nameb[i_s],'.txt', sep=""), quote = TRUE, sep = ",")
        }
        }
}


writeLines(capture.output(sessionInfo()), paste(output_dir,'/step4/info4/sessionInfo.txt', sep=""))

if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}
