#!/usr/bin/env Rscript

####################
# step5 

args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]

.libPaths(r_lib_path)

packages<-c('Seurat','ggplot2', 'dplyr','foreach', 'doParallel')
lapply(packages, library, character.only = TRUE)

sample_name<-list.files(path = paste(output_dir, "/step4/objs4",sep=""),pattern = "*.rds")
if(length(sample_name)<1) {
   print("You do not have any object from step 2 ")
}

for (i in 1:length(sample_name)) {
  if (!grepl(".rds",tolower(sample_name[i]), fixed = TRUE)){
     print(c(sample_name[i],"is not R rds"))
  }
}   

source(paste(output_dir,'/job_info/parameters/step5_par.txt',sep=""))

# library(foreach)
# library(doParallel)
numCores <- detectCores()
cl <- makeCluster(numCores-1)
registerDoParallel(cl) 

seu_list<-list()

seu_list<-foreach (i_s=1:length(sample_name)) %do% {  
    seu<-readRDS(paste(output_dir,'/step4/objs4/',sample_name[i_s], sep=""))
    DefaultAssay(seu) <- par_DefaultAssay
    if (tolower(normlazation_and_scalaing)=='yes'){
                seu <- Seurat::NormalizeData(seu,normalization.method = par_normalization.method,scale.factor =par_scale.factor)
                seu<- FindVariableFeatures(seu, selection.method = par_selection.method, nfeatures = par_nfeatures)
    }
}  

# seu_anchors <- Seurat::FindIntegrationAnchors(object.list = seu_list)
# seu_int <- Seurat::IntegrateData(anchorset = seu_anchors)
# Seurat::DefaultAssay(seu_int) <- "integrated"
# saveRDS(seu_int, paste(output_dir,'/step5/objs5',"/seu_step5.rds", sep=""))
# write.csv(colnames(seu_int[[]]), paste(output_dir,'/step5/info5',"/meta_info.csv", sep=""))
# write.csv(table(seu_int$MULTI_ID), paste(output_dir,'/step5/info5',"/MULTI_IDcount.csv", sep=""))

# if (tolower(Save_RNA)=='yes') {
#     mat <- GetAssayData(object = seu_int, assay = "RNA", slot = "data")
#     writeMM(mat,file= paste(output_dir,'/step5/info5/seu_int',"_RNA.txt", sep=""))
# }

# if (tolower(Save_metadata)=='yes') {
#     write.csv(seu_int[[]], file = paste(output_dir,'/step5/info5/seu_int_MetaData.txt', sep=""), quote = TRUE, sep = ",")
# }

# writeLines(capture.output(sessionInfo()), paste(output_dir,'/step5/info5/sessionInfo.txt', sep=""))




seu_anchors <- Seurat::FindIntegrationAnchors(object.list = seu_list,dims = 1:par_FindIntegrationAnchors_dim)
seu_int <- Seurat::IntegrateData(anchorset = seu_anchors,dims = 1:par_FindIntegrationAnchors_dim)
Seurat::DefaultAssay(seu_int) <- "integrated"


seu_int <- ScaleData(seu_int, verbose = FALSE)
seu_int <- RunPCA(seu_int, npcs = par_RunPCA_npcs, verbose = FALSE)
seu_int <- RunUMAP(seu_int, dims = 1:par_RunUMAP_dims, n.neighbors =par_RunUMAP_n.neighbors)


DimPlot(seu_int, reduction = "pca")
ggsave(paste(output_dir,'/step5/figs5',"/DimPlot_pca.png", sep=""))
ElbowPlot(seu_int)
ggsave(paste(output_dir,'/step5/figs5',"/elbow.png", sep=""))
DimPlot(seu_int, reduction = "umap")
ggsave(paste(output_dir,'/step5/figs5',"/DimPlot_umap.png", sep=""))


saveRDS(seu_int, paste(output_dir,'/step5/objs5',"/seu_step5.rds", sep=""))
write.csv(colnames(seu_int[[]]), paste(output_dir,'/step5/info5',"/meta_info_seu_step5.csv", sep=""))
write.csv(table(seu_int$MULTI_ID), paste(output_dir,'/step5/info5',"/MULTI_IDcount.csv", sep=""))


if (tolower(Save_RNA)=='yes') {
    mat <- GetAssayData(object = seu_int, assay = "RNA", slot = "data")
    #  write.csv(mat, paste(output_dir,'/step2/info2/',sample_name[i],"_RNA.csv", sep=""))
    writeMM(mat,file= paste(output_dir,'/step5/info5/seu_int',"_RNA.txt", sep=""))
}

if (tolower(Save_metadata)=='yes') {
    write.csv(seu_int[[]], file = paste(output_dir,'/step5/info5/seu_int_MetaData.txt', sep=""), quote = TRUE, sep = ",")
}

writeLines(capture.output(sessionInfo()), paste(output_dir,'/step5/info5/sessionInfo.txt', sep=""))

if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}



