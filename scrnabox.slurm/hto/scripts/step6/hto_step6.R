#!/usr/bin/env Rscript

####################
# step6

args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
pipeline_home=args[3]


.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr','cowplot','clustree')
lapply(packages, library, character.only = TRUE)
# library('cowplot')
# library('clustree')

# sample_name<-list.files(path = paste(output_dir, "/step5/objs5",sep=""),pattern = "*.rds")
# if(length(sample_name)<1) {
#    print("You do not have any object from step 2 ")
# }
# for (i in 1:length(sample_name)) {
#   if (!grepl(".rds",tolower(sample_name[i]), fixed = TRUE)){
#      print(c(sample_name[i],"is not R rds"))
#   }
# }   
# seu_int<-readRDS(paste(output_dir,'/step5/objs5/',sample_name, sep=''))

source(paste(output_dir,'/job_info/parameters/step6_par.txt',sep=""))

if (tolower(skip_step5)=='yes') {
    sample_name<-list.files(path = paste(output_dir, "/step4/objs4",sep=""),pattern = "*.rds")
    if(length(sample_name)<1) {
    print("You do not have any object from step 2 ")
    }
    for (i in 1:length(sample_name)) {
    if (!grepl(".rds",tolower(sample_name[i]), fixed = TRUE)){
        print(c(sample_name[i],"is not R rds"))
    }
    }   
    seu_int<-readRDS(paste(output_dir,'/step4/objs4/',sample_name, sep=''))
    par_whatAssay<-'RNA'
} else {
    sample_name<-list.files(path = paste(output_dir, "/step5/objs5",sep=""),pattern = "*.rds")
    if(length(sample_name)<1) {
    print("You do not have any object from step 2 ")
    }
    for (i in 1:length(sample_name)) {
    if (!grepl(".rds",tolower(sample_name[i]), fixed = TRUE)){
        print(c(sample_name[i],"is not R rds"))
    }
    }   
    seu_int<-readRDS(paste(output_dir,'/step5/objs5/',sample_name, sep=''))
    par_whatAssay<-'integrated'
}


# source(paste(output_dir,'/job_info/parameters/step6_par.txt',sep=""))

# DefaultAssay(seu_int) <- "RNA"
# seu_int <- NormalizeData(seu_int)
# seu_int <- ScaleData(seu_int)
# # Do the PCA
# seu_int <- FindVariableFeatures(seu_int, selection.method = par_selection.method, nfeatures = par_nfeatures)
# seu_int <- RunPCA(seu_int)  # should use variable ft by default

# png(paste(output_dir,'/step6/figs6',"/elbow_pca.png", sep=""))
# fig1<-DimPlot(seu_int, reduction = "pca")
# fig2<-ElbowPlot(seu_int)
# print(fig1 + fig2)
# dev.off()



#######
Seurat::DefaultAssay(seu_int) <-  par_whatAssay #"integrated"
# seu_int <- ScaleData(seu_int, verbose = FALSE)
# seu_int <- RunPCA(seu_int, npcs = par_RunPCA_npcs, verbose = FALSE)
# seu_int <- RunUMAP(seu_int, dims = 1:par_RunUMAP_dims, n.neighbors =par_RunUMAP_n.neighbors)
seu_int <- FindNeighbors(seu_int,  dims = 1:par_FindNeighbors_dims, k.param = par_FindNeighbors_k.param, prune.SNN = par_FindNeighbors_prune.SNN)
seu_int <- Seurat::FindClusters(seu_int, resolution = par_FindClusters_resolution)


clustree(seu_int@meta.data, prefix = paste0(par_whatAssay,"_snn_res."))
ggsave(paste(output_dir,"/step6/figs6/clustree_int.png", sep=""))
for (i in par_FindClusters_resolution){
    Seurat::DimPlot(seu_int, group.by = paste((par_whatAssay,"_snn_res.",i,sep=''))
    ggsave(paste(output_dir,'/step6/figs6/',par_whatAssay,"_snn_res.",i,".png",sep=""))
}

########

saveRDS(seu_int, paste(output_dir,'/step6/objs6',"/seu_step6.rds", sep=""))
write.csv(colnames(seu_int[[]]), paste(output_dir,'/step6/info6',"/meta_info.csv", sep=""))

write.table(unique(seu_int$MULTI_ID_Lables),file=paste(output_dir,'/step6/info6',"/class_label_phenotypes.txt", sep=""),col.names=FALSE)


source(paste(pipeline_home,'/general_codes/Rand_index.R',sep=''))

if (tolower(Save_RNA)=='yes') {
    mat <- GetAssayData(object = seu_int, assay = "RNA", slot = "data")
    writeMM(mat,file= paste(output_dir,'/step6/info6/seu',"_RNA.txt", sep=""))
}

if (tolower(Save_metadata)=='yes') {
    write.csv(seu_int[[]], file = paste(output_dir,'/step6/info6/seu_MetaData.txt', sep=""), quote = TRUE, sep = ",")
}

writeLines(capture.output(sessionInfo()), paste(output_dir,'/step6/info6/sessionInfo.txt', sep=""))

if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}