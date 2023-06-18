#!/usr/bin/env Rscript

####################
# step2 

args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
# scrna_method=args[3]
.libPaths(r_lib_path)



packages<-c('Seurat','ggplot2', 'dplyr', 'foreach', 'doParallel')
lapply(packages, library, character.only = TRUE)

# output_dir="/home/samamiri/scratch/Darkgenome_run/pipeline/scrnabox.svn_run/des"
# list<-read.csv(paste(output_dir, "/job_output/configs/sample_dir.list",sep=""),header=FALSE)
# sample_name<-read.csv(paste(output_dir, "/job_output/configs/sample.list",sep=""),header=FALSE)
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
    DefaultAssay(seu) <- "RNA"
    seu <- Seurat::NormalizeData(seu)
    seu <- Seurat::FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
    # aa<-as.character(unique(seu$MULTI_ID))[!as.character(unique(seu$MULTI_ID)) %in% label_dropDN]
    # Idents(seu) <- "MULTI_ID"
    # seu=subset(seu,idents=aa)
}  


# if (dropDN=='yes') {
#     print('The following are deleted')
#     print(label_dropDN)
#     seu_list<-foreach (i_s=1:length(sample_name)) %do% {  
#         seu<-readRDS(paste(output_dir,'/step4/objs/',sample_name[i_s], sep=""))
#         DefaultAssay(seu) <- "RNA"
#         seu <- Seurat::NormalizeData(seu)
#         seu <- Seurat::FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
#         aa<-as.character(unique(seu$MULTI_ID))[!as.character(unique(seu$MULTI_ID)) %in% label_dropDN]
#         Idents(seu) <- "MULTI_ID"
#         seu=subset(seu,idents=aa)
#     }  
# } else {
#     seu_list<-foreach (i_s=1:dim(sample_name)[1]) %do% {  
#         seu<-readRDS(paste(output_dir,'/step4/objs/',sample_name[i_s], sep=""))
#         DefaultAssay(seu) <- "RNA"
#         seu <- Seurat::NormalizeData(seu)
#         seu <- Seurat::FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
#     }  
# }


seu_anchors <- Seurat::FindIntegrationAnchors(object.list = seu_list)
seu_int <- Seurat::IntegrateData(anchorset = seu_anchors)
Seurat::DefaultAssay(seu_int) <- "integrated"
saveRDS(seu_int, paste(output_dir,'/step5/objs5',"/seu_step5.rds", sep=""))
write.csv(colnames(seu_int[[]]), paste(output_dir,'/step5/info5',"/meta_info.csv", sep=""))
if (scrna_method=='HTO') {
    write.csv(table(seu_int$MULTI_ID), paste(output_dir,'/step5/info5',"/MULTI_IDcount.csv", sep=""))
}

if (tolower(Save_RNA)=='yes') {
    mat <- GetAssayData(object = seu_int, assay = "RNA", slot = "data")
    #  write.csv(mat, paste(output_dir,'/step2/info2/',sample_name[i],"_RNA.csv", sep=""))
    writeMM(mat,file= paste(output_dir,'/step5/info5/seu_int',"_RNA.txt", sep=""))
}

if (tolower(Save_metadata)=='yes') {
    write.csv(seu_int[[]], file = paste(output_dir,'/step5/info5/seu_int_MetaData.txt', sep=""), quote = TRUE, sep = ",")
}

writeLines(capture.output(sessionInfo()), paste(output_dir,'/step5/info5/sessionInfo.txt', sep=""))
