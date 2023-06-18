#!/usr/bin/env Rscript

####################
# step2 

args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]

.libPaths(r_lib_path)

packages<-c('Seurat','ggplot2', 'dplyr')
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

library(foreach)
library(doParallel)
numCores <- detectCores()
cl <- makeCluster(numCores-1)
registerDoParallel(cl) 

seu_list<-list()

seu_list<-foreach (i_s=1:length(sample_name)) %do% {  
    seu<-readRDS(paste(output_dir,'/step4/objs4/',sample_name[i_s], sep=""))
    DefaultAssay(seu) <- "RNA"
    seu <- Seurat::NormalizeData(seu)
    seu <- Seurat::FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
}  

seu_anchors <- Seurat::FindIntegrationAnchors(object.list = seu_list)
seu_int <- Seurat::IntegrateData(anchorset = seu_anchors)
Seurat::DefaultAssay(seu_int) <- "integrated"
saveRDS(seu_int, paste(output_dir,'/step5/objs5',"/seu_step5.rds", sep=""))
write.csv(colnames(seu_int[[]]), paste(output_dir,'/step5/info5',"/meta_info.csv", sep=""))
write.csv(table(seu_int$MULTI_ID), paste(output_dir,'/step5/info5',"/MULTI_IDcount.csv", sep=""))

if (tolower(Save_RNA)=='yes') {
    mat <- GetAssayData(object = seu_int, assay = "RNA", slot = "data")
    writeMM(mat,file= paste(output_dir,'/step5/info5/seu_int',"_RNA.txt", sep=""))
}

if (tolower(Save_metadata)=='yes') {
    write.csv(seu_int[[]], file = paste(output_dir,'/step5/info5/seu_int_MetaData.txt', sep=""), quote = TRUE, sep = ",")
}

writeLines(capture.output(sessionInfo()), paste(output_dir,'/step5/info5/sessionInfo.txt', sep=""))

