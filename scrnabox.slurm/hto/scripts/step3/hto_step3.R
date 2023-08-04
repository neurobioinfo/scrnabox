#!/usr/bin/env Rscript

####################
# step3 

args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]

NFRNAL=as.numeric(args[3])
NFRNAU=as.numeric(args[4])
NCRNAL=as.numeric(args[5])
NCRNAU=as.numeric(args[6])
PMTL=as.numeric(args[7])
PMTU=as.numeric(args[8])
GENEUMIL=as.numeric(args[9])
GENEUMIU=as.numeric(args[10])

print(output_dir)
print(NFRNAL)
print(NFRNAU)
print(NCRNAL)
print(NCRNAU)
print(PMTL)
print(PMTU)
print(GENEUMIL)
print(GENEUMIU)

.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr','foreach', 'doParallel','Matrix')
lapply(packages, library, character.only = TRUE)


sample_name<-list.files(path = paste(output_dir, "/step2/objs2",sep=""))
if(length(sample_name)<1) {
   print("You do not have any object from step 2 ")
}

for (i in 1:length(sample_name)) {
  if (!grepl(".rds",tolower(sample_name[i]), fixed = TRUE)){
     print(c(sample_name[i],"is not R rds"))
  }
}   

source(paste(output_dir,'/job_info/parameters/step3_par.txt',sep=""))


set.seed(1234)

numCores <- detectCores()
cl <- makeCluster(numCores-1)
registerDoParallel(cl) 
# i<-1
foreach (i=1:length(sample_name)) %do% {    
    set.seed(1234)
    seu<-readRDS(paste(output_dir,'/step2/objs2/',sample_name[i], sep=""))
    seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")
    if (is.na(NFRNAL)) {
        NFRNAL=min(seu[["nFeature_RNA"]],na.rm =T)
    }
    if (is.na(NFRNAU)) {
        NFRNAU=max(seu[["nFeature_RNA"]],na.rm =T)
    }
    if (is.na(NCRNAL)) {
        NCRNAL=min(seu[["nCount_RNA"]],na.rm =T)
    }
    if (is.na(NCRNAU)) {
        NCRNAU=max(seu[["nCount_RNA"]],na.rm =T)
    }
    if (is.na(PMTL)) {
        PMTL=min(seu[["percent.mt"]],na.rm =T)
    }
    if (is.na(PMTU)) {
        PMTU=max(seu[["percent.mt"]],na.rm =T)
    }
    # seu$log10GenesPerUMI <- log10(seu$nFeature_RNA) / log10(seu$nCount_RNA)
    # if (is.na(GENEUMIL)) {
    #     GENEUMIL=min(seu[["log10GenesPerUMI"]])
    # }
    # if (is.na(GENEUMIU)) {
    #     GENEUMIU=max(seu[["log10GenesPerUMI"]],na.rm =T)
    # }

    seu <- subset(seu, subset = nFeature_RNA > NFRNAL & nFeature_RNA < NFRNAU & nCount_RNA > NCRNAL & nCount_RNA < NCRNAU & PMTL < percent.mt & percent.mt < PMTU )
    # seu <- subset(seu, subset = nFeature_RNA > NFRNAL & nFeature_RNA < NFRNAU & nCount_RNA > NCRNAL & nCount_RNA < NCRNAU & PMTL < percent.mt & percent.mt < PMTU & log10GenesPerUMI > GENEUMIL & log10GenesPerUMI < GENEUMIU)
    # seu <- subset(seu, subset = nFeature_RNA > NFRNAL & nFeature_RNA < NFRNAU & nCount_RNA > NCRNAL & nCount_RNA < NCRNAU & percent.mt < PMTU)
    saveRDS(seu, paste(output_dir,'/step3/objs3',"/seu",i,".rds", sep=""))
    Seurat::VlnPlot(seu, group.by= "orig.ident", features = c("nFeature_RNA","nCount_RNA","percent.mt"), pt.size = 0.1,ncol = 3) + NoLegend()
    ggsave(paste(output_dir,'/step3/figs3/vioplot_seu',i,".png", sep=""))
    write.csv(colnames(seu[[]]), file= paste(output_dir,'/step3/info3/meta_info_seu',i,".txt", sep=""))
    if (tolower(Save_RNA)=='yes') {
       mat <- GetAssayData(object = seu, assay = "RNA", slot = "data")
      writeMM(mat,file= paste(output_dir,'/step3/info3/seu',i,"_RNA.txt", sep=""))
    }
    if (tolower(Save_metadata)=='yes') {
      write.csv(seu[[]], file = paste(output_dir,'/step3/info3/seu_MetaData',i,'.txt', sep=""), quote = TRUE, sep = ",")
    }
    sink(paste(output_dir,'/step3/info3/summary_seu',i,".txt", sep=""))
    cat("Summary of nCount_RNA: \n")
    print(summary(seu$nCount_RNA))
    cat("Summary of nFeature_RNA: \n")
    print(summary(seu$nFeature_RNA))
    cat("Summary of pt_mito: \n")
    print(summary(seu$percent.mt))
    cat("The number of GEM/barcodes: \n")
    print(dim(seu))
    sink()
}

writeLines(capture.output(sessionInfo()), paste(output_dir,'/step3/info3/sessionInfo.txt', sep=""))
