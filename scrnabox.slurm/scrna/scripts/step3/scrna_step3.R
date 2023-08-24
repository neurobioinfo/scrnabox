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



.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr','Matrix', 'foreach', 'doParallel')
lapply(packages, library, character.only = TRUE)
# output_dir="/home/samamiri/scratch/Darkgenome_run/pipeline/scrnabox.svn_run/des"
# list<-read.csv(paste(output_dir, "/job_output/configs/sample_dir.list",sep=""),header=FALSE)
# sample_name<-read.csv(paste(output_dir, "/job_output/configs/sample.list",sep=""),header=FALSE)
sample_name<-list.files(path = paste(output_dir, "/step2/objs2",sep=""))
sample_nameb<-gsub(".rds","",sample_name)
# sample_nameb<-gsub(".RDS","",sample_nameb)
if(length(sample_name)<1) {
   print("You do not have any object from step 2 ")
}

for (i in 1:length(sample_name)) {
  if (!grepl(".rds",tolower(sample_name[i]), fixed = TRUE)){
     print(c(sample_name[i],"is not R rds"))
  }
}   

source(paste(output_dir,'/job_info/parameters/step3_par.txt',sep=""))
# save.image(file = "/lustre03/project/6070393/COMMON/Dark_Genome/analysis_DarkGenome2weeks_SCRNA/temp.RData")


if (exists("nFeature_RNA_L")) NFRNAL=nFeature_RNA_L
if (exists("nFeature_RNA_U")) NFRNAU=nFeature_RNA_U

if (exists("nCount_RNA_L")) NCRNAL=nCount_RNA_L
if (exists("nCount_RNA_U")) NCRNAU=nCount_RNA_U

if (exists("mitochondria_percent_L")) PMTL=mitochondria_percent_L
if (exists("mitochondria_percent_U")) PMTU=mitochondria_percent_U

if (exists("log10GenesPerUMI_L")) GENEUMIL=log10GenesPerUMI_L
if (exists("log10GenesPerUMI_U")) GENEUMIU=log10GenesPerUMI_U

print(output_dir)
print(NFRNAL)
print(NFRNAU)
print(NCRNAL)
print(NCRNAU)
print(PMTL)
print(PMTU)
print(GENEUMIL)
print(GENEUMIU)

set.seed(1234)
# library(foreach)
# library(doParallel)

numCores <- detectCores()
cl <- makeCluster(numCores-1)
registerDoParallel(cl) 
# i<-1
foreach (i=1:length(sample_name)) %do% {    
    set.seed(1234)
    seu<-readRDS(paste(output_dir,'/step2/objs2/',sample_name[i], sep=""))
    # DA.hash<-seu
    # hashtags <- rownames(DA.hash[["HTO"]])
    # seu <-  DA.hash
    seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")
    if (is.na(NFRNAL)) {
        NFRNAL=min(seu[["nFeature_RNA"]])
    }
    if (is.na(NFRNAU)) {
        NFRNAU=max(seu[["nFeature_RNA"]])
    }
    if (is.na(NCRNAL)) {
        NCRNAL=min(seu[["nCount_RNA"]])
    }
    if (is.na(NCRNAU)) {
        NCRNAU=max(seu[["nCount_RNA"]])
    }
    if (is.na(PMTL)) {
        PMTL=min(seu[["percent.mt"]])
    }
    if (is.na(PMTU)) {
        PMTU=max(seu[["percent.mt"]])
    }
    seu$log10GenesPerUMI <- log10(seu$nFeature_RNA) / log10(seu$nCount_RNA)
    if (is.na(GENEUMIL)) {
        GENEUMIL=min(seu[["log10GenesPerUMI"]])
    }
    if (is.na(GENEUMIU)) {
        GENEUMIU=max(seu[["log10GenesPerUMI"]])
    }

    seu <- subset(seu, subset = nFeature_RNA > NFRNAL & nFeature_RNA < NFRNAU & nCount_RNA > NCRNAL & nCount_RNA < NCRNAU & PMTL < percent.mt & percent.mt < PMTU & log10GenesPerUMI > GENEUMIL & log10GenesPerUMI < GENEUMIU)
    # seu <- subset(seu, subset = nFeature_RNA > NFRNAL & nCount_RNA < NFRNAU & percent.mt < PMT)
    # seu <- subset(seu, subset = nFeature_RNA > 300 & nFeature_RNA < 6500 & percent.mt < 25)
    saveRDS(seu, paste(output_dir,'/step3/objs3/',sample_nameb[i],".rds", sep=""))
    Seurat::VlnPlot(seu, group.by= "orig.ident", features = c("nFeature_RNA","nCount_RNA","percent.mt"), pt.size = 0.1,ncol = 3) + NoLegend()
    ggsave(paste(output_dir,'/step3/figs3/vioplot_',sample_nameb[i],".png", sep=""))
    write.csv(colnames(seu[[]]), file= paste(output_dir,'/step3/info3/meta_info_',sample_nameb[i],".txt", sep=""))
    if (tolower(Save_RNA)=='yes') {
       mat <- GetAssayData(object = seu, assay = "RNA", slot = "data")
      #  write.csv(mat, paste(output_dir,'/step2/info2/',sample_name[i],"_RNA.csv", sep=""))
      writeMM(mat,file= paste(output_dir,'/step3/info3/',sample_nameb[i],"_RNA.txt", sep=""))
    }
    if (tolower(Save_metadata)=='yes') {
      write.csv(seu[[]], file = paste(output_dir,'/step3/info3/MetaData',sample_nameb[i],'.txt', sep=""), quote = TRUE, sep = ",")
    }
    sink(paste(output_dir,'/step3/info3/summary_',sample_nameb[i],".txt", sep=""))
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
if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}

