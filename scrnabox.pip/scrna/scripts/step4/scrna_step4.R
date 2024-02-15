#!/usr/bin/env Rscript

##########################################
# step4: Doublet detection
##########################################

## load parameters
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
.libPaths(r_lib_path)

## load library
packages<-c('Seurat','ggplot2', 'dplyr', 'foreach', 'doParallel','DoubletFinder', 'Matrix', 'ggpubr', 'data.table')
lapply(packages, library, character.only = TRUE)


## load parameters
source(paste(output_dir,'/job_info/parameters/step4_par.txt',sep=""))

## load Seurat object
if (exists("par_seurat_object")) {                                                  
    sample_name<-list.files(path = par_seurat_object)
    sample_nameb<-gsub(".rds","",sample_name)
    if(length(sample_name)<1) {
    print("You do not have any existing Seurat object")
    }
for (i in 1:length(sample_name)) {
    if (!grepl(".rds",tolower(sample_name[i]), fixed = TRUE)){
        print(c(sample_name[i],"is not R rds"))
    }
    }  
} else {
    sample_name<-list.files(path = paste(output_dir, "/step3/objs3",sep=""))
    sample_nameb<-gsub(".rds","",sample_name)
    if(length(sample_name)<1) {
    print("You do not have any object from step 3 ")
    }
for (i in 1:length(sample_name)) {
    if (!grepl(".rds",tolower(sample_name[i]), fixed = TRUE)){
        print(c(sample_name[i],"is not R rds"))
    }
    }
}

## identify available cores for parallel processing
numCores <- detectCores()
cl <- makeCluster(numCores-1)
registerDoParallel(cl) 

## make a dataframe for expected doublet rates for each sample -- defined by the user in the parameters
doublet_rate_df <- data.frame(par_sample_names, par_expected_doublet_rate)
doublet_rate_df

## create empty list to be populated with existing Seurat objects
seu_list<-list()

##### if users want to remove doublets from downstream analyses
if (tolower(par_dropDN)=='yes') {
    print('The following are deleted')
    seu_list<-foreach (i_s=1:length(sample_name)) %do% {  
        seu<-readRDS(paste(output_dir,'/step3/objs3/',sample_name[i_s], sep=""))
        print(sample_name[i_s])
        sample_select <- sub(".rds.*", "", sample_name[i_s]) 
        sample_select
        seu[["Sample_ID"]] <- sample_select
        
        ## RunUMAP
        seu <- RunUMAP(seu, dims = 1:par_RunUMAP_dims, n.neighbors =par_RunUMAP_n.neighbors)

        ## parameter sweep 
        sweep.res.list_pbmc <- paramSweep_v3(seu, PCs = 1:par_PCs, sct = par_sct)
        sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
        bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
        
        ## select the pK that corresponds to max bcmvn to optimize doublet detection
        pK <- bcmvn_pbmc %>% filter(BCmetric == max(BCmetric)) %>%select(pK) 
        pK <- as.numeric(as.character(pK[[1]]))
       
        ## Homotypic Doublet Proportion Estimate 
        # annotations <- seu@meta.data$seurat_clusters
        # homotypic.prop <- modelHomotypic(annotations)           
        
        ## Seurat object-specific expected doublet rate
        # nExp_poi <- round(doublet_rate_df$par_expected_doublet_rate[doublet_rate_df$par_sample_names == sample_select]*nrow(seu@meta.data))  
        # nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
        
        ## run doubletFinder 
        seu <- doubletFinder_v3(seu, PCs = 1:par_PCs, pN = par_pN,pK = pK, nExp = par_rate_nExp,reuse.pANN = FALSE, sct = par_sct)
        
        ## set droplet classification 
        DF.classifications=colnames(seu[[]])[which(grepl('DF.classifications', colnames(seu[[]]), fixed=TRUE))]
        DF.classifications_u=eval(parse(text = paste('unique(seu$',DF.classifications,')', sep='')))        
        
        ## print UMAP with doublet/singlet classification
        DimPlot(seu, reduction = 'umap', group.by = DF.classifications)
        ggsave(paste(output_dir,'/step4/figs4/',sample_nameb[i_s],"_DF.classifications.pdf",sep=""))
        
        ##print doublet detection summary plot
        #meta df
        meta_df <- data.frame(seu@meta.data)
        #string detect column
        ## pANN
        df_pANN <- meta_df[,colnames(meta_df) %like% c("pANN")]
        ## classification
        df_class <- meta_df[,colnames(meta_df) %like% c("DF.class")]
        ## bind
        df_bind <- cbind(df_pANN, df_class)
        df_bind <- data.frame(df_bind)
        df_bind$df_pANN <- as.numeric(df_bind$df_pANN )

        ## plot pANN violin plot
        pANN <- ggplot(df_bind, aes(x = df_class, y = df_pANN, fill = df_class)) + 
        geom_violin() + 
        geom_jitter(size=0.001) +
        theme_classic() + 
        xlab("Identity") + 
        ylab("pANN") +
        theme(axis.text = element_text(size = 12))

        ## plot number of doublets and singlets
        count<- ggplot(df_bind, aes(x = df_class, fill = df_class)) + 
        geom_bar() +
        theme_classic() + 
        xlab("Identity") + 
        ylab("Number of droplets") +
        geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) + 
        theme(axis.text = element_text(size = 12))

        ##print plot
        ggarrange(pANN,count, ncol = 2, common.legend = T, legend = "bottom" )
        ggsave(paste(output_dir,'/step4/figs4/',sample_nameb[i_s],"_doublet_summary.pdf",sep=""))

        ## print txt file with number of doublets detected
        doublet_df <- nrow(subset(df_bind, df_class == "Doublet"))
        doublet_df <- data.frame(doublet_df)
        colnames(doublet_df) <- "Number of predicted doublets"
        write.csv(doublet_df, file= paste(output_dir,'/step4/info4/n_predicted_doublets_',sample_nameb[i_s],".txt", sep=""))

        ## remove doublets
        aa<-as.character(DF.classifications_u)[!as.character(DF.classifications_u) %in% "Doublet"]
        Idents(seu) <- DF.classifications
        seu=subset(seu,idents=aa)
        
        ## save individual Seurat objects as RDS
        saveRDS(seu, paste(output_dir,'/step4/objs4/',sample_nameb[i_s],'.rds', sep=""))
        write.csv(colnames(seu[[]]), file= paste(output_dir,'/step4/info4/meta_info_',sample_nameb[i_s],".txt", sep=""))
        
        ## save RNA expression matrix
        if (tolower(par_save_RNA)=='yes') {
        mat <- GetAssayData(object = seu, assay = "RNA", slot = "data")
        writeMM(mat,file= paste(output_dir,'/step4/info4/',sample_nameb[i_s],"_RNA.txt", sep=""))
        }

        ## save metadata dataframe
        if (tolower(par_save_metadata)=='yes') {
        write.csv(seu[[]], file = paste(output_dir,'/step4/info4/MetaData_',sample_nameb[i_s],'.txt', sep=""), quote = TRUE, sep = ",")
        }
    }
}


##### if users do not want to remove doublets from downstream analyses
if (tolower(par_dropDN)=='no') {
    print('The droplet are not deleted')
    seu_list<-foreach (i_s=1:length(sample_name)) %do% {  
        seu<-readRDS(paste(output_dir,'/step3/objs3/',sample_name[i_s], sep=""))
        print(sample_name[i_s])
        sample_select <- sub(".rds.*", "", sample_name[i_s]) 
        seu[["Sample_ID"]] <- sample_select
        
        ## RunUMAP
        seu <- RunUMAP(seu, dims = 1:par_RunUMAP_dims, n.neighbors =par_RunUMAP_n.neighbors)

        ## parameter sweep 
        sweep.res.list_pbmc <- paramSweep_v3(seu, PCs = 1:par_PCs, sct = par_sct)
        sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
        bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
        
        ## select the pK that corresponds to max bcmvn to optimize doublet detection
        pK <- bcmvn_pbmc %>% filter(BCmetric == max(BCmetric)) %>%select(pK) 
        pK <- as.numeric(as.character(pK[[1]]))
        
        ## Homotypic Doublet Proportion Estimate 
        annotations <- seu@meta.data$seurat_clusters
        homotypic.prop <- modelHomotypic(annotations)          
        
        ## Seurat object-specific expected doublet rate
        nExp_poi <- round(doublet_rate_df$par_expected_doublet_rate[doublet_rate_df$par_sample_names == sample_select]*nrow(seu@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
        nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
        
        ## run doubletFinder 
        seu <- doubletFinder_v3(seu, PCs = 1:par_PCs, pN = par_pN,pK = pK, nExp = nExp_poi.adj,reuse.pANN = FALSE, sct = par_sct)
        
        ## set droplet classification 
        DF.classifications=colnames(seu[[]])[which(grepl('DF.classifications', colnames(seu[[]]), fixed=TRUE))]
        DF.classifications_u=eval(parse(text = paste('unique(seu$',DF.classifications,')', sep='')))
        
        ## print UMAP with doublet/singlet classification       
        DimPlot(seu, reduction = 'umap', group.by = DF.classifications) 
        ggsave(paste(output_dir,'/step4/figs4/',sample_nameb[i_s],"_DF.classifications.pdf",sep=""))

        ##print doublet detection summary plot
        #meta df
        meta_df <- data.frame(seu@meta.data)
        #string detect column
        ## pANN
        df_pANN <- meta_df[,colnames(meta_df) %like% c("pANN")]
        ## classification
        df_class <- meta_df[,colnames(meta_df) %like% c("DF.class")]
        ## bind
        df_bind <- cbind(df_pANN, df_class)
        df_bind <- data.frame(df_bind)
        df_bind$df_pANN <- as.numeric(df_bind$df_pANN )

        ## plot pANN violin plot
        pANN <- ggplot(df_bind, aes(x = df_class, y = df_pANN, fill = df_class)) + 
        geom_violin() + 
        geom_jitter(size=0.001) +
        theme_classic() + 
        xlab("Identity") + 
        ylab("pANN") +
        theme(axis.text = element_text(size = 12))

        ## plot number of doublets and singlets
        count<- ggplot(df_bind, aes(x = df_class, fill = df_class)) + 
        geom_bar() +
        theme_classic() + 
        xlab("Identity") + 
        ylab("Number of droplets") +
        geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) + 
        theme(axis.text = element_text(size = 12))

        ##print plot
        ggarrange(pANN,count, ncol = 2, common.legend = T, legend = "bottom" )
        ggsave(paste(output_dir,'/step4/figs4/',sample_nameb[i_s],"_doublet_summary.pdf",sep=""))

        ## print txt file with number of doublets detected
        doublet_df <- nrow(subset(df_bind, df_class == "Doublet"))
        doublet_df <- data.frame(doublet_df)
        colnames(doublet_df) <- "Number of predicted doublets"
        write.csv(doublet_df, file= paste(output_dir,'/step4/info4/n_predicted_doublets_',sample_nameb[i_s],".txt", sep=""))
        
        ## do not remove doublets
        aa<-as.character(DF.classifications_u)
        Idents(seu) <- DF.classifications
        seu=subset(seu,idents=aa)

        ## save individual Seurat objects as RDS
        saveRDS(seu, paste(output_dir,'/step4/objs4/',sample_nameb[i_s],'.rds', sep=""))
        write.csv(colnames(seu[[]]), file= paste(output_dir,'/step4/info4/meta_info_',sample_nameb[i_s],".txt", sep=""))

        ## save RNA expression matrix
        if (tolower(par_save_RNA)=='yes') {
            mat <- GetAssayData(object = seu, assay = "RNA", slot = "data")
            #  write.csv(mat, paste(output_dir,'/step2/info2/',sample_name[i],"_RNA.csv", sep=""))
            writeMM(mat,file= paste(output_dir,'/step4/info4/',sample_nameb[i_s],"_RNA.txt", sep=""))
        }

        ## save metadata dataframe
        if (tolower(par_save_metadata)=='yes') {
            write.csv(seu[[]], file = paste(output_dir,'/step4/info4/MetaData_',sample_nameb[i_s],'.txt', sep=""), quote = TRUE)
        }
        }
}

## save session information
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step4/info4/sessionInfo.txt', sep=""))

if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}