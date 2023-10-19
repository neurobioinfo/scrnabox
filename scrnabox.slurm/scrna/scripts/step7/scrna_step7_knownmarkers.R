#!/usr/bin/env Rscript

###############################################################################
# step7 -- Expression of known marker genes
###############################################################################

## set sample ID metadata column -- this is standard and does not require parameter modification
par_level_genotype <- "Sample_ID"

## load parameters
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
pipeline_home=args[3]

## load library
.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr', 'xlsx', 'Matrix')
lapply(packages, library, character.only = TRUE)

## load parameters
source(paste(output_dir,'/job_info/parameters/step7_par.txt',sep=""))

###############################################################################
# step7 -- module score
###############################################################################

if (tolower(par_run_module_score)=='yes') {

## load existing Seurat objects
sample_name<-list.files(path = paste(output_dir, "/step6/objs6",sep=""),pattern = "*.rds")
seu_int<-readRDS(paste(output_dir,'/step6/objs6/',sample_name, sep=''))

## set cell identity to the clustering resolution defined by the user
Idents(seu_int) <- par_level_cluster

## create directories for module score
## figures 
OUT_DIR_figs <- paste(output_dir,"/step7/figs7",sep='') 
OUT_dir_figs_module_score <- paste(OUT_DIR_figs,"/module_score/",sep='') 
dir.create(OUT_dir_figs_module_score)
## info
OUT_DIR_info <- paste(output_dir,"/step7/info7",sep='') 
OUT_dir_info_module_score <- paste(OUT_DIR_info,"/module_score/",sep='') 
dir.create(OUT_dir_info_module_score)

## set output directory
PWD=OUT_dir_info_module_score
setwd(PWD)

## compute module score
source(paste(pipeline_home,'/general_codes/module_score.R',sep=''))

## save session info
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step7/info7/sessionInfo_module_score.txt', sep=""))
if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}
}

###############################################################################
# step7 -- visualize marker
###############################################################################

if (tolower(par_run_visualize_markers)=='yes') {

## load existing Seurat objects
sample_name<-list.files(path = paste(output_dir, "/step6/objs6",sep=""),pattern = "*.rds")
seu_int<-readRDS(paste(output_dir,'/step6/objs6/',sample_name, sep=''))

## set cell identity to the clustering resolution defined by the user
Idents(seu_int) <- par_level_cluster

## create directories for visualize features annotation method
## figures 
OUT_DIR_figs <- paste(output_dir,"/step7/figs7",sep='') 
OUT_dir_figs_visualize_features <- paste(OUT_DIR_figs,"/visualize_features/",sep='') 
dir.create(OUT_dir_figs_visualize_features)

## set output directory
PWD=OUT_dir_figs_visualize_features
setwd(PWD)

###################
## visulize select features with a list in the parameters
###################
if (exists("par_select_features_list")) {
#define the features
select_features <- par_select_features_list

#violin plot
vln_plt <- VlnPlot(seu_int, features = select_features, group.by = par_level_cluster, pt.size = 0)
ggsave(file = paste(OUT_dir_figs_visualize_features,'list_violin_plot.pdf', sep=''), dpi = 300, height = 15, width = 15, unit = 'in')

#feature plot
Feature_plt <- FeaturePlot(seu_int, features = select_features, raster = FALSE)
ggsave(file = paste(OUT_dir_figs_visualize_features ,'list_feature_plot.pdf', sep=''), dpi = 300, height = 15, width = 20, unit = 'in')

#dotplot
dot_plt <- DotPlot(seu_int, features = select_features, group.by = par_level_cluster) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(file = paste(OUT_dir_figs_visualize_features ,'list_dot_plot.pdf', sep=''))
}

###################
## visulize select features with a csv
###################

if (exists("par_select_features_csv")) {
gene_sets <- read.delim(par_select_features_csv, header = T, sep = ",", na.strings=c("","NA")) 

#convert user inputed dataframe to named list
gene_lists <- list()                   # Create empty list
for(i in 1:ncol(gene_sets)) {             # Using for-loop to add columns to list
  gene_lists[[i]] <- gene_sets[ , i]
}
names(gene_lists) <- colnames(gene_sets)#set names of list

## remove NA from list
gene_lists <- lapply(gene_lists, function(x) x[!is.na(x)])

for (i in 1:length(gene_lists)) {
    # Add module scores for each gene list
    #try(
            #violin plot
            vln_plt <- VlnPlot(seu_int, features = gene_lists[[i]], group.by = par_level_cluster, pt.size = 0)
            ggsave(file = paste(OUT_dir_figs_visualize_features, names(gene_lists[i]),'_violin_plot.pdf', sep=''), dpi = 300, height = 15, width = 15, unit = 'in')

            #feature plot
            Feature_plt <- FeaturePlot(seu_int, features = gene_lists[[i]], raster = FALSE)
            ggsave(file = paste(OUT_dir_figs_visualize_features, names(gene_lists[i]) ,'_feature_plot.pdf', sep=''), dpi = 300, height = 15, width = 20, unit = 'in')

            #dotplot
            dot_plt <- DotPlot(seu_int, features = gene_lists[[i]], group.by = par_level_cluster) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
            ggsave(file = paste(OUT_dir_figs_visualize_features, names(gene_lists[i]),'_dot_plot.pdf', sep=''))
       #)
    }
}
 

## save session info
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step7/info7/sessionInfo_viualize_features.txt', sep=""))
if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
}
}
