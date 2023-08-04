#Rand index
library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(vctrs)
library(fossil)
library(openxlsx)
library(stringr)
library(ggpubr)


#read in seurat object
# seu_int <- readRDS("/Users/mfiorini/Desktop/scRNA pipeline paper/seu_step6.rds") #we will take step 6 seurat object

Seurat::DefaultAssay(seu_int) <- "integrated"
# seu_int <- ScaleData(seu_int, verbose = FALSE)
# seu_int <- RunPCA(seu_int, npcs = 25, verbose = FALSE)
# seu_int <- RunUMAP(seu_int, dims = 1:25, n.neighbors =30)
# seu_int <- FindNeighbors(seu_int,  dims = 1:25, k.param = 30, prune.SNN = 1/15)
# seu_int <- Seurat::FindClusters(seu_int, resolution = par_FindClusters_resolution)

#### THE CODE ABOVE ALREADY EXISTS IN THE STEP 6 CODE


#### THE CODE BELOW IS NOT YET IMPLEMNETED, AND COULD COME after the above code. The code does not change the existing Seurat object, it only produces a plot and a summary
# par_FindClusters_resolution = seq(0.1, 0.3, by=0.1) # this is define in the parameters of step 6
OUT_DIR <- paste(output_dir,"/step6/ARI",sep='')  #this is defined in the parameters of step 6
dir.create(OUT_DIR)



#function to run louvain clustering multiple times and calculate ARI
calculate_rand_index <- function(seu_object, par_FindClusters_resolution, OUT_DIR) {
  mybiglist <- list()
  means <- list()
  standard_devs <- list()
  result <- list()
  reps = par_RI_reps
  for(i in par_FindClusters_resolution) { 
    for(n in 1:reps) { 
    seu_object <- Seurat::FindClusters(seu_int, resolution = i, random.seed = n) 
    #retrieve the clustering
    column <- paste0("integrated_snn_res.",i)
    df <- data.frame(seu_object@meta.data)
    df1 <- df %>% select(column)
    result[[n]] <- df1[,1]
    }
    # result 
    # print(result)
    mybiglist <- list()
    for (f in 1:length(result)){
      # inner loop
      for (j in 1:length(result)){
        if (f==j) {
          name <- paste0(f, ",", j)
          tmp <- -5  #place holder to remove contrasts between the same reps
          mybiglist[[name]] <- tmp
          mybiglist
        } else {
          name <- paste0(f, ",", j)
          tmp <- adj.rand.index(result[[f]], result[[j]])
          mybiglist[[name]] <- tmp
          mybiglist
        }
      }
    }
    #remove values from list if they are equal to -5 (i.e those computing the adj rand index between the same clusters)
    mybiglist_test <-  mybiglist[mybiglist != -5]
    # print(mybiglist_test)
    mean <- mean(unlist((mybiglist_test)))
    sd <- sd(unlist((mybiglist_test)))
    means[[paste0("mean_",i)]] <- mean                           
    standard_devs[[paste0("sd_",i)]] <- sd                            
    # print(means)
    # print(standard_devs)
  }
  # standard_devs
  # #plot standard dev
  sd_df<- data.frame(standard_devs)
  sd_df<- data.frame(t(sd_df))
  sd_df$resolution <- rownames(sd_df)
  colnames(sd_df) <- c("sd", "Resolution")
  sd_df$group = "group"
  sd_plot <- ggplot(sd_df, aes(x = Resolution, y = sd, group = group)) +
    geom_line() + geom_point() + theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank())
  # means
  #plot means
  sd_mean<- data.frame(means)
  sd_mean<- data.frame(t(sd_mean))
  sd_mean$resolution <- rownames(sd_mean)
  colnames(sd_mean) <- c("mean", "Resolution")
  sd_mean$group = "group"
  mean_plot <- ggplot(sd_mean, aes(x = Resolution, y = mean, group = group)) +
    geom_line() + geom_point() + theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank())
  #arrange plots
  ggarrange(mean_plot, sd_plot, ncol = 1, nrow = 2)
  #combine sd and means list
  list3 <- c(means,standard_devs)
  list3_df <- data.frame(list3)
  list3_df
  #write excel file
  write.xlsx(list3_df, paste(OUT_DIR, "/test.xlsx", sep='')) #output dir must be changed to the output dir of the pipeline
  plots <- ggarrange(mean_plot, sd_plot, ncol = 1, nrow = 2, align = "hv")
  plots
}  

test_plot <- calculate_rand_index(seu_int,par_FindClusters_resolution, OUT_DIR ) #this must be run inorder to prepare for below code



#function to plot ARI and cluster together
calculate_clusters<- function(seu_int, par_FindClusters_resolution, test_plot, OUT_DIR) {
  mybiglist <- list()
  for(i in par_FindClusters_resolution) { 
    column <- paste0("integrated_snn_res.",i)
    df <- data.frame(seu_int@meta.data)
    df1 <- df %>% select(column)
    n_levels <-length(unique(df1[,1]))
    mybiglist[[column]] <- n_levels
  }
  mybiglist
  sd_df<- data.frame(mybiglist)
  sd_df<- data.frame(t(sd_df))
  sd_df$resolution <- rownames(sd_df)
  colnames(sd_df) <- c("cluster", "Resolution")
  sd_df$group = "group"
  cluster_plot <- ggplot(sd_df, aes(x = Resolution, y = cluster, group = group)) +
    geom_line() +
    geom_point() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggarrange(test_plot, cluster_plot, ncol = 1, nrow = 2, heights = c(1.25,1))
  ggsave(file = paste(OUT_DIR, "/ARI.png", sep='')) #change this to piepline directory
}  

calculate_clusters(seu_int, par_FindClusters_resolution, test_plot, OUT_DIR) #this must be run

