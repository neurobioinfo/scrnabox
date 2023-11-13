# Step 7: Cluster annotation
In Step 7, cluster annotation is performed to define the cell types comprising the clusters identified in Step 6. ScRNAbox provides three tools to identify cell types comprising the clusters:

[Tool 1: Cluster marker gene identification and gene set enrichment analysis (GSEA)](#tool-1-cluster-marker-gsea)<br /> 
Seurat's _FindAllMarkers_ function is used to identify differentially expressed marker genes (DEG) by the Wilcoxon rank-sum test ([Macosko et al. 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4481139/)). DEGs in the positive direction     (Log2 fold-change > 0.00) are then tested for enrichment across user-defined gene set libraries that define cell types using the EnrichR tool ([Chen et al. 2013](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-128)). <br />

[Tool 2: Expression profiling of cell type markers and module scores](#tool-2-expression-profiling-of-known-marker-genes)   
Users can visualize the expression of individual genes and the aggregated expression of multiple genes. For each gene in a user-defined list, plots are produced to visualize its expression at the cluster or cell level. The aggregated expression of genes in a user-defined list are calculated using the Seurat AddModuleScore function ([Tirosh et al. 2016](https://pubmed.ncbi.nlm.nih.gov/27124452/). <br />

[Tool 3: Cell type predictions based on reference data](#tool-3-reference-based-annotation) 
Seurat's _FindTransferAnchors_ and _TransferData_ functions are used to leverage cell-type annotations from a reference Seurat object and generate annotation predictions for the query dataset ([Butlet et al. 2019](https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867419305598%3Fshowall%3Dtrue)). 
<br />

Additionally, users can [add cluster annotations](#adding-annotations) to the Seurat object.
 - - - -

The following parameters are adjustable for Step 7 (`~/working_directory/job_info/parameters/step7_par.txt`):

|Annotation tool|Parameter|Default|Description|
|:--|:--|:--|:--|
|**General**|par_save_RNA| Yes| Whether or not to export an RNA expression matrix|
|**General**|par_save_metadata| Yes|Whether or not to export a metadata dataframe|
|**General**|par_seurat_object| NULL |If users already have a Seurat object, they may provide the path to the Seurat object to initiate the pipeline at Step 7|
|**General**|par_level_cluster| integrated_snn_res.0.75| The cluster resolution that you want to annotate. If you skipped integration in Step 5, use par_level_cluster='RNA_snn_res.0.7', if you want to proceed with a clustering resolution of 0.7.|
|**Tool 1**|par_run_find_marker|Yes|Whether or not to find marker genes for each cluster|
|**Tool 1**|par_run_enrichR|No|Whether or not to run gene set enrichment analysis (GSEA) on the marker genes for each cluster using the EnrichR tools. Note that the HPC must have access to the internet to run GSEA.|
|**Tool 1**|par_top_sel|5|Number of top markers to identify based on avg_log2FC|
|**Tool 1**|par_db|Descartes_Cell_Types_and_Tissue_2021,<br /> CellMarker_Augmented_2021,<br />Azimuth_Cell_Types_2021|Character vector of EnrichR databases that define cell types. The top marker genes for each cluster will be tested for enrichment across these databases.|
|**Tool 2**|par_run_module_score|Yes|Whether or not to compute module score for aggregated expression |
|**Tool 2**|par_run_visualize_markers|Yes|Whether or not to visualize the expression of individual genes|
|**Tool 2**|par_module_score|NULL|Path to the csv file containing the gene sets for the module score|
|**Tool 2**|par_select_features_list|NULL|List of genes whose expression will be visualized individually|
|**Tool 2**|par_select_features_csv|NULL|If you want to define multiple lists of features to visualize individually, you can do so with a csv file. The header should contain the list names and all features belonging to the same list should be in the same column.|
|**Tool 3**|par_reference|NULL| Path defining the location of the reference Seurat object|
|**Tool 3**|par_reference_name|Reference| An arbitrary name for the reference object. This will be used to name the metadata slot.|
|**Tool 3**|par_level_celltype|NULL|The name of the metadata column in the reference Seurat object that defines cell types|
|**Tool 3**|par_FindTransferAnchors_dim|50| Number of dimensions from linear dimensional reduction used to find transfer anchors between the reference and query Seurat objects|
|**Tool 3**|par_futureglobalsmaxSize|60000 * 1024^2|This will increase your RAM usage so set this number mindfully|
|**Annotate**|par_annotate_resolution|integrated_snn_res.0.75| Which clustering resolution you want to annotate|
|**Annotate**|par_name_metadata|Celltypes1| The name of the metadata slot that will contain the annotations|
|**Annotate**|par_annotate_labels|NULL| A list of cluster labels. There must as many labels as clusters at the defined clustering resolution. Please refrain from using "_" when annotating.|

 - - - -

## Tool 1: Cluster marker GSEA
To run cluster marker GSEA, use the following command:
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--markergsea T
```
The resulting output files are deposited into `~/working_directory/step7/*/marker`. For a description of the outputs see [here](outputs.md).
 - - - -
 
**Note:**  In order to test the cluster marker genes for enrichment across EnrichR libraries, the HPC must have access to the internet. If your HPC cannot access the internet, you must set `par_run_enrichR = "no"`. The pipeline will still run differential gene expression and find the markers for each cluster. You can then take the pipeline output and run the enrichment step on your local machine directly in R. To do so, begin by downloading the `ClusterMarkers.csv` file obtained from running the above command with `par_run_find_marker = "yes"` to your computer:

```
scp username@beluga.computecanada.ca:~/working_directory/step7/info7/marker/ClusterMarkers.csv ~/Desktop/working_directory
```
Then run the following code in R:
```
#set up the environment
library(Seurat)
library(dplyr)
library(ggplot2)
library(enrichR)


## set up the parameters
PWD <-'/path/directory/where/outputs/will/be/deposited'
cluster_marker <- '/path/to/ClusterMarkers.csv'
db <- c('Descartes_Cell_Types_and_Tissue_2021','CellMarker_Augmented_2021','Azimuth_Cell_Types_2021')

## set up the function
annotation<-function(PWD,cluster_marker,db) {
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(enrichR)
  setwd(PWD)

  cluster_marker <- read.delim(cluster_marker, header = T, sep = ",") 
  
  dir.create("annot_enrich") 
  for (i in unique(cluster_marker$cluster)) {
    dir.create(paste0(PWD,"/annot_enrich", "/clust",i))
  }

  for (i in unique(cluster_marker$cluster)) {
    for (j in 1:length(db)) { 
    setwd(PWD)
    setwd(paste0('./annot_enrich/clust',i))
    N1.c0 <- cluster_marker %>% filter(cluster == i & avg_log2FC > 0)
    genes <- N1.c0$gene
    N1.c0.Er <- enrichr(genes, databases = db[j])
    if(is.null(N1.c0.Er)) next
    plotEnrich(N1.c0.Er[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value") +
      ggtitle(paste0(db[j], " cluster ", i))
    ggsave(file = paste0("plotenrich_clust_",i, "_", j,".pdf"))
    N1.Er.genes.1 <- N1.c0.Er[[1]] %>% dplyr::select(Term, Genes, Combined.Score)
    write.csv(N1.Er.genes.1,paste0("Er_genes_clust_",i,"_",db[j],".csv"))
  }
}
}

## perform enrichment
annotation(PWD,cluster_marker,db)
```
**Note**: Users can define whichever libraries they want. For more information regarding the available libraries, see [here](https://maayanlab.cloud/Enrichr/).
 - - - -

## Tool 2: Expression profiling of known marker genes
To profile the expression known marker genes, use the following command:
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--knownmarkers T
```
The resulting output files are deposited into `~/working_directory/step7/*/visualize_features` for individual expression or `~/working_directory/step7/*/module_score` for aggregated expression. For a description of the outputs see [here](outputs.md).

 - - - -
 

## Tool 3: Reference-based annotation
To perform reference-based annotation, use the following command:
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--referenceannotation T
```
The resulting output files are deposited into `~/working_directory/step7/*/reference_based_annotation`. For a description of the outputs see [here](outputs.md).
 - - - -

## Adding annotations
To add cluster annotations to the Seurat object's metadata, use the following command:
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--annotate T
```
The resulting output files are deposited into `~/working_directory/step7/*/annotate`. For a description of the outputs see [here](outputs.md).

**Note:** The cluster annotations from each iteration of the step will be retained, allowing users to define multiple clustering levels.

