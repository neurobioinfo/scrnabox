# Step 6: Clustering
In Step 6, clustering is performed to define groups of cells with similar expression profiles using the Seurat implementation of the Louvain network detection with PCA dimensionality reduction as input  ([Macosko et al. 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4481139/)).
 - - - -
The following parameters are adjustable for Step 6 (`~/working_directory/job_info/parameters/step6_par.txt`):

|Parameter|Default|Description|
|:--|:--|:--|
|par_save_RNA| No| Whether or not to export an RNA expression matrix|
|par_save_metadata| No|Whether or not to export a metadata dataframe|
|par_seurat_object| NULL |If users already have a Seurat object, they may provide the path to the Seurat object to initiate the pipeline at Step 6|
|par_skip_integration|No|Whether or not the user skipped integration in Step 5| 
|par_FindNeighbors_dims|30| Number of dimensions from linear dimensional reduction used as input to identify neighbours. Can be informed by the elbow and Jackstraw plots produced in Step 5|
|par_FindNeighbors_k.param|60|Defines k for the k-nearest neighbor algorithm|
|par_FindNeighbors_prune.SNN|1/15|Sets the cutoff for acceptable Jaccard index when computing the neighborhood overlap for the shared nearest-neighbour (SNN) construction
|par_FindClusters_resolution|0, 0.05, 0.2, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0|Value of the clustering resolution parameter. You may provide multiple resolution values|
|par_compute_ARI|Yes| Whether or not you want to compute the Adjusted Rand Index (ARI) between clusters at a given clustering resolution|
|par_RI_reps|25|Number of iterations for clustering the data at a given resolution in order to calculate the ARI|
 
 - - - -
To run Step 6, use the following command:
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 6 
```
The resulting output files are deposited into `~/working_directory/step6`. For a description of the outputs see [here](outputs.md).
