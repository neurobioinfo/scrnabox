# Step 5: Creation of a single Seurat object from all samples 
In Step 5, individual Seurat objects from each sample are combined to enable the joint analysis across samples. Users can either merge or integrate their Seurat objects ([Stuart et al. 2019](https://pubmed.ncbi.nlm.nih.gov/31178118/)). Alternatively, if experiments are limited to a single sequencing run, merging/integration can be bypassed; however, Step 5 must still be run because normalization, scaling, and linear dimensional reduction is then performed to inform the optimal parameters for clustering in Step 6.

**Note:** For more information regarding the difference between merging and integration, please see our pre-print manuscript [here]().
 - - - -

The following parameters are adjustable for Step 5 (`~/working_directory/job_info/parameters/step5_par.txt`):

|Parameter|Default|Description|
|:--|:--|:--|
|par_save_RNA| No| Whether or not to export an RNA expression matrix|
|par_save_metadata| No|Whether or not to export a metadata dataframe|
|par_seurat_object| NULL |If users already have a Seurat object(s), they may provide the path to a directory that contains an existing Seurat object(s) to initiate the pipeline at Step 5|
|par_one_seurat| No| Whether or not the experiment comprises of only one sequencing run. If this parameter is set to "Yes", set par_integrate_seurat and par_merge_seurat to "No".|
|par_integrate_seurat| Yes| Whether or not to integrate the samples. If "Yes", par_merge_seurat must be "No". |
|par_merge_seurat| No| Whether or not to merge the samples. If "Yes", par_integrate_seurat must be "No". |
|par_DefaultAssay|RNA|The assay to perform normalization, scaling, and linear dimensiona reduction on. For most use cases this will be RNA.|
|par_normalization.method|LogNormalize|Method to use for normalization|
|par_scale.factor|10000|Scale factor for scaling the data|
|par_selection.method|vst|Method for detecting top variable features|
|par_nfeatures|2500|Number of features to select as top variable features|
|par_FindIntegrationAnchors_dim|25|Which dimensions to use from the canonical correlation analysis (CCA) to specify the neighbor search space|
|par_RunPCA_npcs|30| Total Number of principal components to compute and store for principal component analysis (PCA)|
|par_RunUMAP_dims|10| Number of dimensions to use as input features for uniform manifold approximation and projection (UMAP)|
|par_RunUMAP_n.neighbors|65|Number of neighboring points used in local approximations of manifold structure|
|par_compute_jackstraw |No|Whether or not to perform JackStraw computation. This computation takes a long time.|

 - - - -

To run Step 5, use the following command:
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 5 
```
The resulting output files are deposited into `~/working_directory/step5`. For a description of the outputs see [here](outputs.md).

