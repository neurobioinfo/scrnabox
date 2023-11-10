# Step 4: Demultiplexing and doublet detection (HTO track)
In Step 4 of the HTO track, Seurat’s implementation (_MULTIseqDemux_) of the tag assignment algorithm outlined in Multi-seq is used to demultiplex pooled samples and identify doublets according to the expression matrices of the sample-specific barcodes ([McGinnis et al 2019](https://pubmed.ncbi.nlm.nih.gov/31209384/)).

 - - - -

The following parameters are adjustable for Step 4 (`~/working_directory/job_info/parameters/step4_par.txt`):

|Parameter|Default|Description|
|:--|:--|:--|
|par_save_RNA| Yes| Whether or not to export an RNA expression matrix|
|par_save_metadata| Yes|Whether or not to export a metadata dataframe|
|par_seurat_object| NULL |If users already have a Seurat object(s), they may provide the path to a directory that contains an existing Seurat object(s) to initiate the pipeline at Step 4|
|par_normalization.method|CLR|Method for normalizing the HTO assay|
|par_scale.factor| 1000|Scale factor for scaling the HTO assay|
|par_selection.method|vst|Method for selecting the most variable features in the HTO assay|
|par_nfeatures|5|Number of features to select as top variable features for the HTO assay. This value is dependent on the number of sample specific barcodes used in the experiment|
|par_dims_umap|5|Number of dimensions to use as input features for uniform manifold approximation and projection (UMAP) of HTO assay|
|par_n.neighbor|65|Number of neighboring points to use in local approximations of manifold structure|
|par_dimensionality_reduction|Yes|Whether or not to perform linear dimensionality reduction on the HTO assay|
|par_npcs_pca|30|Total Number of principal components to compute and store for principal component analysis (PCA) of HTO assay|
|par_dropDN|Yes|Whether or not to remove predicted doublets and negatives from downstream analyses|
|par_label_dropDN|Doublet, Negative| Labels used to identify doublet and negative droplets|
|par_quantile|0.9|The quantile to use for droplet classification using _MULTIseqDemux_|
|par_autoThresh|TRUE| Whether or not to perform automated threshold finding to define the best quantile for droplet classification using _MULTIseqDemux_|
|par_maxiter|5|Maximum number of iterations to use if autoThresh = TRUE|
|par_RidgePlot_ncol|3|Number of columns used to display RidgePlots, which visualizes the enrichment of barcode labels across samples|
|par_old_antibody_label|NULL| If you wish to rename the barcode labels, first list the existing barcode labels in this parameter. old antibody labels can be identified in the "_old_antibody_label_MULTIseqDemuxHTOcounts" file produced by running Step 4 msd|
|par_new_antibody_label|NULL|If you wish to rename the barcode labels, list the new labels corresponding to the old labels listed in the parameter above|

 - - - -

To demultiplex the samples and identify doublets, the first step is to obtain the barcode labels used in the analysis by running the following command:

```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 4 \
--msd T 
```
**Note:** This step will produce the `old_antibody_label_MULTIseqDemuxHTOcounts.csv` file, which contains the names of the old HTO labels. The names of the HTO labels can be revised to be more descriptive in the execution parameters of this step (par_old_antibody_label; par_new_antibody_label)

Next, demultiplex the samples and identify doublets by running the following command:

```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 4 
```
The resulting output files are deposited into `~/working_directory/step4`. For a description of the outputs see [here](outputs.md).

