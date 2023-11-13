# Step 3: Quality control and generation of filtered data objects
In Step 3, low quality cells are filtered based on the user-defined thresholds for:

- the number of unique transcripts (genes; nFeaturesRNA);
- the total number of transcripts (nCountsRNA);
- the percentage of mitochondrial-encoded transcripts; 
- the percentage of ribosome gene transcripts.

In addition, users can  remove or regress a custom gene list from the dataset. Finally, normalization and scaling is performed on the filtered Seurat objects. <br />

 - - - -

The following parameters are adjustable for Step 3 (`~/working_directory/job_info/parameters/step3_par.txt`):

|Parameter|Default|Description|
|:--|:--|:--|
|par_save_RNA| Yes| Whether or not to export an RNA expression matrix|
|par_save_metadata| Yes|Whether or not to export a metadata dataframe|
|par_seurat_object| NULL |If users already have a Seurat object(s), they may provide the path to a directory that contains an existing Seurat object(s) to initiate the pipeline at Step 3|
|par_nFeature_RNA_L|300 |Only retain cells expressing a minimum number of unique RNA transcripts|
|par_nFeature_RNA_U|10000 |Only retain cells expressing a maximum number of unique RNA transcripts|
|par_nCount_RNA_L|300 |Only retain cells with a minimum number of total RNA transcripts|
|par_nCount_RNA_U|20000 |Only retain cells with a maximum number of total RNA transcripts|
|par_mitochondria_percent_L|0 | Only retain cells with a minimum percentage of mitochondrial-encoded genes|
|par_mitochondria_percent_U|20 |Only retain cells with a maximum percentage of mitochondrial-encoded genes|
|par_ribosomal_percent_L|0 |Only retain cells with a minimum percentage of ribosome genes|
|par_ribosomal_percent_U|100 |Only retain cells with a maximum percentage of ribosome genes|
|par_remove_mitochondrial_genes|No| Whether or not to remove mitochondrial genes|
|par_remove_ribosomal_genes|No| Whether or not to remove ribosomal genes|
|par_remove_genes|NULL|If users want to remove specific genes from their data, they may define a list of gene identifiers|
|par_regress_cell_cycle_genes|No|Whether or not to regress cell cycle genes|
|par_regress_custom_genes|No|Whether or not to regress a custom list of genes|
|par_regress_genes|NULL|List of custom genes to regress|
|par_normalization.method|LogNormalize|Method to use for normalization|
|par_scale.factor|10000|Scale factor for scaling the data|
|par_selection.method|vst|Method for choosing the top variable features|
|par_nfeatures|2500|Number of features to select as top variable features|
|par_top|10|Number of most variable features to be reported in the csv file|
|par_npcs_pca|30|Total Number of principal components to compute and store for principal component analysis (PCA)|

 - - - -
To run Step 3, use the following command:
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 3
```

The resulting output files are deposited into `~/working_directory/step3`. For a description of the outputs see [here](outputs.md).

