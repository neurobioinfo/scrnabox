# Step 2: Create Seurat object and remove ambient RNA
In Step 2, the CellRanger outputs generated in Step 1 (expression matrix, features, and barcodes) are used to create a Seurat object for each sample. The ambient RNA quantity is estimated and there is an option to correct gene expression profiles for RNA contamination using SoupX ([Young et al. 2020](https://academic.oup.com/gigascience/article/9/12/giaa151/6049831)). Then, CellRanger (if not removing ambient RNA) or SoupX (if removing ambient RNA) feature-barcode expression matrices are transformed into Seurat objects. Quality control measures are then computed to inform filtering in Step 3, including:

- the number of unique transcripts (genes);
- the total number of transcripts;
- the percentage of mitochondrial-encoded transcripts; 
- the percentage of ribosome gene transcripts.

Normalization and scaling is then performed on the individual Seurat objects prior to cell-cycle scoring.
 - - - -
The following parameters are adjustable for Step 2 (`~/working_directory/job_info/parameters/step2_par.txt`):

|Parameter|Default|Description|
|:--|:--|:--|
|par_save_RNA| Yes| Whether or not to export an RNA expression matrix|
|par_save_metadata| Yes|Whether or not to export a metadata dataframe|
|par_ambient_RNA| Yes|Whether or not to correct the feature-barcode expression matrices for ambient RNA contamination|
|par_min.cells_L| 3|Only retain genes expressed in a minimum number of cells|
|par_normalization.method|LogNormalize|Method to use for normalization|
|par_scale.factor|10000|Scale factor for scaling the data|
|par_selection.method|vst|Method for choosing the top variable features|
|par_nfeatures|2500|Number of features to select as top variable features|

 - - - -

To run Step 2, use the following command:
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 2
```

The resulting output files are deposited into `~/working_directory/step2`. For a description of the outputs see [here](outputs.md).

