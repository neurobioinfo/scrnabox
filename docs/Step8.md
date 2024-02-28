# Step 8: Differential gene expression (DGE) analysis
In Step 8, DGE analysis is computed to identify differentially expressed genes (DEG) between two conditions. Prior to computing DGE, users can [add metdata](#add-metadata) containing phenotypic and experimental data to the Seurat object, which can then be used to define the groups used for DGE analysis. In order to define the contrasts used in the DGE analysis, users must modify the [contrast matrices ](#contrast-matrices) prior to submitting the command to [compute DGE](#computing-dge). ScRNAbox can compute DGE between conditions using all cell types or cell type groups. Furthermore, scRNAbox provides two frameworks for computing DGE: <br />

**1) Cell-based DGE**<br />
Cells are used as replicates and DGE is computed using the Seurat _FindMarkers_ ([Macosko et al. 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4481139/)). While _FindMarkers_ supports several statistical frameworks to compute DGE, we set the default method in our implementation to MAST, which is tailored for scRNAseq data ([Finak et al. 2015](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0844-5))<br />

**2) Sample-based DGE**<br />
Samples are used as replicates by applying a pseudo-bulk analysis. The Seurat _AggregateExpression_ function is used to compute the sum of RNA counts for each gene across all cells from a particular sample ([Cao et al. 2022](https://academic.oup.com/nar/article/50/21/e121/6709246)). The DESq2 statistical framework is then used to compute DGE between conditions using the aggregated counts. ([Love et al. 2014](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8))

 - - - -
The following parameters are adjustable for Step 8:

|DGE method|Parameter|Default|Description|
|:--|:--|:--|:--|
|**General**|par_save_RNA| Yes| Whether or not to export an RNA expression matrix|
|**General**|par_save_metadata| Yes|Whether or not to export a metadata dataframe|
|**General**|par_seurat_object| NULL |If users already have a Seurat object, they may provide the path to the Seurat object to initiate the pipeline at Step 7|
|**Add metadata**|par_merge_meta|orig.ident|The column from the Seurat metdata that will be used to merge the new metadata. This column must also exist in the submitted csv file contaning new metadata.|
|**Add metadata**|par_metadata|NULL|csv file containing metadata to be added to the Seurat object|
|**Cell-based DGE with all cells**|par_run_cell_based_all_cells|Yes|Whether or not to compute cell-based DGE with all cells |
|**Cell-based DGE with cell type groups**|par_run_cell_based_cell_type_groups|Yes|Whether or not to compute cell-based DGE with cell type groups|
|**Sample-based DGE with all cells**|par_run_sample_based_all_cells|Yes|Whether or not to compute sample-based DGE with all cells|
|**Sample-based DGE with cell type groups**|par_run_sample_based_cell_type_groups|Yes|Whether or not to compute sample-based DGE with cell type groups|
|**Cell-based DGE**|par_statistical_method|MAST| Which statistical framework to use for computing cell-based DGE|

 - - - -

## Add metadata
To add metadata to the Seurat object,  use the following command:
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--addmeta T
```
An example of a metadata csv file is available [here](https://github.com/neurobioinfo/scrnabox/blob/main/tutorial/Midbrain_dataset_example_files/metadata.csv). <br />
The resulting output files are deposited into `~/working_directory/step8`. For a description of the outputs see [here](outputs.md).

 - - - -

## Contrast matrices
#### Cell-based DGE using all cells
To perform cell-based DGE using all cells, users must fill in the `step8_contrast_cell_based_all_cells.txt` file located in `~/working_directory/job_info/parameters`. The contrast matrix contains the following columns:

1. **contast_name:** An informative name for the contrast. This will appear as the name of the output spreadsheet. 
2. **meta_data_variable:** The metadata slot containing the Sample IDs defined in group1 and group2 
3. **group1:** A list of sample IDs to be contrasted against the sample IDs listed in group2
4. **group2:**A list of sample IDs to be contrasted against the sample IDs listed in group1

Multiple contrasts can be defined in the same file. In addition, multiple samples can be listed under group1 and group 2. For example: 
```
contrast_name meta_data_variable group1 group2
Design1 orig.ident Control1,Control2,Control3 Case1,Case2,Case3
Design3 DiseaseStatus HealthyControl Disease
```
 - - - -

#### Cell-based DGE using cell type groups
To perform cell-based DGE using cell type groups, users must fill in the `step8_contrast_cell_based_celltype_groups.txt` file located in `~/working_directory/job_info/parameters`. The contrast matrix contains the following columns:

1. **contast_name:** An informative name for the contrast. This will appear as the name of the output spreadsheet. 
2. **meta_data_celltype:** The metadata slot containing cell type annotations
3. **cell_type:** The cell type used to compute DGE 
2. **meta_data_variable:** The metadata slot containing the Sample IDs defined in group1 and group2 
3. **group1:** A list of sample IDs to be contrasted against the sample IDs listed in group2
4. **group2:**A list of sample IDs to be contrasted against the sample IDs listed in group1

Multiple contrasts can be defined in the same file. In addition, multiple samples can be listed under group1 and group 2. For example: 
```
contrast_name meta_data_celltype cell_type meta_data_variable group1 group2
Design1 Annotation1 Neuron orig.ident Control1,Control2,Control3, Case1,Case2,Case3,
Design2 Annotation2 Microglia DiseaseStatus HealthyControl Disease
```
 - - - -

#### Sample-based DGE using all cells
To perform sample-based DGE using all cells, users must fill in the `step8_contrast_sample_based_all_cells.txt` file located in `~/working_directory/job_info/parameters`. The contrast matrix contains the following columns:

1. **ContrastName:** An informative name for the contrast. This will appear as the name of the output spreadsheet. 
2. **MainContrast:** The metadata slot containing the two groups used for the main contrast (e.g. case and control)
3. **Sample_ID:** The metadata slot containing the Sample IDs of the individual subjects (e.g. sample 1, sample 2, etc.)

```
ContrastName MainContrast SampleID
Design DiseaseStatus orig.ident
```

In addition, users may add additional columns if they want to further group their samples. For example, users may wich to group samples by experimental batch:

```
ContrastName MainContrast SampleID Batch
Design DiseaseStatus orig.ident Batch_Id
```
In this case, **Batch** is arbitrary, but **Batch_ID** must be a metadata slot. 
 - - - -

#### Sample-based DGE using cell type groups
To perform sample-based DGE using all cells, users must fill in the `step8_contrast_sample_based_celltype_groups.txt` file located in `~/working_directory/job_info/parameters`. The contrast matrix contains the following columns:

1. **ContrastName:** An informative name for the contrast. This will appear as the name of the output spreadsheet. 
2. **CellType:** The metadata slot containing cell type annotations
3. **MainContrast:** The metadata slot containing the two groups used for the main contrast (e.g. case and control)
4. **Sample_ID:** The metadata slot containing the Sample IDs of the individual subjects (e.g. sample 1, sample 2, etc.)

```
ContrastName CellType MainContrast SampleID
Design Annotation1 DiseaseStatus orig.ident
```

In addition, users may add additional columns if they want to further group their samples. For example, users may wich to group samples by experimental batch:

```
ContrastName CellType MainContrast SampleID Batch
Design Annotation1 DiseaseStatus orig.ident Batch_ID
```
In this case, **Batch** is arbitrary, but **Batch_ID** must be a metadata slot. 
 - - - -

## Compute DGE

To compute DGE, use the following command:
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--rundge T
```

The resulting output files are deposited into `~/working_directory/step8`. For a description of the outputs see [here](outputs.md).
