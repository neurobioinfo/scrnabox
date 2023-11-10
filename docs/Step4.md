# Step 4: Doublet removal (standard track)
In Step 4 of the standard analysis track, doublets (barcodes produced by sequencing two or more cells) are identified and optionally removed from downstream analysis using the DoubletFinder tool ([McGinnis et al. 2019](https://www.cell.com/cell-systems/pdfExtended/S2405-4712(19)30073-0)). 
 - - - -

The following parameters are adjustable for Step 4 of the standard track (`~/working_directory/job_info/parameters/step4_par.txt`):

|Parameter|Default|Description|
|:--|:--|:--|
|par_save_RNA| Yes| Whether or not to export an RNA expression matrix|
|par_save_metadata| Yes|Whether or not to export a metadata dataframe|
|par_seurat_object| NULL |If users already have a Seurat object(s), they may provide the path to a directory that contains an existing Seurat object(s) to initiate the pipeline at Step 4|
|par_RunUMAP_dims|25| Number of dimensions to use as input features for uniform manifold approximation and projection (UMAP)|
|par_RunUMAP_n.neighbors|45|Number of neighboring points used in local approximations of manifold structure|
|par_dropDN| Yes| Whether or not to remove predicted doublets from downstream analyses|
|par_PCs|25| The number of statistically significant principal components. Can be informed by elbow plot produced in Step 3|
|par_pN|0.25| The number of artificial doublets to generate. DoubletFinderr is largely invariant to this parameter. We suggest keeping 0.25|
|par_sct|FALSE|Logical representing whether SCTransform was used during original Seurat object pre-processing|
|par_sample_names|NULL| A list of sample names for each sample in the experiement, corresponding to the expected doublet rates listed in the parameter below. Sample names should be the same as those used to produce the `samples_info` folder during the setup procedures.|
|par_expected_doublet_rate|NULL| A vector of expected doublet rates for each sample (e.g. for a 5% expected doublet rate, write 0.05). The expected doublet rates for each sample should be listed in the same order as the sample names in the above parameter. Make sure to have as many expected doublet rates listed as you have samples.|

**Note:** For more information regarding the expected doublet rates, please see the 10X Genomics [documentation](https://kb.10xgenomics.com/hc/en-us/articles/360059124751-Why-is-the-multiplet-rate-different-for-the-Next-GEM-Single-Cell-3-LT-v3-1-assay-compared-to-other-single-cell-applications-).

 - - - -

To run Step 4, use the following command:
```
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 4 
```
The resulting output files are deposited into `~/working_directory/step4`. For a description of the outputs see [here](outputs.md).

