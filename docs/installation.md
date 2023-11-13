# Installation
To use the scRNAbox pipeline, the folowing must be installed on your High-Performance Computing (HPC) system:

- [scrnabox.slurm](#scrnaboxslurm-installation)
- [CellRanger](#cellranger-installation)
- [R and R packages](#r-library-preparation-and-r-package-installation)

 - - - -

### scrnabox.slurm installation

`scrnabox.slurm` is written in bash and can be used with any Slurm system. To download the latest version of `scrnabox.slurm` (v0.1.35) run the following command: 
```
wget https://github.com/neurobioinfo/scrnabox/releases/download/v0.1.35/scrnabox.slurm.zip
unzip scrnabox.slurm.zip
```

For a description of the options for running `scrnabox.slurm` run the following command:
```
bash /pathway/to/scrnabox.slurm/launch_scrnabox.sh -h 
```

If the `scrnabox.slurm` has been installed properly, the above command should return the folllowing:
```
        mandatory arguments:
                -d  (--dir)  = Working directory (where all the outputs will be printed) (give full path)
                --steps  =  Specify what steps, e.g., 2 to run just step 2, 2-4, run steps 2 through 4)

        optional arguments:
                -h  (--help)  = See helps regarding the pipeline arguments. 
                --method  = Select your preferred method: HTO and SCRNA for hashtag, and Standard scRNA, respectively. 
                --msd  = You can get the hashtag labels by running the following code 
                --markergsea  = Identify marker genes for each cluster and run marker gene set enrichment analysis (GSEA) using EnrichR libraries. 
                --knownmarkers  = Run module score and visualize the expression of known cell type marker genes. 
                --referenceannotation  = Run module score and visualize the expression of known cell type marker genes. 
                --annotate  = Run module score and visualize the expression of known cell type marker genes. 
                --addmeta  = Add metadata columns to the Seurat object 
                --rundge  = Perform differential gene expression contrasts 
                --seulist  = You can directly call the list of seurat objects to the pipeline.  
 
 ------------------- 
 For a comprehensive help, visit https://github.com/neurobioinfo/scrnabox for documentation.
```
 - - - -

### CellRanger installation

For information regarding the installation of `CellRanger`, please visit the 10X Genomics [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation). If CellRanger is already installed on your HPC system, you may skip the CellRanger installation procedures.

 - - - -

### R library preparation and R package installation
R must be installed on the HPC system. If R is not already installed, the user should contact their system administrator. If R is installed on the HPC system, users must load `R` (v4.2 or later) into their working environment: 

```
# load R
module load r/4.2.1
```
Then, users must run the following installation code, which will create a directory where the R packages will be loaded and will install the required R packages:

```
# Folder for R packages 
R_PATH=~/path/to/R/library
mkdir -p $R_PATH

# Install package
Rscript ./scrnabox.slurm/soft/R/install_packages_scrnabox.R $R_PATH
```
Alternatively, users can install the packages manually. The R packages required for scRNAbox are shown [here](https://github.com/neurobioinfo/scrnabox/blob/main/scrnabox.slurm/soft/R/R.library.ini).
 - - - -
Upon completing the installation procedures, users can proceed with the scRNAbox pipeline.



