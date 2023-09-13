# Installation
`scrnabox.slurm` is written in bash and can be used with any Slurm system. To download the latest version of `scrnabox.slurm` (v0.1.35) run the following command: 
```
wget https://github.com/neurobioinfo/scrnabox/releases/download/v0.1.35/scrnabox.slurm.zip
unzip scrnabox.slurm.zip
```

For a description of the options for running `scrnabox.slurm` run the following command:
```
bash ./scrnabox.slurm/launch_scrnabox.sh -h 
```

`scrnabox.slurm` requires that `R` and `cellranger` are also installed on the HPC system. In addition, the following R packages must be loaded: `'Seurat','ggplot2', 'dplyr', 'foreach', 'doParallel', 'Matrix', 'DoubletFinder','cowplot','clustree'`. Then, install the `'scrnaboxR'` R package by running the following command: 
```
devtools::install_github("neurobioinfo/scrnabox/scrnaboxR")
```
`'scrnaboxR'` provides a collection of functions for conducting enrichment analysis and other analyses associated with scRNAseq data. It serves as a companion to scRNAbox, offering a range of tools and functionalities to enhance scRNAseq data analysis. 

Please note that all R packages must be loaded into a common R library. For further instructions regarding the preparation of an R library please visit scRNAbox's [documentation](https://neurobioinfo.github.io/scrnabox/site/). Users must then define the location of their R library (`R_LIB_PATH=`), their version of R (`R_VERSION=`), and the location of CellRanger (`MODULECELLRANGER=`) in the `scrnabox_config.ini` file which is deposited into the working directory upon running the pipeline initation Step:
```
bash ./scrnabox.slurm/launch_scrnabox.sh \
-d ./working_directory \
--steps 0 \
--method SCRNA
