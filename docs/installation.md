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

`scrnabox.slurm` requires that `cellranger` and `R` are also installed on the HPC system. In addition, the following R packages must be loaded: `'Seurat','ggplot2', 'dplyr', 'foreach', 'doParallel', 'Matrix', 'DoubletFinder','cowplot','clustree'`. Then, install the `'scrnaboxR'` R package by running the following command: 
```
devtools::install_github("neurobioinfo/scrnabox/scrnaboxR")
```

Please note that all R packages must be loaded to a common R library. Shown below is an example of how to load packages into a common library in R.
```
R_LIB_PATH=“Path_to_R_library”
.libPaths(R_LIB_PATH)
devtools::install_github("neurobioinfo/scrnabox/scrnaboxR")
```

Upon installing `scrnabox.slurm`,`cellranger`, `R`, and the required R packages, users can run the pipeline initiation Step and define their desired Analysis Track (SCRNA or HTO for Standard scRNAseq or Cell Hashtag scRNAseq, respectively) using the following command:
```
bash ./scrnabox.slurm/launch_scrnabox.sh \
-d ./working_directory \
--steps 0 \
--method SCRNA
```

After initiating the pipeline, the structure of the working directory should be as follows:

```
├── working_directory
    ├── job_info
        ├── configs
        ├── logs
        ├── parameters
```
- The `configs/` directory contains the `scrnabox_config.ini` file which allows users to specify their job allocations (memory, threads, and walltime) for each Analytical Step using the Slurm Workload Manager <br /> 
- The `logs/` directory records the events of each Analytical Step <br />
- The `parameters/` directory contains adjustable, Step-specific text files which allow users to define the execution parameters for each Analytical Step <br />

Users must then navigate to the `scrnabox_config.ini` file in `~/working_directory/job_info/configs` to define the location of their R library (`R_LIB_PATH=`), their version of R (`R_VERSION=`), and the location of CellRanger (`MODULECELLRANGER=`). For example: 

```
MODULECELLRANGER=mugqic/cellranger/5.0.1
R_VERSION=4.2.1
R_LIB_PATH=Path_to_R_library
```

Upon completing the installation procedures, users can proceed with the scRNAbox pipeline using either the [Standard scRNAseq Analysis Track](SCRNA.md) or [Cell Hashtag scRNAseq Analysis Track](HTO.md). 



