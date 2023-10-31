# Step 0: scRNAbox pipeline set up
Step 0 initiates the pipeline and sets up the working directory by depositing the required files.
 - - - -
## Running Step 0
Create a dedicated folder for the analysis (hereafter referred to as the working directory). Then, export the path to the working directory and the path to `scrnabox.slurm`:
```
mkdir working_directory
cd /pathway/to/working_directory

export SCRNABOX_HOME=/pathway/to/scrnabox.slurm
export SCRNABOX_PWD=/pathway/to/working_directory
```

Next, run Step 0 and choose whether to use the standard analysis track (`--method SCRNA`) or HTO analysis track (`--method HTO`) using the following command from the working directory:
```
cd /pathway/to/working_directory 

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 0 \
--method SCRNA
```

After running Step 0, the structure of the working directory should be:
```
 working_directory
 └── job_info
     ├── configs
     ├── logs
     └── parameters
```
- The `configs/` directory contains the `scrnabox_config.ini` file which allows users to specify their job allocations (memory, threads, and walltime) for each analytical step using the Slurm Workload Manager; <br /> 
- The `logs/` directory records the events of each analytical step; <br />
- The `parameters/` directory contains adjustable, step-specific text files which allow users to define the execution parameters for each analytical step. <br />

Next, navigate to the `scrnabox_config.ini` file in `~/working_directory/job_info/configs` to define the path to the R library, the version of R, and the path to CellRanger:

```
MODULECELLRANGER=mugqic/cellranger/5.0.1
R_VERSION=4.2.1
R_LIB_PATH=/path/to/R/library
```
**Note:** For more information, please see the [Job cofigurations](config.md) sections of the scRNAbox documentation.
 - - - -
Upon completing the setup procedures, users can run their analysis using the scRNAbox pipeline. 




