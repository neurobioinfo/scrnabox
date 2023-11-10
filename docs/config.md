## Job configurations for the scRNAbox pipeline
- [Introduction](#introduction)
- [General configurations](#general-configurations)
- [Step configurations](#step-configurations)
- - - -
## Introduction
Upon running Step 0, the `scrnabox_config.ini` file is automatically deposited into `~/working_directory/job_info/configs`:

```
ACCOUNT=account-name
MODULEUSE=/path/to/environmental/module (e.g. /cvmfs/soft.mugqic/CentOS6/modulefiles)
CELLRANGER=/path/to/cellranger/from/module/directory (e.g. mugqic/cellranger)
CELLRANGER_VERSION=5.0.1
R_VERSION=4.2.1
R_LIB_PATH=/path/to/R/library
SCRNA_METHOD=SCRNA
############# 

############# [step2]
# THREADS_ARRAY["step_2"]=10
# MEM_ARRAY["step_2"]=16g
# WALLTIME_ARRAY["step_2"]=00-05:00 

############# [step3]
# THREADS_ARRAY["step_3"]=10
# MEM_ARRAY["step_3"]=16g
# WALLTIME_ARRAY["step_3"]=00-05:00 

############# [step4]
# THREADS_ARRAY["step_4"]=10
# MEM_ARRAY["step_4"]=45g
# WALLTIME_ARRAY["step_4"]=00-05:00 

############# [step5]
# THREADS_ARRAY["step_5"]=10
# MEM_ARRAY["step_5"]=16g
# WALLTIME_ARRAY["step_5"]=00-05:00 

#############[step6]
# THREADS_ARRAY["step_6"]=10
# MEM_ARRAY["step_6"]=16g
# WALLTIME_ARRAY["step_6"]=00-05:00 

############# [step7]
# THREADS_ARRAY["step_7_markergsea"]=4
# MEM_ARRAY["step_7_markergsea"]=16g
# WALLTIME_ARRAY["step_7_markergsea"]=00-05:00

# THREADS_ARRAY["step_7_knownmarkers"]=4
# MEM_ARRAY["step_7_knownmarkers"]=16g
# WALLTIME_ARRAY["step_7_knownmarkers"]=00-05:00

# THREADS_ARRAY["step_7_referenceannotation"]=10
# MEM_ARRAY["step_7_referenceannotation"]=150g
# WALLTIME_ARRAY["step_7_referenceannotation"]=00-12:00

# THREADS_ARRAY["step_7_annotate"]=4
# MEM_ARRAY["step_7_annotate"]=16g
# WALLTIME_ARRAY["step_7_annotate"]=00-05:00

############# [step8]
# THREADS_ARRAY["step_8_addmeta"]=4
# MEM_ARRAY["step_8_addmeta"]=16g
# WALLTIME_ARRAY["step_8_addmeta"]=00-12:00 

# THREADS_ARRAY["step_8_rundge"]=4
# MEM_ARRAY["step_8_rundge"]=40g
# WALLTIME_ARRAY["step_8_rundge"]=00-12:00 

############# [integrate]
# THREADS_ARRAY["integrate"]=10
# MEM_ARRAY["integrate"]=40g
# WALLTIME_ARRAY["integrate"]=00-04:00 

```

The `scrnabox_config.ini` defines the parameters for each Job submission to the [Slurm Workload Manager](https://slurm.schedmd.com/). Prior to performing their analysis with the scRNAbox pipeline, users should modify this file to adjust the [General configurations](#general-configurations) and [Analytical Step configurations](#analytical-step-configurations). 

- - - -

## General configurations
The following parameters must be adjusted by the user:

1) **ACCOUNT**: the name of the account holder for the HPC system; <br />
2) **MODULEUSE**: the path to the environmental module; <br />
3) **CELLRANGER**: the path to CellRanger from the environmental module directory; <br />
4) **CELLRANGER_VERSION**: the version of CellRanger; <br />
5) **R_VERSION**: the version of R installed on the HPC system (must be v4.2 or greater); <br />
6) **R_LIB_PATH**: the path to the R library containing the required R packages for running scRNAbox; <br />
7) **SCRNA_METHOD**: The analysis track to use (SCRNA or HTO). This will be automatically defined upon running Step 0. <br />

- - - -

## Step configurations
Each step of the scRNAbox pipeline has three configuartion parameters:

1) **THREADS_ARRAY**: number of CPUs for the job; <br />
2) **MEM_ARRAY**: amount of memory (RAM) for the job; <br />
3) **WALLTIME_ARRAY**: amount of time for the job.  <br />

In the original `scrnabox_config.ini` file, the configuration parameters for each step will be commented out and the **default** configurations will be used. If users need to change the configurations for any step, they can uncomment the line of code and define the parameter:

```
# orginal scrnabox_config.ini
#THREADS_ARRAY["step_2"]=4
#MEM_ARRAY["step_2"]=16g
#WALLTIME_ARRAY["step_2"]=00-05:00 

# modified scrnabox_config.ini
THREADS_ARRAY["step_2"]=8
MEM_ARRAY["step_2"]=24g
WALLTIME_ARRAY["step_2"]=00-10:00 
```