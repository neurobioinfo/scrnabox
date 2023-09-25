## Job configurations for the scRNAbox pipeline
- [Introduction](#introduction)
- [General configurations](#general-configurations)
- [Analytical Step configurations](#analytical-step-configurations)
- - - -
## Introduction
Upon running the pipeline initiation Step (Step 0), the `scrnabox_config.ini` file is automatically deposited into `~/working_directory/job_info/configs`:

```
ACCOUNT=HPC_account_owner_name
MODULEUSE=/cvmfs/soft.mugqic/CentOS6/modulefiles
MODULECELLRANGER=mugqic/cellranger/5.0.1
R_VERSION=4.2.1
R_LIB_PATH=/path/to/R/library
SCRNA_METHOD=SCRNA
############# 

############# [step2]
#THREADS_ARRAY["step_2"]=4
#MEM_ARRAY["step_2"]=16g
#WALLTIME_ARRAY["step_2"]=00-05:00 

############# [step3]
#THREADS_ARRAY["step_3"]=4
#MEM_ARRAY["step_3"]=16g
#WALLTIME_ARRAY["step_3"]=00-05:00 

############# [step4]
#THREADS_ARRAY["step_4"]=4
#MEM_ARRAY["step_4"]=45g
#WALLTIME_ARRAY["step_4"]=00-05:00 

############# [step5]
#THREADS_ARRAY["step_5"]=4
#MEM_ARRAY["step_5"]=45g
#WALLTIME_ARRAY["step_5"]=00-05:00 

#############[step6]
#THREADS_ARRAY["step_6"]=4
#MEM_ARRAY["step_6"]=16g
#WALLTIME_ARRAY["step_6"]=00-05:00 

############# [step7]
#THREADS_ARRAY["step_7marker"]=4
#MEM_ARRAY["step_7marker"]=40g
#WALLTIME_ARRAY["step_7marker"]=00-1:00 

#THREADS_ARRAY["step_7enrich"]=10
#MEM_ARRAY["step_7enrich"]=15g
#WALLTIME_ARRAY["step_7enrich"]=00-01:00 

#THREADS_ARRAY["step_7fta"]=4
#MEM_ARRAY["step_7fta"]=150g
#WALLTIME_ARRAY["step_7fta"]=00-09:00 

############# [step8]
#THREADS_ARRAY["step_8_dgelist"]=4
#MEM_ARRAY["step_8_dgelist"]=40g
#WALLTIME_ARRAY["step_8_dgelist"]=00-12:00 

############# [step8]
#THREADS_ARRAY["step_8_cont"]=10
#MEM_ARRAY["step_8_cont"]=40g
#WALLTIME_ARRAY["step_8_cont"]=00-12:00 
```

The `scrnabox_config.ini` defines the parameters for each Job submission to the [Slurm Workload Manager](https://slurm.schedmd.com/). Prior to performing their analysis with the scRNAbox pipeline, users should modify the this file to adjust the [General configurations](#general-configurations) and [Analytical Step configurations](#analytical-step-configurations). 

- - - -

## General configurations
The following **general configuration** parameters must be adjusted by the user:

1) **ACCOUNT**: the name of the account for the HPC system; <br />
2) **MODULEUSE**: <br />
3) **MODULECELLRANGER**:  <br />
4) **R_VERSION**: the version of R installed on the HPC system; <br />
5) **R_LIB_PATH**: the path to the R library containing the required R packages for running scRNAbox; <br />
6) **SCRNA_METHOD**: The Analytical Track used for analysis. This will be automatically defined upon running Step 0. <br />

- - - -

## Analytical Step configurations
Each Analytical Step of the scRNAbox pipeline has three configuartion parameters:

1) **THREADS_ARRAY**: number of CPUs for the Job; <br />
2) **MEM_ARRAY**: Amount of memory (RAM) for the Job; <br />
3) **WALLTIME_ARRAY**: Amount of time for the Job.  <br />

In the original `scrnabox_config.ini` file, the configuration parameters for each Analytical Step will be commented out and the **default** configurations will be used. If users need to change the configurations for any Analytical Step, they may uncomment the line of code and define the parameter:

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