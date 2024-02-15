#!/bin/bash

# The pipeline is done as part of Project Dark Genome 
# Copyright belong MNI BIOINFO CORE (https://github.com/neurobioinfo)
# The pipeline is written by Saeid Amiri (saeid.amiri@mcgill.ca)

VERSION=0.1.53.31
DATE0=2024-02-13
echo -e "scrnabox pipeline version $VERSION"

# ===============================================
# default variables values
# ===============================================
unset SAMPLE OUTPUT_DIR PIPELINE_HOME QUEUE ACCOUNT STEP8rundge STEP8addmeta STEP8i

#Assuming script is in root of PIPELINE_HOME
#set queue-specific values, depending on system check
# submit_cmd="bash"
# if [  -z "${submit_cmd}"  ]; then
#    submit_cmd="bash"
# fi
# # echo $submit_cmd

# if [[ $submit_cmd =~ sbatch ]]; then
#    QUEUE="sbatch"          # default job scheduler: qsub
# elif [[ $submit_cmd =~ bash ]]; then
#     QUEUE="bash"
# else
#     echo "The pipeline is testing under slurm system,"
#     echo "please choose sbatch for slurm."   
#     exit 42
# fi


PIPELINE_HOME0=`realpath ${BASH_SOURCE[0]}`
export PIPELINE_HOME=$(cd $(dirname $PIPELINE_HOME0) && pwd -P)

TIMESTAMP=`date +%FT%H.%M.%S`


# create function to handle error messages
# ===============================================
Usage() {
	echo
  echo "------------------- " 
	echo -e "Usage:\t$0 [arguments]"
	echo -e "\tmandatory arguments:\n" \
          "\t\t-d  (--dir)  = Working directory (where all the outputs will be printed) (give full path)\n" \
          "\t\t--steps  =  Specify what steps, e.g., 2 to run step 2.\n" 
	echo -e "\toptional arguments:\n " \
          "\t\t-h  (--help)  = See helps regarding the pipeline arguments. \n" \
          "\t\t--method  = Select your preferred method: HTO and SCRNA for hashtag, and Standard scRNA, respectively. \n" \
          "\t\t--msd  = You can get the hashtag labels by running the following code (HTO Step 4). \n" \
          "\t\t--markergsea  = Identify marker genes for each cluster and run marker gene set enrichment analysis (GSEA) using EnrichR libraries (Step 7). \n" \
          "\t\t--knownmarkers  = Profile the individual or aggregated expression of known marker genes. \n" \
          "\t\t--referenceannotation  = Generate annotation predictions based on the annotations of a reference Seurat object (Step 7). \n" \
          "\t\t--annotate  = Add clustering annotations to Seurat object metadata (Step 7). \n" \
          "\t\t--addmeta  = Add metadata columns to the Seurat object (Step 8). \n" \
          "\t\t--rundge  = Perform differential gene expression contrasts (Step 8). \n" \
          "\t\t--seulist  = You can directly call the list of Seurat objects to the pipeline. \n" \
          "\t\t--rcheck  = You can identify which libraries are not installed.  \n \n" \
          "------------------- \n" \
          "For a comprehensive help, visit  https://neurobioinfo.github.io/scrnabox/site/ for documentation. "

echo 
}


# ===============================================
# PARSING ARGUMENTS
# ===============================================
if ! options=$(getopt --name pipeline --alternative --unquoted --options hs:d:t:m:vw:f:S:c:a:x: --longoptions dir:,steps:,method:,jobmode:,sinfo:,markergsea:,knownmarkers:,referenceannotation:,annotate:,addmeta:,rundge:,dgelist:,genotype:,celltype:,msd:,cont:,seulist:,verbose,help,rcheck -- "$@")
then
    # something went wrong, getopt will put out an error message for us
    echo "Error processing options."
    exit 42
fi

# ===============================================
# LOAD & OVERRIDE EXTRA CONFIG FILE FOR PROJECT
# ===============================================
set -- $options

while [ $# -gt 0 ]

do
    case $1 in
    -x| --extra) 
      EXTRA_CONF="$2" ;
      if [ -f $EXTRA_CONF ]; then
        echo "* LOADING EXTRA CONFIG FILE $EXTRA_CONF";
        . $EXTRA_CONF
      else
        echo "ERROR: invalid EXTRA CONFIG file: $EXTRA_CONF";
        echo "Please check options and try again"; exit 42;
      fi
    esac
    shift
done

if [[ -n $ACCOUNT ]]; then ACCOUNT="-A $ACCOUNT"; fi
if [[ $MODULEUSE ]]; then module use $MODULEUSE ; fi
if [[ ${CELLRANGER} ]]; then module load ${CELLRANGER}/${CELLRANGER_VERSION} ; fi

# ===============================================
# LOAD ALL OTHER OPTIONS
# ===============================================
set -- $options

while [ $# -gt 0 ]
do
    case $1 in
    -h| --help) Usage; exit 0;;
    -d| --dir) OUTPUT_DIR="$2" ; shift ;;
    -v| --verbose) VERBOSE=1 ;; 
    --steps) MODE="$2"; shift ;; 
    --method) SCRNA_METHOD0="$2"; shift ;; 
    --jobmode) JOB_MODE="$2"; shift ;;     
    --markergsea) STEP7markergsea="$2"; shift ;;
    --knownmarkers) STEP7knownmarkers="$2"; shift ;;
    --referenceannotation) STEP7referenceannotation="$2"; shift ;;
    --annotate) STEP7annotate="$2"; shift ;;
    --addmeta) STEP8addmeta="$2"; shift ;;
    --rundge) STEP8rundge="$2"; shift ;;
    --msd) MSD="$2"; shift ;;
    --sinfo) SINFO="$2"; shift ;;
    --seulist) SEULIST="$2"; shift ;;    
    --rcheck) RCHECK="$2"; shift ;;  
    (--) shift; break;;
    (-*) echo "$0: error - unrecognized option $1" 1>&2; exit 42;;
    (*) break;;
    esac
    shift
done


if [[ -z $STEP8rundge ]]; then  STEP8rundge="F"; fi
if [[ -z $STEP8addmeta ]]; then  STEP8addmeta="F"; fi
if [[ -z $MSD ]]; then  MSD="F"; fi


source $PIPELINE_HOME/tools/utils.sh
FOUND_ERROR=0

# ===============================================
# CHECKING VARIABLES
# ===============================================
#check to ensure all mandatory arguments have been entered

if [ -z $OUTPUT_DIR ]; then echo "ERROR: missing mandatory option: -d (--dir) must be specified"; FOUND_ERROR=1; fi

if [ $RCHECK ]; then 
  echo -e "You pass << $R_LIB_PATH >> to pipeline"
  echo "The pipeline is verifying the installation status of all R libraries."
  source $OUTPUT_DIR/job_info/configs/scrnabox_config.ini 
  SUMMARY_R_CHECK=$OUTPUT_DIR/job_info/summary_r_check.txt
  if [[ -f $SUMMARY_R_CHECK ]]; then 
    rm $SUMMARY_R_CHECK
  fi
  touch $SUMMARY_R_CHECK
  echo "It is the result of verifying the R library." >>$SUMMARY_R_CHECK 
  echo "The list of libraries not installed will be displayed below." >>$SUMMARY_R_CHECK 
  if [[ $JOB_MODE =~ slurm ]]; then
    module load r/$R_VERSION
  fi
  Rscript ${PIPELINE_HOME}/soft/Rcheck/rcheck.R $OUTPUT_DIR $R_LIB_PATH $PIPELINE_HOME/soft/Rcheck/list_R_packages.ini $SUMMARY_R_CHECK
  exit
fi

# if [ -z $MODE ]; thexn echo "ERROR: missing mandatory option: --steps (2 to run step 2) must be specified"; FOUND_ERROR=1; fi

if (( $FOUND_ERROR )); then echo "Please check options and try again"; exit 42; fi

if [[ -z $SINFO ]]; then  SINFO="F"; fi

# if [ $MODE == 'ALL' ]; then  MODE00=`echo {2..10}`; MODE0=`eval echo $MODE00`; fi 
if [[ "$MODE" == *"-"* ]]; then
  MODE00=`echo $MODE | sed  "s/-/../g"  | awk '{print "{"$0"}"}'`
  MODE0=`eval echo $MODE00`
else 
  MODE0=$MODE
fi


# STEP 0: RUN setting 
# ===============================================
#


if [[ ${MODE0[@]} == 0 ]]; then 
  JOB_OUTPUT_DIR=$OUTPUT_DIR/job_info
  if [[ -s $OUTPUT_DIR/job_info/configs ]]; then
    # echo "config file already exists in $OUTPUT_DIR. Overwrite? y|n"
    # read answer
    read -p "config file already exists in $OUTPUT_DIR. Overwrite? y|n: "$'\n'  answer
    # answer=${answer,,}; answer=${answer:0:1}
    # echo $answer
    if [[ $answer =~ y ]]; then 
      echo "NOTE: the config file and parameters are Overwritten."
        if [ ! -d $JOB_OUTPUT_DIR ]; then   
        mkdir -p $OUTPUT_DIR/job_info
        mkdir -p $OUTPUT_DIR/job_info/parameters
        mkdir -p $OUTPUT_DIR/job_info/logs
        mkdir -p $OUTPUT_DIR/job_info/configs
        mkdir -p $OUTPUT_DIR/job_info/.tmp        
        fi
        if [ -d $OUTPUT_DIR/job_info ]; then
            rm -r $OUTPUT_DIR/job_info
            mkdir -p $OUTPUT_DIR/job_info
            mkdir -p $OUTPUT_DIR/job_info/parameters
            mkdir -p $OUTPUT_DIR/job_info/logs
            mkdir -p $OUTPUT_DIR/job_info/configs
            mkdir -p $OUTPUT_DIR/job_info/.tmp
        fi
        touch $JOB_OUTPUT_DIR/summary_report.txt
        if [[  ${SCRNA_METHOD0} =~ SCRNA ]] ; then
          cp $PIPELINE_HOME/scrna/configs/scrnabox_config.ini $OUTPUT_DIR/job_info/configs/ 
          cp $PIPELINE_HOME/scrna/pars/* $OUTPUT_DIR/job_info/parameters/ 
        elif [[  ${SCRNA_METHOD0} =~ HTO ]] ; then
          cp $PIPELINE_HOME/hto/configs/scrnabox_config.ini $OUTPUT_DIR/job_info/configs/ 
          cp $PIPELINE_HOME/hto/pars/* $OUTPUT_DIR/job_info/parameters/ 
        elif [[ -z  ${SCRNA_METHOD0}  ]]; then
          cp $PIPELINE_HOME/scrna/configs/scrnabox_config.ini $OUTPUT_DIR/job_info/configs/ 
          cp $PIPELINE_HOME/scrna/pars/* $OUTPUT_DIR/job_info/parameters/ 
          echo "NOTE: You didn't specify the method, so the pipeline choose << SCRNA seq >>"
        else
          echo "Error: the pipeline is not developed for ${SCRNA_METHOD0}"
          exit 42
          fi
    else
     echo "NOTE: the pipeline is using the existing config file and parameters."
    fi
  else 
    echo "The configuration files do not exist, the pipeline will create them during execution."
        if [ ! -d $JOB_OUTPUT_DIR ]; then 
        mkdir -p $OUTPUT_DIR/job_info
        mkdir -p $OUTPUT_DIR/job_info/parameters
        mkdir -p $OUTPUT_DIR/job_info/logs
        mkdir -p $OUTPUT_DIR/job_info/configs
        mkdir -p $OUTPUT_DIR/job_info/.tmp                
        fi
        touch $JOB_OUTPUT_DIR/summary_report.txt
        if [[  ${SCRNA_METHOD0} =~ SCRNA ]] ; then
          cp $PIPELINE_HOME/scrna/configs/scrnabox_config.ini $OUTPUT_DIR/job_info/configs/
          cp $PIPELINE_HOME/scrna/pars/* $OUTPUT_DIR/job_info/parameters/ 
        elif [[  ${SCRNA_METHOD0} =~ HTO ]] ; then
          cp $PIPELINE_HOME/hto/configs/scrnabox_config.ini $OUTPUT_DIR/job_info/configs/
          cp $PIPELINE_HOME/hto/pars/* $OUTPUT_DIR/job_info/parameters/ 
        else
          echo "Error: the pipeline is not developed for ${SCRNA_METHOD}"
          exit 42
          fi
  fi
  if [[   ${SINFO}  =~  T ]]; then
    if  [ ! -d $OUTPUT_DIR/samples_info ]; then   echo 'ERROR: The pipeline can not find samples_info directory,  '; exit 42 ; fi 
    if  [ -d ${OUTPUT_DIR}/step1 ]; then 
      rm ${OUTPUT_DIR}/job_info/.tmp/sample_dir.list
      rm ${OUTPUT_DIR}/job_info/.tmp/sample.list
      search_dir=${OUTPUT_DIR}/step1
      for entry in "$search_dir"/*
        do
          echo "$entry" >> ${OUTPUT_DIR}/job_info/.tmp/sample_dir.list
          echo $(basename  "$entry") >> ${OUTPUT_DIR}/job_info/.tmp/sample.list
        done
   fi
  fi 
        EXPECTED_DONE_FILES=$JOB_OUTPUT_DIR/summary_report.txt
        echo -e  "-------------------------------------------" >> $EXPECTED_DONE_FILES  
        echo -e  "--------Pipeline is set up for ${SCRNA_METHOD0}-----------------" $VERSION >> $EXPECTED_DONE_FILES
        if [[ $JOB_MODE =~ slurm ]]; then
            echo -e  "--------running with scheduler-----------------" >> $EXPECTED_DONE_FILES
            echo JOB_MODE=slurm >> $OUTPUT_DIR/job_info/configs/scrnabox_config.ini
        elif [[ $JOB_MODE =~ local ]]; then
            echo -e  "--------running without scheduler-----------------" >> $EXPECTED_DONE_FILES
            echo JOB_MODE=local >> $OUTPUT_DIR/job_info/configs/scrnabox_config.ini
            remove_argument $OUTPUT_DIR/job_info/configs/scrnabox_config.ini
        else
            JOB_MODE=local
            echo -e  "--------running without scheduler-----------------" >> $EXPECTED_DONE_FILES
            echo JOB_MODE=local >> $OUTPUT_DIR/job_info/configs/scrnabox_config.ini            
            remove_argument $OUTPUT_DIR/job_info/configs/scrnabox_config.ini
        fi 
        echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
        echo "The Output is under ${OUTPUT_DIR}/job_info/" >> $EXPECTED_DONE_FILES
        echo "NOTE: the scrnbox is setup with  ${JOB_MODE} mode"
#  exit 0
fi 

declare -A THREADS_ARRAY
declare -A  WALLTIME_ARRAY
declare -A  MEM_ARRAY
source $OUTPUT_DIR/job_info/configs/scrnabox_config.ini

export JOB_OUTPUT_DIR=$OUTPUT_DIR/job_info
export EXPECTED_DONE_FILES=$JOB_OUTPUT_DIR/summary_report.txt
chmod 775 $EXPECTED_DONE_FILES
if [[ $MODULEUSE ]]; then export MODULEUSE=$MODULEUSE ; fi
export CELLRANGER=${CELLRANGER}
export ACCOUNT=$ACCOUNT
export OUTPUT_DIR=$OUTPUT_DIR

echo -e "NOTE: the pipeline is running  << $SCRNA_METHOD seq >>"

TEMPCONFIG=$OUTPUT_DIR/job_info/.tmp/temp_config.ini
touch $TEMPCONFIG
if [[ -f "TEMPCONFIG" ]]; then 
  rm $TEMPCONFIG
fi

if [  -z "${JOB_MODE}"  ]; then
   JOB_MODE=local
fi

if [[ $JOB_MODE =~ slurm ]]; then
   QUEUE="sbatch"
elif [[ $JOB_MODE =~ local ]]; then
    QUEUE="bash"
else
    echo "The pipeline is tested under slurm system, and linux"
    echo "please choose choose slurm or local"   
    exit 42
fi

echo " # IT IS A temp FILE. DO NOT EDIT THIS FILE DIRECTLY."  > $TEMPCONFIG
echo STEP8addmeta=$STEP8addmeta  >> $TEMPCONFIG
echo STEP8rundge=$STEP8rundge   >> $TEMPCONFIG
echo MSD=$MSD  >> $TEMPCONFIG
echo SEULIST=$SEULIST   >> $TEMPCONFIG
echo STEP7markergsea=$STEP7markergsea   >> $TEMPCONFIG
echo STEP7knownmarkers=$STEP7knownmarkers   >> $TEMPCONFIG
echo STEP7referenceannotation=$STEP7referenceannotation   >> $TEMPCONFIG
echo STEP7annotate=$STEP7annotate   >> $TEMPCONFIG
echo SINFO=$SINFO  >> $TEMPCONFIG
echo OUTPUT_DIR=$OUTPUT_DIR >> $TEMPCONFIG
echo JOB_OUTPUT_DIR=$JOB_OUTPUT_DIR  >> $TEMPCONFIG
echo EXPECTED_DONE_FILES=$EXPECTED_DONE_FILES   >> $TEMPCONFIG
echo MODE=$MODE >> $TEMPCONFIG
echo QUEUE=$QUEUE >> $TEMPCONFIG
echo VERSION=$VERSION >> $TEMPCONFIG
if [[ $MODULEUSE ]]; then 
    echo MODULEUSE=$MODULEUSE >> $TEMPCONFIG
fi
if [[ $CELLRANGER ]]; then 
    echo CELLRANGER=$CELLRANGER >> $TEMPCONFIG
fi
if [[ ${CELLRANGER_VERSION} ]]; then 
    echo CELLRANGER_VERSION=${CELLRANGER_VERSION} >> $TEMPCONFIG
fi

# if [  -z "${submit_cmd}"  ]; then
  #  submit_cmd="bash"
# fi
# echo $submit_cmd

# if [[ $submit_cmd =~ sbatch ]]; then
#    QUEUE="sbatch"          # default job scheduler: qsub
# elif [[ $submit_cmd =~ bash ]]; then
#     QUEUE="bash"
# else
#     echo "The pipeline is testing under slurm system,"
#     echo "please choose sbatch for slurm."   
#     exit 42
# fi

if    [[  ${MODE0[@]}  =~  integrate  ]] && [[ ${QUEUE} =~ sbatch ]] ; then
   bash ${PIPELINE_HOME}/launch/launch_scrnabox_integ_slurm.sh
elif  [[  ${MODE0[@]}  =~  integrate  ]] && [[ ${QUEUE} =~ bash   ]] ; then
   bash ${PIPELINE_HOME}/launch/launch_scrnabox_integ_local.sh   
elif  [[  ${SCRNA_METHOD} =~ SCRNA    ]] && [[ ${QUEUE} =~ sbatch ]] ; then
   bash ${PIPELINE_HOME}/launch/launch_scrnabox_scrna_slurm.sh
elif  [[  ${SCRNA_METHOD} =~ SCRNA    ]] && [[ ${QUEUE} =~ bash   ]]; then
   bash ${PIPELINE_HOME}/launch/launch_scrnabox_scrna_local.sh
elif  [[  ${SCRNA_METHOD} =~ HTO      ]] && [[ ${QUEUE} =~ sbatch ]] ; then
   bash ${PIPELINE_HOME}/launch/launch_scrnabox_hto_slurm.sh
elif  [[  ${SCRNA_METHOD} =~ HTO      ]] && [[ ${QUEUE} =~ bash   ]] ; then
   bash ${PIPELINE_HOME}/launch/launch_scrnabox_hto_local.sh   
fi
echo -e " \n"
exit 0
