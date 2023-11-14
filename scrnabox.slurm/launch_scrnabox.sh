#!/bin/bash

# The pipeline is done as part of Project Dark Genome 
# Copyright belong MNI BIOINFO CORE (https://github.com/neurobioinfo)
# The pipeline is written by Saeid Amiri (saeid.amiri@mcgill.ca)

VERSION=0.1.52;
DATE0=2023-11-12
echo -e "scrnabox pipeline version $VERSION"
# updated on $DATE0"

# ===============================================
# default variables values
# ===============================================
unset SAMPLE OUTPUT_DIR PIPELINE_HOME QUEUE ACCOUNT STEP8rundge STEP8addmeta STEP8i

#Assuming script is in root of PIPELINE_HOME
#set queue-specific values, depending on system check
# submit_cmd="bash"
if [  -z "${submit_cmd}"  ]; then
   submit_cmd="sbatch"
fi
# echo $submit_cmd

if [[ $submit_cmd =~ sbatch ]]; then
   QUEUE="sbatch"          # default job scheduler: qsub
elif [[ $submit_cmd =~ bash ]]; then
    QUEUE="bash"
else
    echo "The pipeline is testing under slurm system,"
    echo "please choose sbatch for slurm."   
    exit 42
fi


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
          "\t\t--steps  =  Specify what steps, e.g., 2 to run just step 2, 2-4, run steps 2 through 4)\n" 
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
if ! options=$(getopt --name pipeline --alternative --unquoted --options hs:d:t:m:vw:f:S:c:a:x: --longoptions dir:,steps:,method:,sinfo:,markergsea:,knownmarkers:,referenceannotation:,annotate:,addmeta:,rundge:,dgelist:,genotype:,celltype:,msd:,cont:,seulist:,verbose,help,rcheck -- "$@")
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
# MODULEUSE=/cvmfs/soft.mugqic/CentOS6/modulefiles
# MODULELOAD=mugqic/cellranger/5.0.1
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
if [[ -z $CONT ]]; then  CONT="F"; fi


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
  module load r/$R_VERSION
  Rscript ${PIPELINE_HOME}/soft/Rcheck/rcheck.R $OUTPUT_DIR $R_LIB_PATH $PIPELINE_HOME/soft/Rcheck/list_R_packages.ini $SUMMARY_R_CHECK
  exit
fi
if [ -z $MODE ]; then echo "ERROR: missing mandatory option: --steps ('ALL' to run all, 2 to run step 2, step 2-4, run steps 2 through 4) must be specified"; FOUND_ERROR=1; fi

if (( $FOUND_ERROR )); then echo "Please check options and try again"; exit 42; fi

if [[ -z $SINFO ]]; then  SINFO="F"; fi

if [ $MODE == 'ALL' ]; then  MODE00=`echo {2..10}`; MODE0=`eval echo $MODE00`; fi 
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
    echo "config file already exists in $OUTPUT_DIR. Overwrite? y|n"
    read answer
    answer=${answer,,}; answer=${answer:0:1}
    if [[ $answer == "y" ]]; then 
      echo "NOTE: the config file and parameters are Overwritten."
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
        if [[ $submit_cmd =~ sbatch ]]; then
            echo -e  "--------running with scheduler-----------------" >> $EXPECTED_DONE_FILES
        elif [[ $submit_cmd =~ bash ]]; then
            echo -e  "--------running with the container-----------------" >> $EXPECTED_DONE_FILES
        fi 
        echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
        echo "The Output is under ${OUTPUT_DIR}/job_info/" >> $EXPECTED_DONE_FILES
fi 

declare -A THREADS_ARRAY
declare -A  WALLTIME_ARRAY
declare -A  MEM_ARRAY
source $OUTPUT_DIR/job_info/configs/scrnabox_config.ini

export JOB_OUTPUT_DIR=$OUTPUT_DIR/job_info
export EXPECTED_DONE_FILES=$JOB_OUTPUT_DIR/summary_report.txt
# export STEP8_CON=$PIPELINE_HOME/configs/step8_contrast.txt
chmod 775 $EXPECTED_DONE_FILES    
# export MODULEUSE=$MODULEUSE
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

echo " # IT IS A temp FILE. DO NOT EDIT THIS FILE DIRECTLY."  > $TEMPCONFIG
echo STEP8addmeta=$STEP8addmeta  >> $TEMPCONFIG
echo STEP8rundge=$STEP8rundge   >> $TEMPCONFIG
echo MSD=$MSD  >> $TEMPCONFIG
echo CONT=$CONT  >> $TEMPCONFIG
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
if [[ $MODULEUSE ]]; then 
    echo CELLRANGER=$CELLRANGER >> $TEMPCONFIG
fi
if [[ $MODULEUSE ]]; then 
    echo CELLRANGER_VERSION=${CELLRANGER_VERSION} >> $TEMPCONFIG
fi


if  [[  ${MODE0[@]}  =~  integrate  ]]  ; then
   bash ${PIPELINE_HOME}/launch/launch_scrnabox_integ.sh
elif  [[  ${SCRNA_METHOD} =~ SCRNA ]] ; then
   bash ${PIPELINE_HOME}/launch/launch_scrnabox_scrna.sh
elif [[  ${SCRNA_METHOD} =~ HTO ]] ; then
   bash ${PIPELINE_HOME}/launch/launch_scrnabox_hto.sh
fi


exit 0
