#!/bin/bash

# The pipeline is done as part of Project Dark Genome 
# Copyright belong MNI BIOINFO CORE (https://github.com/neurobioinfo)
# The pipeline is written by Saeid Amiri (saeid.amiri@mcgill.ca)



VERSION=0.1.41;
DATE0=2023-10-06
echo -e "scrnabox pipeline version $VERSION, dev"
# updated on $DATE0"

# ===============================================
# default variables values
# ===============================================
unset SAMPLE OUTPUT_DIR PIPELINE_HOME QUEUE ACCOUNT STEP8dgelist STEP8m STEP8i

#Assuming script is in root of PIPELINE_HOME
submit_cmd=sbatch


PIPELINE_HOME0=`realpath ${BASH_SOURCE[0]}`
export PIPELINE_HOME=$(cd $(dirname $PIPELINE_HOME0) && pwd -P)


# echo $PIPELINE_HOME
# echo $0
# echo $PIPELINE_HOME0
# SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
# echo $SCRIPTPATH
# echo "[$0] vs. [${BASH_SOURCE[0]}]"
# cd -- "$(dirname "./launch_scrnabox.sh")" >/dev/null 2>&1 ; pwd -P 

# SCRIPT=$(realpath "$0")
# SCRIPTPATH=$(dirname "$SCRIPT")
# echo $SCRIPTPATH

TIMESTAMP=`date +%FT%H.%M.%S`


# create function to handle error messages
# ===============================================
Usage() {
	echo
	echo -e "Usage:\t$0 [arguments]"
	echo -e "\tmandatory arguments:\n" \
          "\t\t-d  (--dir)  = Working directory (where all the outputs will be printed) (give full path)\n" \
          "\t\t--steps  =  Specify what steps, e.g., 2 to run just step 2, 2-4, run steps 2 through 4)\n" 
	echo -e "\toptional arguments:\n " \
          "\t\t-h  (--help)  = See helps regarding the pipeline arguments. \n" \
          "\t\t--method  = Select your preferred method: HTO and SCRNA for hashtag, and Standard scRNA, respectively. \n" \
          "\t\t--msd  = You can get the hashtag labels by running the following code \n" \
          "\t\t--marker  = Find marker. \n" \
          "\t\t--sinfo  = Do you need sample info? \n" \
          "\t\t--fta  = FindTransferAnchors \n" \
          "\t\t--enrich  = Annotation \n" \
          "\t\t--dgelist  = creates a DGEListobject from a table of counts obtained from seurate objects. \n" \
          "\t\t--genotype  = Run the genotype contrast. \n" \
          "\t\t--celltype  = Run the Genotype-cell contrast. \n" \
          "\t\t--cont  = You can directly call the contrast to the pipeline.  \n" \
          "\t\t--seulist                = You can directly call the list of seurat objects to the pipeline.  \n" 
echo 
}
          # "\t\t-v  (--verbose)  = set verbosity level [CURRENT \"$VERBOSE\"]\n" 


# ===============================================
# PARSING ARGUMENTS
# ===============================================
if ! options=$(getopt --name pipeline --alternative --unquoted --options hs:d:t:m:vw:f:S:c:a:x: --longoptions dir:,steps:,method:,marker:,sinfo:,fta:,enrich:,dgelist:,genotype:,celltype:,msd:,cont:,seulist:,verbose,help -- "$@")
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
if [[ $MODULECELLRANGER ]]; then module load $MODULECELLRANGER ; fi


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
    --sinfo) SINFO="$2"; shift ;;
    --marker) STEP7marker="$2"; shift ;;
    --fta) STEP7fta="$2"; shift ;;
    --enrich) STEP7enrich="$2"; shift ;;    
    --dgelist) STEP8dgelist="$2"; shift ;;
    --genotype) STEP8m="$2"; shift ;;
    --celltype) STEP8i="$2"; shift ;;
    --msd) MSD="$2"; shift ;;
    --cont) CONT="$2"; shift ;;
    --seulist) SEULIST="$2"; shift ;;    
    (--) shift; break;;
    (-*) echo "$0: error - unrecognized option $1" 1>&2; exit 42;;
    (*) break;;
    esac
    shift
done



if [[ -z $STEP8m ]]; then  STEP8m="F"; fi
if [[ -z $STEP8dgelist ]]; then  STEP8dgelist="F"; fi
if [[ -z $STEP8i ]]; then  STEP8i="F"; fi
if [[ -z $MSD ]]; then  MSD="F"; fi
if [[ -z $CONT ]]; then  CONT="F"; fi


# if [[ -z $NFRNAL ]]; then  NFRNAL="F"; fi
# if [[ -z $NFRNAU ]]; then  NFRNAU="F"; fi
# if [[ -z $NCRNAL ]]; then  NCRNAL="F"; fi
# if [[ -z $NCRNAU ]]; then  NCRNAU="F"; fi
# if [[ -z $PMTL ]]; then  PMTL="F"; fi
# if [[ -z $PMTU ]]; then  PMTU="F"; fi
# if [[ -z $GENEUMIL ]]; then  GENEUMIL="F"; fi
# if [[ -z $GENEUMIU ]]; then  GENEUMIU="F"; fi


FOUND_ERROR=0


# ===============================================
# CHECKING VARIABLES
# ===============================================
#check to ensure all mandatory arguments have been entered

if [ -z $OUTPUT_DIR ]; then echo "ERROR: missing mandatory option: -d (--dir) must be specified"; FOUND_ERROR=1; fi
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
          # cp $PIPELINE_HOME/scrna/pars/step2_par.txt $OUTPUT_DIR/job_info/parameters/
          # cp $PIPELINE_HOME/scrna/pars/step3_par.txt $OUTPUT_DIR/job_info/parameters/  
          # cp $PIPELINE_HOME/scrna/pars/step4_par.txt $OUTPUT_DIR/job_info/parameters/ 
          # cp $PIPELINE_HOME/scrna/pars/step5_par.txt $OUTPUT_DIR/job_info/parameters/ 
          # cp $PIPELINE_HOME/scrna/pars/step6_par.txt $OUTPUT_DIR/job_info/parameters/ 
          # cp $PIPELINE_HOME/scrna/pars/step7_par.txt $OUTPUT_DIR/job_info/parameters/ 
          # cp $PIPELINE_HOME/scrna/pars/step8_par.txt $OUTPUT_DIR/job_info/parameters/   
          # cp $PIPELINE_HOME/scrna/pars/step8_contrast_genotype.txt $OUTPUT_DIR/job_info/parameters/ 
          # cp $PIPELINE_HOME/scrna/pars/step8_contrast_celltype.txt $OUTPUT_DIR/job_info/parameters/ 
          # # cp $PIPELINE_HOME/scrna/pars/step8_clus_label.txt $OUTPUT_DIR/job_info/parameters/ 
          # cp $PIPELINE_HOME/scrna/pars/stepint_par.txt $OUTPUT_DIR/job_info/parameters/ 
        elif [[  ${SCRNA_METHOD0} =~ HTO ]] ; then
          cp $PIPELINE_HOME/hto/configs/scrnabox_config.ini $OUTPUT_DIR/job_info/configs/ 
          cp $PIPELINE_HOME/hto/pars/* $OUTPUT_DIR/job_info/parameters/ 
          # cp $PIPELINE_HOME/hto/pars/step2_par.txt $OUTPUT_DIR/job_info/parameters/
          # cp $PIPELINE_HOME/hto/pars/step3_par.txt $OUTPUT_DIR/job_info/parameters/  
          # cp $PIPELINE_HOME/hto/pars/step4_par.txt $OUTPUT_DIR/job_info/parameters/ 
          # # cp $PIPELINE_HOME/hto/pars/step4_antibody_label.txt $OUTPUT_DIR/job_info/parameters/ 
          # cp $PIPELINE_HOME/hto/pars/step5_par.txt $OUTPUT_DIR/job_info/parameters/ 
          # cp $PIPELINE_HOME/hto/pars/step6_par.txt $OUTPUT_DIR/job_info/parameters/ 
          # cp $PIPELINE_HOME/hto/pars/step7_par.txt $OUTPUT_DIR/job_info/parameters/ 
          # cp $PIPELINE_HOME/hto/pars/step8_par.txt $OUTPUT_DIR/job_info/parameters/   
          # cp $PIPELINE_HOME/hto/pars/step8_contrast_genotype.txt $OUTPUT_DIR/job_info/parameters/ 
          # cp $PIPELINE_HOME/hto/pars/step8_contrast_celltype.txt $OUTPUT_DIR/job_info/parameters/ 
          # cp $PIPELINE_HOME/hto/pars/step8_clus_label.txt $OUTPUT_DIR/job_info/parameters/ 
          # cp $PIPELINE_HOME/hto/pars/stepint_par.txt $OUTPUT_DIR/job_info/parameters/ 
        elif [[ -z  ${SCRNA_METHOD0}  ]]; then
          cp $PIPELINE_HOME/scrna/configs/scrnabox_config.ini $OUTPUT_DIR/job_info/configs/ 
          cp $PIPELINE_HOME/scrna/pars/* $OUTPUT_DIR/job_info/parameters/ 
          # cp $PIPELINE_HOME/scrna/pars/step2_par.txt $OUTPUT_DIR/job_info/parameters/
          # cp $PIPELINE_HOME/scrna/pars/step3_par.txt $OUTPUT_DIR/job_info/parameters/  
          # cp $PIPELINE_HOME/scrna/pars/step4_par.txt $OUTPUT_DIR/job_info/parameters/ 
          # cp $PIPELINE_HOME/scrna/pars/step5_par.txt $OUTPUT_DIR/job_info/parameters/ 
          # cp $PIPELINE_HOME/scrna/pars/step6_par.txt $OUTPUT_DIR/job_info/parameters/ 
          # cp $PIPELINE_HOME/scrna/pars/step7_par.txt $OUTPUT_DIR/job_info/parameters/ 
          # cp $PIPELINE_HOME/scrna/pars/step8_par.txt $OUTPUT_DIR/job_info/parameters/   
          # cp $PIPELINE_HOME/scrna/pars/step8_contrast_genotype.txt $OUTPUT_DIR/job_info/parameters/ 
          # cp $PIPELINE_HOME/scrna/pars/step8_contrast_celltype.txt $OUTPUT_DIR/job_info/parameters/ 
          # # cp $PIPELINE_HOME/scrna/pars/step8_clus_label.txt $OUTPUT_DIR/job_info/parameters/ 
          # cp $PIPELINE_HOME/scrna/pars/stepint_par.txt $OUTPUT_DIR/job_info/parameters/ 
              echo "NOTE: You didn't specify the method, so the pipeline selected"
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
          # cp $PIPELINE_HOME/scrna/pars/step2_par.txt $OUTPUT_DIR/job_info/parameters/
          # cp $PIPELINE_HOME/scrna/pars/step3_par.txt $OUTPUT_DIR/job_info/parameters/  
          # cp $PIPELINE_HOME/scrna/pars/step4_par.txt $OUTPUT_DIR/job_info/parameters/ 
          # cp $PIPELINE_HOME/scrna/pars/step5_par.txt $OUTPUT_DIR/job_info/parameters/ 
          # cp $PIPELINE_HOME/scrna/pars/step6_par.txt $OUTPUT_DIR/job_info/parameters/ 
          # cp $PIPELINE_HOME/scrna/pars/step7_par.txt $OUTPUT_DIR/job_info/parameters/ 
          # cp $PIPELINE_HOME/scrna/pars/step8_par.txt $OUTPUT_DIR/job_info/parameters/   
          # cp $PIPELINE_HOME/scrna/pars/step8_contrast_genotype.txt $OUTPUT_DIR/job_info/parameters/ 
          # cp $PIPELINE_HOME/scrna/pars/step8_contrast_celltype.txt $OUTPUT_DIR/job_info/parameters/ 
          # # cp $PIPELINE_HOME/scrna/pars/step8_clus_label.txt $OUTPUT_DIR/job_info/parameters/ 
          # cp $PIPELINE_HOME/scrna/pars/stepint_par.txt $OUTPUT_DIR/job_info/parameters/ 
        elif [[  ${SCRNA_METHOD0} =~ HTO ]] ; then
          cp $PIPELINE_HOME/hto/configs/scrnabox_config.ini $OUTPUT_DIR/job_info/configs/
          cp $PIPELINE_HOME/hto/pars/* $OUTPUT_DIR/job_info/parameters/ 
          # cp $PIPELINE_HOME/hto/pars/step2_par.txt $OUTPUT_DIR/job_info/parameters/
          # cp $PIPELINE_HOME/hto/pars/step3_par.txt $OUTPUT_DIR/job_info/parameters/  
          # cp $PIPELINE_HOME/hto/pars/step4_par.txt $OUTPUT_DIR/job_info/parameters/ 
          # # cp $PIPELINE_HOME/hto/pars/step4_antibody_label.txt $OUTPUT_DIR/job_info/parameters/ 
          # cp $PIPELINE_HOME/hto/pars/step5_par.txt $OUTPUT_DIR/job_info/parameters/ 
          # cp $PIPELINE_HOME/hto/pars/step6_par.txt $OUTPUT_DIR/job_info/parameters/ 
          # cp $PIPELINE_HOME/hto/pars/step7_par.txt $OUTPUT_DIR/job_info/parameters/ 
          # cp $PIPELINE_HOME/hto/pars/step8_par.txt $OUTPUT_DIR/job_info/parameters/   
          # cp $PIPELINE_HOME/hto/pars/step8_contrast_genotype.txt $OUTPUT_DIR/job_info/parameters/ 
          # cp $PIPELINE_HOME/hto/pars/step8_contrast_celltype.txt $OUTPUT_DIR/job_info/parameters/ 
          # # cp $PIPELINE_HOME/hto/pars/step8_clus_label.txt $OUTPUT_DIR/job_info/parameters/ 
          # cp $PIPELINE_HOME/hto/pars/stepint_par.txt $OUTPUT_DIR/job_info/parameters/ 
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
        echo -e  "--------Pipeline is set up-----------------" $VERSION >> $EXPECTED_DONE_FILES
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
export MODULEUSE=$MODULEUSE
export MODULECELLRANGER=$MODULECELLRANGER
export ACCOUNT=$ACCOUNT
export OUTPUT_DIR=$OUTPUT_DIR

echo -e "NOTE: the pipeline is running  << $SCRNA_METHOD seq >>"


#set queue-specific values, depending on system check
if [[ $submit_cmd =~ sbatch ]]; then
   QUEUE="sbatch"          # default job scheduler: qsub
elif [[ $submit_cmd =~ bash ]]; then
    QUEUE=$PIPELINE_HOME/soft/queueInterpreter/QueueInterpreter
else
    echo "The pipeline is testing under slurm system and ubuntu,"
    echo "please choose sbatch for slurm, or bash for ubuntu"   
    exit 42
fi


# ===============================================
# STEP 1: RUN cellranger 
# ===============================================
#

TEMPCONFIG=$OUTPUT_DIR/job_info/.tmp/temp_config.ini
touch $TEMPCONFIG
if [[ -f "TEMPCONFIG" ]]; then 
  rm $TEMPCONFIG
fi

echo " # IT IS A temp FILE. DO NOT EDIT THIS FILE DIRECTLY."  > $TEMPCONFIG
echo STEP8m=$STEP8m  >> $TEMPCONFIG
echo STEP8dgelist=$STEP8dgelist   >> $TEMPCONFIG
echo STEP8i=$STEP8i  >> $TEMPCONFIG
echo MSD=$MSD  >> $TEMPCONFIG
echo CONT=$CONT  >> $TEMPCONFIG
# echo NFRNAL=$NFRNAL  >> $TEMPCONFIG
# echo NFRNAU=$NFRNAU  >> $TEMPCONFIG
# echo NCRNAL=$NCRNAL  >> $TEMPCONFIG
# echo NCRNAU=$NCRNAU  >> $TEMPCONFIG
# echo PMTL=$PMTL  >> $TEMPCONFIG
# echo PMTU=$PMTU  >> $TEMPCONFIG
# echo GENEUMIL=$GENEUMIL >> $TEMPCONFIG
# echo GENEUMIU=$GENEUMIU >> $TEMPCONFIG
echo SEULIST=$SEULIST   >> $TEMPCONFIG
echo STEP7enrich=$STEP7enrich   >> $TEMPCONFIG
echo STEP7fta=$STEP7fta  >> $TEMPCONFIG
echo STEP7marker=$STEP7marker  >> $TEMPCONFIG
echo SINFO=$SINFO  >> $TEMPCONFIG
echo OUTPUT_DIR=$OUTPUT_DIR >> $TEMPCONFIG
echo JOB_OUTPUT_DIR=$JOB_OUTPUT_DIR  >> $TEMPCONFIG
echo EXPECTED_DONE_FILES=$EXPECTED_DONE_FILES   >> $TEMPCONFIG
echo MODE=$MODE >> $TEMPCONFIG
echo QUEUE=$QUEUE >> $TEMPCONFIG
echo VERSION=$VERSION >> $TEMPCONFIG



# $OUTPUT_DIR/job_info/configs/temp
if  [[  ${MODE0[@]}  =~  integrate  ]]  ; then
   bash ${PIPELINE_HOME}/launch/launch_scrnabox_integ.sh
elif  [[  ${SCRNA_METHOD} =~ SCRNA ]] ; then
   bash ${PIPELINE_HOME}/launch/launch_scrnabox_scrna.sh
elif [[  ${SCRNA_METHOD} =~ HTO ]] ; then
   bash ${PIPELINE_HOME}/launch/launch_scrnabox_hto.sh
fi

exit 0
