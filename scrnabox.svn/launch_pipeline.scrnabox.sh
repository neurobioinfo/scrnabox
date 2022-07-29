#!/bin/bash

# The pipeline is done as part of Project ....
# Copyright belong .......  
# Written by Saeid Amiri (saeid.amiri@mcgill.ca) with associate with Rhalena Thomas, Roxanne...  
VERSION=0.1.0; echo " ongoing scrna pipeline version $VERSION"
#last updated version 2022-07-27
# test on ...

# ===============================================
# default variables values
# ===============================================
unset SAMPLE OUTPUT_DIR PIPELINE_HOME QUEUE ACCOUNT STEP7m STEP7i

#Assuming script is in root of PIPELINE_HOME
export PIPELINE_HOME=$(cd $(dirname $0) && pwd -P)

# export CONFIG_FILES=$PIPELINE_HOME/configs/scrnabox.config.ini
TIMESTAMP=`date +%FT%H.%M.%S`

#set queue-specific values, depending on system check
QUEUE="sbatch"                                                                                # default job scheduler: qsub

# create function to handle error messages
# ===============================================
Usage() {
	echo
	echo -e "Usage:\t$0 [arguments]"
	echo -e "\tmandatory arguments:\n" \
          "\t\t-s  (--sample)            = folder includes sample information (give full path) \n" \
          "\t\t-d  (--dir)               = run directory (where all the outputs will be printed) (give full path)\n" \
          "\t\t--steps                   = 'ALL' to run all steps, or specify what steps, e.g., 2 to run just step 2, 2-4, run steps 2 through 4)\n" 
	echo -e "\toptional arguments:\n " \
          "\t\t-h  (--help)              = get the program options and exit\n" \
          "\t\t--nFRNAl                = set account on cluster (to be used with qsub with \n" \
          "\t\t--nFRNAu                = set account on cluster (to be used with qsub with \n" \
          "\t\t--main                = set account on cluster (to be used with qsub with \n" \
          "\t\t--inte                = set account on cluster (to be used with qsub with \n" \
          "\t\t--pmt                = set account on cluster (to be used with qsub with \n" \
          "\t\t--msd                = set account on cluster (to be used with qsub with \n" \
          "\t\t--cont                = set account on cluster (to be used with qsub with \n" \
          "\t\t--account                = set account on cluster (to be used with qsub with -A [CURRENT \"$ACCOUNT\"]\n" \
          "\t\t-m  (--max_mem)           = max memory on a node (in gigs) to use [CURRENT \"$MAX_MEM\"]\n" \
          "\t\t-t  (--threads)           = number of threads to use [CURRENT \"$THREADS\"]\n" \
          "\t\t-w  (--walltime)          = wall time for the scheduler (format HH:MM:SS) [CURRENT \"$WALLTIME\"]\n" \
          "\t\t-v  (--verbose)           = set verbosity level [CURRENT \"$VERBOSE\"]\n" 
        echo
}

# nFRNAl > 300 & nFRNAu < 6500 & pmt


# ===============================================
# step 0 create folder 
# ===============================================

# ===============================================
# PARSING ARGUMENTS
# ===============================================
# options may be followed by one colon to indicate they have a required argument
if ! options=$(getopt --name pipeline --alternative --unquoted --options hs:d:t:m:vw:f:S:c:a:x: --longoptions sample:,dir:,steps:,nFRNAl:,nFRNAu:,main:,inte:,pmt:,msd:,cont:,verbose,help -- "$@")
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

# ===============================================
# LOAD ALL OTHER OPTIONS
# ===============================================
set -- $options

while [ $# -gt 0 ]
#echo -e "$1 $2"
do
    case $1 in
    -h| --help) Usage; exit 0;;
    # for options with required arguments, an additional shift is required
    -s| --sample) SAMPLE="$2" ; shift ;;
    -d| --dir) OUTPUT_DIR="$2" ; shift ;;
    -v| --verbose) VERBOSE=1 ;; 
    --steps) MODE="$2"; shift ;; #SA
    --nFRNAl) NFRNAL="$2"; shift ;; #SA
    --nFRNAu) NFRNAU="$2"; shift ;; #SA
    --pmt) PMT="$2"; shift ;; #SA
    --main) STEP7m="$2"; shift ;; #SA
    --inte) STEP7i="$2"; shift ;; #SA
    --msd) MSD="$2"; shift ;; #SA
    --cont) CONT="$2"; shift ;; #SA
    (--) shift; break;;
    (-*) echo "$0: error - unrecognized option $1" 1>&2; exit 42;;
    (*) break;;
    esac
    shift
done

if [[ -z $STEP7m ]]; then  STEP7m="F"; fi
if [[ -z $STEP7i ]]; then  STEP7i="F"; fi
if [[ -z $MSD ]]; then  MSD="F"; fi
if [[ -z $CONT ]]; then  CONT="F"; fi

# echo $CONT


if [[ "$MODE" == *"-"* ]]; then
  abc0=`echo $MODE | sed  "s/-/../g"  | awk '{print "{"$0"}"}'`
  MODE0=`eval echo $abc0`
else 
  MODE0=$MODE
fi

FOUND_ERROR=0

# ===============================================
# CHECKING VARIABLES
# ===============================================
#check to ensure all mandatory arguments have been entered

# if [ -z $SAMPLE ]; then echo "ERROR: missing mandatory option: -s (sample folder) must be specified"; FOUND_ERROR=1; fi
if [ -z $OUTPUT_DIR ]; then echo "ERROR: missing mandatory option: -d (--dir) must be specified"; FOUND_ERROR=1; fi
if [ -z $MODE ]; then echo "ERROR: missing mandatory option: --steps ('ALL' to run all, 2 to run step 2, step 2-4, run steps 2 through 4) must be specified"; FOUND_ERROR=1; fi

if (( $FOUND_ERROR )); then echo "Please check options and try again"; exit 42; fi
if [ $MODE == 'ALL' ]; then  MODE0=`echo {2..10}`; fi 


# STEP 0: RUN setting 
# ===============================================
#

if [[ ${MODE0[@]} == 0 ]]; then 
  # if  [ ! -d $OUTPUT_DIR ]; then   mkdir -p $OUTPUT_DIR  ; fi 
  JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
  if [ ! -d $JOB_OUTPUT_DIR ]; then 
    mkdir -p $OUTPUT_DIR/job_output
    mkdir -p $OUTPUT_DIR/job_output/temp
  fi
  # EXPECTED_DONE_FILES=$JOB_OUTPUT_DIR/expected.done.files.txt
  # if [ ! -z $EXPECTED_DONE_FILES ]; then 
    # touch $JOB_OUTPUT_DIR/expected.done.files.txt
  # fi
  if [ ! -z $EXPECTED_DONE_FILES ]; then 
    touch $JOB_OUTPUT_DIR/expected.done.files.txt
  else
    rm -f  $JOB_OUTPUT_DIR/expected.done.files.txt
    touch $JOB_OUTPUT_DIR/expected.done.files.txt
  fi
  if [ ! -f $OUTPUT_DIR/job_output/scrnabox.config.ini ]; then 
     cp $PIPELINE_HOME/configs/scrnabox.config.ini $OUTPUT_DIR/job_output/ 
  else
    rm -f $OUTPUT_DIR/job_output/scrnabox.config.ini
    cp $PIPELINE_HOME/configs/scrnabox.config.ini $OUTPUT_DIR/job_output/ 
  fi
  if [ ! -f $OUTPUT_DIR/job_output/step4_par.txt ]; then 
     cp $PIPELINE_HOME/configs/step4_par.txt $OUTPUT_DIR/job_output/ 
  else
    rm -f $OUTPUT_DIR/job_output/step4_par.txt
    cp $PIPELINE_HOME/configs/step4_par.txt $OUTPUT_DIR/job_output/ 
  fi
  if [ ! -f $OUTPUT_DIR/job_output/step7_contrast_main.txt ]; then 
     cp $PIPELINE_HOME/configs/step7_contrast_main.txt $OUTPUT_DIR/job_output/ 
  else
    rm -f $OUTPUT_DIR/job_output/step7_contrast_main.txt
    cp $PIPELINE_HOME/configs/step7_contrast_main.txt $OUTPUT_DIR/job_output/ 
  fi
  if [ ! -f $OUTPUT_DIR/job_output/step7_contrast_inte.txt ]; then 
     cp $PIPELINE_HOME/configs/step7_contrast_inte.txt $OUTPUT_DIR/job_output/ 
  else
    rm -f $OUTPUT_DIR/job_output/step7_contrast_inte.txt
    cp $PIPELINE_HOME/configs/step7_contrast_inte.txt $OUTPUT_DIR/job_output/ 
  fi
  if [ ! -f $OUTPUT_DIR/job_output/step7_clus_label.txt ]; then 
     cp $PIPELINE_HOME/configs/step7_clus_label.txt $OUTPUT_DIR/job_output/ 
  else
    rm -f $OUTPUT_DIR/job_output/step7_clus_label.txt
    cp $PIPELINE_HOME/configs/step7_clus_label.txt $OUTPUT_DIR/job_output/ 
  fi
  # export EXPECTED_DONE_FILES=$JOB_OUTPUT_DIR/expected.done.files.txt
  # if  [ ! -f $OUTPUT_DIR/job_output/scrnabox.config.ini ]; then   cp $PIPELINE_HOME/configs/scrnabox.config.ini $OUTPUT_DIR/job_output/ ; fi 
  if  [ ! -d $OUTPUT_DIR/samples_info ]; then   echo 'ERROR: The pipline can not findsamples_info directory,  '; exit 42 ; fi 
fi 


# source $PIPELINE_HOME/configs/scrnabox.config.ini
source $OUTPUT_DIR/job_output/scrnabox.config.ini
export JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
export EXPECTED_DONE_FILES=$JOB_OUTPUT_DIR/expected.done.files.txt
export STEP4_PAR=$OUTPUT_DIR/job_output/step4_par.txt
export STEP7_CON=$PIPELINE_HOME/configs/step7_contrast.txt
chmod 775 $EXPECTED_DONE_FILES    


# ===============================================
# STEP 1: RUN cellranger 

# ===============================================
#
STEP=step_1

if [[  ${MODE0[@]}  =~  1  ]]; then
  if [ ! -d $OUTPUT_DIR/step1 ]; then 
    mkdir -p $OUTPUT_DIR/step1 
    # if [ -d $SAMPLE ]; then
      # cp -r  $SAMPLE/* $OUTPUT_DIR/step1 
    # else 
    cp -r  $OUTPUT_DIR/samples_info/* $OUTPUT_DIR/step1 
    # fi
  fi
  if  [ -f $JOB_OUTPUT_DIR/sample.list ]; then rm $JOB_OUTPUT_DIR/sample.list; touch $JOB_OUTPUT_DIR/sample.list ; else touch $JOB_OUTPUT_DIR/sample.list; fi 
  if  [ -f $JOB_OUTPUT_DIR/sample_dir.list ]; then rm $JOB_OUTPUT_DIR/sample_dir.list; touch $JOB_OUTPUT_DIR/sample_dir.list ; else touch $JOB_OUTPUT_DIR/sample_dir.list; fi  
  search_dir=${OUTPUT_DIR}/step1
  for entry in "$search_dir"/*
  do
    echo "$entry" >> ${OUTPUT_DIR}/job_output/sample_dir.list
    echo $(basename  "$entry") >> ${OUTPUT_DIR}/job_output/sample.list
  done
  while read item
  do
    cp ${PIPELINE_HOME}/scripts/step1/slurm.template $item
    cp ${PIPELINE_HOME}/scripts/step1/launch_cellranger.hashtagged.slurm.sh $item
  done < ${OUTPUT_DIR}/job_output/sample_dir.list
  while read item
  do
    cd ${item}; sh  ${item}/launch_cellranger.hashtagged.slurm.sh -r ouput_folder &
  done < ${OUTPUT_DIR}/job_output/sample_dir.list
  wait
  echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  echo "--------"  >> $EXPECTED_DONE_FILES
  echo "step 1 submited " >> $EXPECTED_DONE_FILES
fi 


if  [ ! -f $JOB_OUTPUT_DIR/sample.list ]; then
  if  [ ! -f $JOB_OUTPUT_DIR/sample.list ]; then
    search_dir=${OUTPUT_DIR}/step1
    for entry in "$search_dir"/*
    do
      echo "$entry" >> ${OUTPUT_DIR}/job_output/sample_dir.list
      echo $(basename  "$entry") >> ${OUTPUT_DIR}/job_output/sample.list
    done
  fi 
fi 

SAMPLE_SIZE=`wc -l < ${OUTPUT_DIR}/job_output/sample.list`
# ===============================================
# STEP 2: Generate the Seurat onject
# ===============================================
#
STEP=step_2
export THREADS=$((SAMPLE_SIZE*2)) #${THREADS_ARRAY[$STEP]}
export WALLTIME=${WALLTIME_ARRAY[$STEP]}
export MEM=${MEM_ARRAY[$STEP]}
if [[  ${MODE0[@]}  =~  2 ]]  &&  [[  ${MODE0[@]} =~ 1 ]] ; then
  if [ ! -d $OUTPUT_DIR/step2 ]; then 
    mkdir -p $OUTPUT_DIR/step2/objs 
  fi
  echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  echo "--------"  >> $EXPECTED_DONE_FILES
  echo "STEP 2 submitted following step 1"  >> $EXPECTED_DONE_FILES
  step_2="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem=${MEM} \
    --time=${WALLTIME} \
    --job-name $STEP \
    $DEPEND_CELLRANGER \
    --export OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH} \
    --output $JOB_OUTPUT_DIR/temp/%x.o%j \
    $PIPELINE_HOME/scripts/step2/pipeline.step2.qsub"
  step_2=$($step_2 | grep -oP "\d+")
  echo "[Q] STEP 2         : $step_2 " >> $EXPECTED_DONE_FILES
  DEPEND_step_2="--dependency=afterok:$step_2"; echo $DEPEND_step_2
elif [[  ${MODE0[@]}  =~  2  ]]  &&  [[  ${MODE0[@]} != 1 ]]; then
  if [ ! -d $OUTPUT_DIR/step2 ]; then 
    mkdir -p $OUTPUT_DIR/step2/objs 
  fi
  echo "just step 2 at" $TIMESTAMP >> $EXPECTED_DONE_FILES
  echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  echo "--------"  >> $EXPECTED_DONE_FILES
  echo "STEP 2 submitted (without following  step 1) "  >> $EXPECTED_DONE_FILES
  step_2="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem=${MEM} \
    --time=${WALLTIME} \
    --job-name $STEP \
    --export OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH} \
    --output $JOB_OUTPUT_DIR/temp/%x.o%j \
    $PIPELINE_HOME/scripts/step2/pipeline.step2.qsub"
  step_2=$($step_2 | grep -oP "\d+")
  echo "[Q] STEP 2         : $step_2 " >> $EXPECTED_DONE_FILES 
  DEPEND_step_2="--dependency=afterok:$step_2"; echo $DEPEND_step_2
fi 


# ===============================================
# STEP 3: 
# ===============================================
#
STEP=step_3
export THREADS=$((SAMPLE_SIZE*2)) #${THREADS_ARRAY[$STEP]}
export WALLTIME=${WALLTIME_ARRAY[$STEP]}
export MEM=${MEM_ARRAY[$STEP]}
if [[  ${MODE0[@]}  =~  3 ]]  &&  [[  ${MODE0[@]} =~ 2 ]] ; then
  if [ ! -d $OUTPUT_DIR/step3 ]; then 
    mkdir -p $OUTPUT_DIR/step3/objs 
    # mkdir -p $OUTPUT_DIR/step3/figs     
  fi
  echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  echo "--------"  >> $EXPECTED_DONE_FILES
  echo "STEP 3 submitted following step 2"  >> $EXPECTED_DONE_FILES
  step_3="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem=${MEM} \
    --time=${WALLTIME} \
    --job-name $STEP \
    $DEPEND_step_2 \
    --export OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},NFRNAL=${NFRNAL},NFRNAU=${NFRNAU},PMT=${PMT} \
    --output $JOB_OUTPUT_DIR/temp/%x.o%j \
    $PIPELINE_HOME/scripts/step3/pipeline.step3.qsub"
  step_3=$($step_3 | grep -oP "\d+")
  echo "[Q] STEP 3         : $step_3 " >> $EXPECTED_DONE_FILES
  DEPEND_step_3="--dependency=afterok:$step_3"; echo $DEPEND_step_3
elif [[  ${MODE0[@]}  =~  3  ]]  &&  [[  ${MODE0[@]} != 2 ]]; then
  if [ ! -d $OUTPUT_DIR/step3 ]; then 
    mkdir -p $OUTPUT_DIR/step3/objs 
    mkdir -p $OUTPUT_DIR/step3/figs 
  fi
  echo "just step 3 at" $TIMESTAMP >> $EXPECTED_DONE_FILES
  echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  echo "--------"  >> $EXPECTED_DONE_FILES
  echo "STEP 3 submitted (without following  step 2) "  >> $EXPECTED_DONE_FILES
  step_3="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem=${MEM} \
    --time=${WALLTIME} \
    --job-name $STEP \
    --export OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},NFRNAL=${NFRNAL},NFRNAU=${NFRNAU},PMT=${PMT} \
    --output $JOB_OUTPUT_DIR/temp/%x.o%j \
    $PIPELINE_HOME/scripts/step3/pipeline.step3.qsub"
  step_3=$($step_3 | grep -oP "\d+")
  echo "[Q] STEP 3         : $step_3 " >> $EXPECTED_DONE_FILES 
  DEPEND_step_3="--dependency=afterok:$step_3"; echo $DEPEND_step_3
fi 

# ===============================================
# STEP 4: 
# ===============================================
#
STEP=step_4
export THREADS=$((SAMPLE_SIZE*2)) #${THREADS_ARRAY[$STEP]}
export WALLTIME=${WALLTIME_ARRAY[$STEP]}
export MEM=${MEM_ARRAY[$STEP]}
if [[  ${MODE0[@]}  =~  4 ]]  &&  [[  ${MODE0[@]} =~ 3 ]] && [[  ${MSD}  =~  F  ]] ; then
  if [ ! -d $OUTPUT_DIR/step4 ]; then 
    mkdir -p $OUTPUT_DIR/step4/objs 
    mkdir -p $OUTPUT_DIR/step4/figs     
  fi
  echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  echo "--------"  >> $EXPECTED_DONE_FILES
  echo "STEP 4 submitted following step 3"  >> $EXPECTED_DONE_FILES
  step_4="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem=${MEM} \
    --time=${WALLTIME} \
    --job-name $STEP \
    $DEPEND_step_3 \
    --export=OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},STEP4_PAR=${STEP4_PAR} \
    --output $JOB_OUTPUT_DIR/temp/%x.o%j \
    $PIPELINE_HOME/scripts/step4/pipeline.step4.qsub"
  step_4=$($step_4 | grep -oP "\d+")
  echo "[Q] STEP 4         : $step_4 " >> $EXPECTED_DONE_FILES
  DEPEND_step_4="--dependency=afterok:$step_4"; echo $DEPEND_step_4
elif [[  ${MODE0[@]}  =~  4  ]]  &&  [[  ${MODE0[@]} != 3 ]] && [[  ${MSD}  =~  F  ]]; then
  if [ ! -d $OUTPUT_DIR/step4 ]; then 
    mkdir -p $OUTPUT_DIR/step4/objs 
    mkdir -p $OUTPUT_DIR/step4/figs 
  fi
  echo "just step 4 at" $TIMESTAMP >> $EXPECTED_DONE_FILES
  echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  echo "--------"  >> $EXPECTED_DONE_FILES
  echo "STEP 4 submitted (without following  step 3) "  >> $EXPECTED_DONE_FILES
  step_4="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem=${MEM} \
    --time=${WALLTIME} \
    --job-name $STEP \
    --export=OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},STEP4_PAR=${STEP4_PAR} \
    --output $JOB_OUTPUT_DIR/temp/%x.o%j \
    $PIPELINE_HOME/scripts/step4/pipeline.step4.qsub"
  step_4=$($step_4 | grep -oP "\d+")
  echo "[Q] STEP 4         : $step_4 " >> $EXPECTED_DONE_FILES 
  DEPEND_step_4="--dependency=afterok:$step_4"; echo $DEPEND_step_4
fi 

if [[  ${MSD}  =~  T  ]] ; then
  if [ ! -d $OUTPUT_DIR/step4 ]; then 
    mkdir -p $OUTPUT_DIR/step4/objs 
    mkdir -p $OUTPUT_DIR/step4/figs 
  fi
  echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  echo "--------"  >> $EXPECTED_DONE_FILES
  echo "STEP 4-MSD "  >> $EXPECTED_DONE_FILES
  step_4="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem-per-cpu=6g \
    --time=${WALLTIME} \
    --job-name ${STEP}_MSD \
    --export=OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},STEP4_PAR=${STEP4_PAR} \
    --output $JOB_OUTPUT_DIR/temp/%x.o%j \
    $PIPELINE_HOME/scripts/step4/pipeline.step4_msd.qsub"
  step_4=$($step_4 | grep -oP "\d+")
  echo "[Q] STEP 4         : $step_4 " >> $EXPECTED_DONE_FILES 
  DEPEND_step_4="--dependency=afterok:$step_4"; echo $DEPEND_step_4
fi 

# ===============================================
# STEP 5: 
# ===============================================
#
STEP=step_5
export THREADS=$((SAMPLE_SIZE*2)) #${THREADS_ARRAY[$STEP]}
export WALLTIME=${WALLTIME_ARRAY[$STEP]}
export MEM=${MEM_ARRAY[$STEP]}
if [[  ${MODE0[@]}  =~  5 ]]  &&  [[  ${MODE0[@]} =~ 4 ]] ; then
  if [ ! -d $OUTPUT_DIR/step5 ]; then 
    mkdir -p $OUTPUT_DIR/step5/objs    
  fi
  echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  echo "--------"  >> $EXPECTED_DONE_FILES
  echo "STEP 5 submitted following step 4"  >> $EXPECTED_DONE_FILES
  step_5="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem=${MEM} \
    --time=${WALLTIME} \
    --job-name $STEP \
    $DEPEND_step_4 \
    --export=OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH} \
    --output $JOB_OUTPUT_DIR/temp/%x.o%j \
    $PIPELINE_HOME/scripts/step5/pipeline.step5.qsub"
  step_5=$($step_5 | grep -oP "\d+")
  echo "[Q] STEP 5         : $step_5 " >> $EXPECTED_DONE_FILES
  DEPEND_step_5="--dependency=afterok:$step_5"; echo $DEPEND_step_5
elif [[  ${MODE0[@]}  =~  5  ]]  &&  [[  ${MODE0[@]} != 4 ]]; then
  if [ ! -d $OUTPUT_DIR/step5 ]; then 
    mkdir -p $OUTPUT_DIR/step5/objs 
  fi
  echo "just step 5 at" $TIMESTAMP >> $EXPECTED_DONE_FILES
  echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  echo "--------"  >> $EXPECTED_DONE_FILES
  echo "STEP 5 submitted (without following  step 4) "  >> $EXPECTED_DONE_FILES
  step_5="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem=${MEM} \
    --time=${WALLTIME} \
    --job-name $STEP \
    --export=OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH} \
    --output $JOB_OUTPUT_DIR/temp/%x.o%j \
    $PIPELINE_HOME/scripts/step5/pipeline.step5.qsub"
  step_5=$($step_5 | grep -oP "\d+")
  echo "[Q] STEP 5         : $step_5 " >> $EXPECTED_DONE_FILES 
  DEPEND_step_5="--dependency=afterok:$step_5"; echo $DEPEND_step_5
fi 






# ===============================================
# STEP 6: 
# ===============================================
#
STEP=step_6
export THREADS=$((SAMPLE_SIZE*2)) #${THREADS_ARRAY[$STEP]}
export WALLTIME=${WALLTIME_ARRAY[$STEP]}
export MEM=${MEM_ARRAY[$STEP]}
if [[  ${MODE0[@]}  =~  6 ]]  &&  [[  ${MODE0[@]} =~ 5 ]] ; then
  if [ ! -d $OUTPUT_DIR/step6 ]; then 
    mkdir -p $OUTPUT_DIR/step6/objs    
    mkdir -p $OUTPUT_DIR/step6/figs    
  fi
  echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  echo "--------"  >> $EXPECTED_DONE_FILES
  echo "STEP 6 submitted following step 5"  >> $EXPECTED_DONE_FILES
  step_6="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem=${MEM} \
    --time=${WALLTIME} \
    --job-name $STEP \
    $DEPEND_step_5 \
    --export=OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH} \
    --output $JOB_OUTPUT_DIR/temp/%x.o%j \
    $PIPELINE_HOME/scripts/step6/pipeline.step6.qsub"
  step_6=$($step_6 | grep -oP "\d+")
  echo "[Q] STEP 6         : $step_6 " >> $EXPECTED_DONE_FILES
  DEPEND_step_6="--dependency=afterok:$step_6"; echo $DEPEND_step_6
elif [[  ${MODE0[@]}  =~  6  ]]  &&  [[  ${MODE0[@]} != 5 ]]; then
  if [ ! -d $OUTPUT_DIR/step6 ]; then 
    mkdir -p $OUTPUT_DIR/step6/objs 
    mkdir -p $OUTPUT_DIR/step6/figs 
  fi
  echo "just step 6 at" $TIMESTAMP >> $EXPECTED_DONE_FILES
  echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  echo "--------"  >> $EXPECTED_DONE_FILES
  echo "STEP 6 submitted (without following  step 5) "  >> $EXPECTED_DONE_FILES
  step_6="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem=${MEM} \
    --time=${WALLTIME} \
    --job-name $STEP \
    --export=OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH} \
    --output $JOB_OUTPUT_DIR/temp/%x.o%j \
    $PIPELINE_HOME/scripts/step6/pipeline.step6.qsub"
  step_6=$($step_6 | grep -oP "\d+")
  echo "[Q] STEP 6         : $step_6 " >> $EXPECTED_DONE_FILES 
  DEPEND_step_6="--dependency=afterok:$step_6"; echo $DEPEND_step_6
fi 



# echo ${STEP7i}
# echo ${STEP7m}

# ===============================================
# STEP 7: 
# ===============================================
#
STEP=step_7
THRm=`wc -l < ${OUTPUT_DIR}/job_output/step7_contrast_main.txt`
THRi=`wc -l < ${OUTPUT_DIR}/job_output/step7_contrast_inte.txt`
export THREADS=$((SAMPLE_SIZE*3)) #${THREADS_ARRAY[$STEP]}
export WALLTIME=${WALLTIME_ARRAY[$STEP]}
export MEM=${MEM_ARRAY[$STEP]}
if [[  ${MODE0[@]}  =~  7 ]]  &&  [[  ${MODE0[@]} =~ 6 ]] &&  [[   ${STEP7m}  =~  T ]]; then
  if [ ! -d $OUTPUT_DIR/step7 ]; then 
    # mkdir -p $OUTPUT_DIR/step7/objs    
    # mkdir -p $OUTPUT_DIR/step7/figs    
    mkdir -p $OUTPUT_DIR/step7/cont_main
    mkdir -p $OUTPUT_DIR/step7/cont_inte
  fi
  echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  echo "--------"  >> $EXPECTED_DONE_FILES
  echo "STEP 7 submitted following step 6"  >> $EXPECTED_DONE_FILES
  step_6="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem=${MEM}  \
    --time=${WALLTIME} \
    --job-name ${STEP}_main \
    $DEPEND_step_6 \
    --export OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},CONT=${CONT} \
    --output $JOB_OUTPUT_DIR/temp/%x.o%j \
    $PIPELINE_HOME/scripts/step7/pipeline.step7m.qsub"
  step_7=$($step_7 | grep -oP "\d+")
  echo "[Q] STEP 7         : $step_7 " >> $EXPECTED_DONE_FILES
  DEPEND_step_7="--dependency=afterok:$step_7"; echo $DEPEND_step_7
elif [[  ${MODE0[@]}  =~  7  ]]  &&  [[  ${MODE0[@]} != 6 ]] &&  [[   ${STEP7m}  =~  T ]]; then
  if [ ! -d $OUTPUT_DIR/step7 ]; then 
    # mkdir -p $OUTPUT_DIR/step7/objs 
    # mkdir -p $OUTPUT_DIR/step7/figs 
    mkdir -p $OUTPUT_DIR/step7/cont_main
    mkdir -p $OUTPUT_DIR/step7/cont_inte
  fi
  echo "just step 7 at" $TIMESTAMP >> $EXPECTED_DONE_FILES
  echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  echo "--------"  >> $EXPECTED_DONE_FILES
  echo "STEP 7 submitted (without following  step 6) "  >> $EXPECTED_DONE_FILES
  step_7="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem=${MEM} \
    --time=${WALLTIME} \
    --job-name ${STEP}_main \
    --export OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},CONT=${CONT} \
    --output $JOB_OUTPUT_DIR/temp/%x.o%j \
    $PIPELINE_HOME/scripts/step7/pipeline.step7m.qsub"
  step_7=$($step_7 | grep -oP "\d+")
  echo "[Q] STEP 7         : $step_7 " >> $EXPECTED_DONE_FILES 
  DEPEND_step_7="--dependency=afterok:$step_7"; echo $DEPEND_step_7
fi 
if [[  ${MODE0[@]}  =~  7 ]]  &&  [[  ${MODE0[@]} =~ 6 ]] &&  [[   ${STEP7i}  =~  T ]]; then
  if [ ! -d $OUTPUT_DIR/step7 ]; then 
    # mkdir -p $OUTPUT_DIR/step7/objs    
    # mkdir -p $OUTPUT_DIR/step7/figs   
    mkdir -p $OUTPUT_DIR/step7/cont_main
    mkdir -p $OUTPUT_DIR/step7/cont_inte 
  fi
  echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  echo "--------"  >> $EXPECTED_DONE_FILES
  echo "STEP 7 submitted following step 6"  >> $EXPECTED_DONE_FILES
  step_6="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem=${MEM} \
    --time=${WALLTIME} \
    --job-name ${STEP}_inte \
    $DEPEND_step_6 \
    --export OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},CONT=${CONT} \
    --output $JOB_OUTPUT_DIR/temp/%x.o%j \
    $PIPELINE_HOME/scripts/step7/pipeline.step7i.qsub"
  step_7=$($step_7 | grep -oP "\d+")
  echo "[Q] STEP 7         : $step_7 " >> $EXPECTED_DONE_FILES
  DEPEND_step_7="--dependency=afterok:$step_7"; echo $DEPEND_step_7
elif [[  ${MODE0[@]}  =~  7  ]]  &&  [[  ${MODE0[@]} != 6 ]] &&  [[   ${STEP7i}  =~  T ]]; then
  if [ ! -d $OUTPUT_DIR/step7 ]; then 
    # mkdir -p $OUTPUT_DIR/step7/objs 
    # mkdir -p $OUTPUT_DIR/step7/figs 
    mkdir -p $OUTPUT_DIR/step7/cont_main
    mkdir -p $OUTPUT_DIR/step7/cont_inte
  fi
  echo "just step 7 at" $TIMESTAMP >> $EXPECTED_DONE_FILES
  echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  echo "--------"  >> $EXPECTED_DONE_FILES
  echo "STEP 7 submitted (without following  step 6) "  >> $EXPECTED_DONE_FILES
  step_7="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem=${MEM} \
    --time=${WALLTIME} \
    --job-name ${STEP}_inte \
    --export OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},CONT=${CONT} \
    --output $JOB_OUTPUT_DIR/temp/%x.o%j \
    $PIPELINE_HOME/scripts/step7/pipeline.step7i.qsub"
  step_7=$($step_7 | grep -oP "\d+")
  echo "[Q] STEP 7         : $step_7 " >> $EXPECTED_DONE_FILES 
  DEPEND_step_7="--dependency=afterok:$step_7"; echo $DEPEND_step_7
fi 




exit 0

