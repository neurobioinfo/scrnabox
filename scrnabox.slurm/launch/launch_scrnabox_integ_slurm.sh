#!/bin/bash


# The pipeline is done as part of Project Dark Genome 
# Copyright belong MNI BIOINFO CORE (https://github.com/neurobioinfo)
# The pipeline is written by Saeid Amiri (saeid.amiri@mcgill.ca)

declare -A THREADS_ARRAY
declare -A  WALLTIME_ARRAY
declare -A  MEM_ARRAY
source $OUTPUT_DIR/job_info/configs/scrnabox_config.ini
source $OUTPUT_DIR/job_info/.tmp/temp_config.ini
source $PIPELINE_HOME/tools/utils.sh
# ===============================================
# Integrate: Add the Seurat object
# ===============================================
#
STEP=integrate 

TOTALSIZE=0
while read -r item
do
  if [ -z "$item" ]
  then
    continue
  fi 
  temp=`stat -c%s  $item`
  TOTALSIZE=$(($TOTALSIZE + $temp))
done < $SEULIST

y=10**9
TOTALSIZE=` printf "%.0f\n" $((10**6 * TOTALSIZE/y))e-6`
if (( TOTALSIZE == 0 )); then
  TOTALSIZE=1
fi

# if [[ $QUEUE =~ sbatch ]] && [[  ${MODE0[@]}  =~  integrate  ]] ; then
  if [ -z "$SAMPLE_SIZE" ]; then
    SAMPLE_SIZE=4
  fi
  if [ -z "${THREADS_ARRAY[$STEP]}" ]; then
      THREADS=$((SAMPLE_SIZE*4)) 
  else
      THREADS=${THREADS_ARRAY[$STEP]}
  fi
  if [ -z "${MEM_ARRAY[$STEP]}" ]; then
      MEM=$((SAMPLE_SIZE*10))g 
  else 
      MEM=${MEM_ARRAY[$STEP]}
  fi
  if [ -z "${WALLTIME_ARRAY[$STEP]}" ]; then
      WALLTIME=$((SAMPLE_SIZE*400)) 
  else
      WALLTIME=${WALLTIME_ARRAY[$STEP]}
  fi
  export THREADS=$THREADS
  export WALLTIME=$WALLTIME
  export MEM=$MEM

  echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  echo "Done using pipeline" $VERSION >> $EXPECTED_DONE_FILES
  echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  echo "STEP integrate submitted "  >> $EXPECTED_DONE_FILES
  step_integrate="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem=${MEM} \
    --time=${WALLTIME} \
    --job-name $STEP \
    --export OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},R_VERSION=${R_VERSION},SEULIST=${SEULIST},QUEUE=${QUEUE} \
    --output $JOB_OUTPUT_DIR/logs/%x.o%j \
    $PIPELINE_HOME/integrate/pipeline_integrate.sh"
  step_integrate=$($step_integrate | grep -oP "\d+")
  echo "[Q] STEP integrate         : $step_integrate " >> $EXPECTED_DONE_FILES
  echo "THREADS: ${THREADS} " >> $EXPECTED_DONE_FILES
  echo "WALLTIME: ${WALLTIME} minutes" >> $EXPECTED_DONE_FILES
  echo "MEM: ${MEM} " >> $EXPECTED_DONE_FILES  
    if [ -z "${THREADS_ARRAY[$STEP]}" ]; then
        echo -e "NOTE: THREADS of $STEP is empty, so pipeline assigns it based on #samples,\n you can add THREADS_ARRAY in scrnabox_config.ini"   >> $EXPECTED_DONE_FILES
  fi
  if [ -z "${MEM_ARRAY[$STEP]}" ]; then
      echo -e "NOTE: MEM of $STEP is empty, so pipeline assigns it based on #samples,\n you can add MEM_ARRAY in scrnabox_config.ini"   >> $EXPECTED_DONE_FILES
  fi
  if [ -z "${WALLTIME_ARRAY[$STEP]}" ]; then
      echo  -e "NOTE: WALLTIME of $STEP is empty, so pipeline assigns it based on  #samples,\n you can add WALLTIME_ARRAY in scrnabox_config.ini"   >> $EXPECTED_DONE_FILES
  fi   
  DEPEND_step_integrate="--dependency=afterok:$step_integrate"
  echo_general="STEP integrate: Job Number:$step_integrate"; echo -e $echo_general
# fi 

exit 0