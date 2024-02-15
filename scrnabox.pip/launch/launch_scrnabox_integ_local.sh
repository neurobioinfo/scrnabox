#!/bin/bash


# The pipeline is done as part of Project Dark Genome 
# Copyright belong MNI BIOINFO CORE (https://github.com/neurobioinfo)
# The pipeline is written by Saeid Amiri (saeid.amiri@mcgill.ca)

source $OUTPUT_DIR/job_info/configs/scrnabox_config.ini
source $OUTPUT_DIR/job_info/.tmp/temp_config.ini
source $PIPELINE_HOME/tools/utils.sh
# ===============================================
# Integrate: Add the Seurat object
# ===============================================
#

export QUEUE=${QUEUE}
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

# if [[ $QUEUE =~ bash]] && [[  ${MODE0[@]}  =~  integrate  ]] ; then
  if [ -z "$SAMPLE_SIZE" ]; then
    SAMPLE_SIZE=4
  fi

  echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  echo "Done using pipeline" $VERSION >> $EXPECTED_DONE_FILES
  echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  echo "STEP integrate submitted "  >> $EXPECTED_DONE_FILES
  $QUEUE $PIPELINE_HOME/integrate/pipeline_integrate.sh --export=OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},R_VERSION=${R_VERSION},SEULIST=${SEULIST} &> $JOB_OUTPUT_DIR/logs/${STEP}.$(date +%FT%H.%M.%S) 
   echo "-------------------------------------------" >> $EXPECTED_DONE_FILES  

# fi 

exit 0

