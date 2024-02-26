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

# if [ $MODE == 'ALL' ]; then  MODE00=`echo {2..8}`; MODE0=`eval echo $MODE00`; fi 
if [[ "$MODE" == *"-"* ]]; then
  MODE00=`echo $MODE | sed  "s/-/../g"  | awk '{print "{"$0"}"}'`
  MODE0=`eval echo $MODE00`
else 
  MODE0=$MODE
fi


export QUEUE=$QUEUE  
# ===============================================
# STEP 1: 
# ===============================================
#
STEP=step_1

if [[ $QUEUE =~ sbatch ]] && [[  ${MODE0[@]}  =~  1  ]]; then
  export ACCOUNT=$ACCOUNT
  if [ ! -d $OUTPUT_DIR/step1 ]; then 
    mkdir -p $OUTPUT_DIR/step1 
  fi
  step1_par_auto=`grep 'par_automated_library_prep=' ${OUTPUT_DIR}/job_info/parameters/step1_par.txt`
  step1_par_auto1=`echo ${step1_par_auto//[[:blank:]]/} | tr 'A-Z' 'a-z'`
  if [[ "${step1_par_auto1}" == "par_automated_library_prep=\"yes\"" ]]; then
        module load r/$R_VERSION 
        Rscript ${PIPELINE_HOME}/scrna/scripts/step1/scrna_step1_auto.R $OUTPUT_DIR  $R_LIB_PATH 
        echo "Note: You are generating samples_info automatically"
  fi

  cp -r  $OUTPUT_DIR/samples_info/* $OUTPUT_DIR/step1 
  if  [ -f $JOB_OUTPUT_DIR/.tmp/sample.list ]; then rm $JOB_OUTPUT_DIR/.tmp/sample.list; touch $JOB_OUTPUT_DIR/.tmp/sample.list ; else touch $JOB_OUTPUT_DIR/.tmp/sample.list; fi 
  if  [ -f $JOB_OUTPUT_DIR/.tmp/sample_dir.list ]; then rm $JOB_OUTPUT_DIR/.tmp/sample_dir.list; touch $JOB_OUTPUT_DIR/.tmp/sample_dir.list ; else touch $JOB_OUTPUT_DIR/.tmp/sample_dir.list; fi  
  if [ -f $OUTPUT_DIR/job_info/.tmp/step1_par.txt ]; then
    rm $OUTPUT_DIR/job_info/.tmp/step1_par.txt
  fi
  
  grep "par_ref_dir_grch=" $OUTPUT_DIR/job_info/parameters/step1_par.txt | sed 's/\"//g' | sed "s/'//g" | sed "s/[[:blank:]]//g" > $OUTPUT_DIR/job_info/.tmp/step1_par.txt
  arr=(par_include_introns par_mempercode par_r1_length par_r2_length par_no_target_umi_filter par_expect_cells par_force_cells par_no_bam par_no_libraries) 
  for item in ${arr[@]}; do
      grep ${item}=  $OUTPUT_DIR/job_info/parameters/step1_par.txt | sed 's/\"//g' | sed "s/'//g" | sed "s/[[:blank:]]//g" | sed 's/[A-Z]/\L&/g' >> $OUTPUT_DIR/job_info/.tmp/step1_par.txt
  done

  search_dir=${OUTPUT_DIR}/step1
  for entry in "$search_dir"/*
  do
    echo "$entry" >> ${OUTPUT_DIR}/job_info/.tmp/sample_dir.list
    echo $(basename  "$entry") >> ${OUTPUT_DIR}/job_info/.tmp/sample.list
  done
  
  source $OUTPUT_DIR/job_info/.tmp/step1_par.txt
  bash ${PIPELINE_HOME}/scrna/scripts/step1/create_cellranger_scrna.sh $OUTPUT_DIR/job_info/.tmp/ scrna_cellranger.slurm.sh $OUTPUT_DIR/job_info/.tmp/step1_par.txt

  while read item
  do
    cp ${PIPELINE_HOME}/scrna/scripts/step1/slurm.template $item
    sed -i $item/slurm.template -e 's/account0/'$ACCOUNT'/'
    cp $OUTPUT_DIR/job_info/.tmp/scrna_cellranger.slurm.sh $item
  done < ${OUTPUT_DIR}/job_info/.tmp/sample_dir.list

  while read item
  do
      cd ${item}; bash  ${item}/scrna_cellranger.slurm.sh -r ouput_folder &
  done < ${OUTPUT_DIR}/job_info/.tmp/sample_dir.list

  wait
  # echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  echo -e  "-------------------------------------------" $VERSION >> $EXPECTED_DONE_FILES  
  echo -e  "--------Job submitted using pipeline-------" $VERSION >> $EXPECTED_DONE_FILES
  echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  echo "step 1  submited " >> $EXPECTED_DONE_FILES
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES
  echo "The Output is under ${OUTPUT_DIR}/step1/" >> $EXPECTED_DONE_FILES
fi 

# ===============================================
# ===============================================

if  [ ! -f $JOB_OUTPUT_DIR/.tmp/sample.list ]; then
    search_dir=${OUTPUT_DIR}/step1
    for entry in "$search_dir"/*
    do
      echo $(basename  "$entry") >> ${OUTPUT_DIR}/job_info/.tmp/sample.list
    done
fi


if  [ ! -f $JOB_OUTPUT_DIR/.tmp/sample_dir.list ]; then
    search_dir=${OUTPUT_DIR}/step1
    for entry in "$search_dir"/*
    do
      echo "$entry" >> ${OUTPUT_DIR}/job_info/.tmp/sample_dir.list
    done
fi 

SAMPLE_SIZE=`wc -l < ${OUTPUT_DIR}/job_info/.tmp/sample.list`


# ===============================================
# STEP 2:
# ===============================================
#
STEP=step_2

if [[ $QUEUE =~ sbatch ]] && [[  ${MODE0[@]}  =~  2 ]]  ; then
  echo -e "\n\n-----------------------------------------------------------" >> $EXPECTED_DONE_FILES
  echo -e  "--------Job submitted using pipeline version $VERSION--------"  >> $EXPECTED_DONE_FILES
  echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  if [ -z "$SAMPLE_SIZE" ]; then
    SAMPLE_SIZE=4
  fi
  if [ -z "${THREADS_ARRAY[$STEP]}" ]; then
      THREADS=$((SAMPLE_SIZE*4)) 
  else
      THREADS=${THREADS_ARRAY[$STEP]}
  fi
  if [ -z "${MEM_ARRAY[$STEP]}" ]; then
      MEM=$((SAMPLE_SIZE*4))g 
  else 
      MEM=${MEM_ARRAY[$STEP]}
  fi
  if [ -z "${WALLTIME_ARRAY[$STEP]}" ]; then
      WALLTIME=$((SAMPLE_SIZE*80)) 
  else
      WALLTIME=${WALLTIME_ARRAY[$STEP]}
  fi
  export THREADS=$THREADS
  export WALLTIME=$WALLTIME
  export MEM=$MEM
  if [  -d $OUTPUT_DIR/step2/objs2 ]; then 
    rm -rf  $OUTPUT_DIR/step2/objs2 ; mkdir -p $OUTPUT_DIR/step2/objs2 &
  else
    mkdir -p $OUTPUT_DIR/step2/objs2 &
  fi 
  if [  -d $OUTPUT_DIR/step2/figs2 ]; then  
    rm -rf  $OUTPUT_DIR/step2/figs2 ; mkdir -p $OUTPUT_DIR/step2/figs2  &
    else 
    mkdir -p $OUTPUT_DIR/step2/figs2 
  fi
  if [  -d $OUTPUT_DIR/step2/info2 ]; then    
    rm -rf  $OUTPUT_DIR/step2/info2 ; mkdir -p $OUTPUT_DIR/step2/info2 &  
    else 
    mkdir -p $OUTPUT_DIR/step2/info2   
  fi
fi 

if [[ $QUEUE =~ sbatch ]] && [[  ${MODE0[@]}  =~  2 ]]  &&  [[  ${MODE0[@]} =~ 1 ]] ; then
  # echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  # echo "Done using pipeline" $VERSION >> $EXPECTED_DONE_FILES
  # echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  echo "STEP 2 submitted following step 1"  >> $EXPECTED_DONE_FILES
  step_2="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem=${MEM} \
    --time=${WALLTIME} \
    --job-name $STEP \
    $DEPEND_CELLRANGER \
    --export OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},R_VERSION=${R_VERSION},SCRNA_METHOD=${SCRNA_METHOD},QUEUE=${QUEUE} \
    --output $JOB_OUTPUT_DIR/logs/%x.o%j \
    $PIPELINE_HOME/scrna/scripts/step2/pipeline_step2.sh"
  step_2=$($step_2 | grep -oP "\d+")
  echo "STEP 2: $step_2 is submitted"
  echo "THREADS: ${THREADS} " >> $EXPECTED_DONE_FILES
  echo "WALLTIME: ${WALLTIME} " >> $EXPECTED_DONE_FILES
  echo "MEM: ${MEM} " >> $EXPECTED_DONE_FILES  
  echo "[Q] STEP 2         : $step_2 " >> $EXPECTED_DONE_FILES
  DEPEND_step_2="--dependency=afterok:$step_2"; echo $DEPEND_step_2
  echo_general="STEP 2: Job Number:$step_2"; echo -e $echo_general
  # echo -e "\n------Parameters used to run this step-----" >> $EXPECTED_DONE_FILES
  # cat  $OUTPUT_DIR/job_info/parameters/step2_par.txt >> $EXPECTED_DONE_FILES
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES
  echo "The result is under ${OUTPUT_DIR}/step2" >> $EXPECTED_DONE_FILES
elif [[ $QUEUE =~ sbatch ]] && [[  ${MODE0[@]}  =~  2  ]]  &&  [[  ${MODE0[@]} != 1 ]]; then
  # echo "just step 2 at" $TIMESTAMP >> $EXPECTED_DONE_FILES
  # echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  # echo "Done using pipeline" $VERSION >> $EXPECTED_DONE_FILES
  # echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  echo "STEP 2 submitted (not dependent on other job) "  >> $EXPECTED_DONE_FILES
  step_2="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem=${MEM} \
    --time=${WALLTIME} \
    --job-name $STEP \
    --export OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},R_VERSION=${R_VERSION},SCRNA_METHOD=${SCRNA_METHOD},QUEUE=${QUEUE} \
    --output $JOB_OUTPUT_DIR/logs/%x.o%j \
    $PIPELINE_HOME/scrna/scripts/step2/pipeline_step2.sh"
  step_2=$($step_2 | grep -oP "\d+")
  echo "[Q] STEP 2         : $step_2 " >> $EXPECTED_DONE_FILES 
  echo "THREADS: ${THREADS} " >> $EXPECTED_DONE_FILES
  echo "WALLTIME: ${WALLTIME} " >> $EXPECTED_DONE_FILES
  echo "MEM: ${MEM} " >> $EXPECTED_DONE_FILES  
  DEPEND_step_2="--dependency=afterok:$step_2"
  echo_general="STEP 2: Job Number:$step_2"; echo -e $echo_general
  # echo -e "\n------Parameters used to run this step-----" >> $EXPECTED_DONE_FILES
  # cat  $OUTPUT_DIR/job_info/parameters/step2_par.txt >> $EXPECTED_DONE_FILES
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES
  echo "The output is under ${OUTPUT_DIR}/step2" >> $EXPECTED_DONE_FILES  
fi 


# ===============================================
# STEP 3: 
# ===============================================
#
STEP=step_3
if [[ $QUEUE =~ sbatch ]] && [[  ${MODE0[@]}  =~  3 ]]  ; then
  echo -e "\n\n-----------------------------------------------------------" >> $EXPECTED_DONE_FILES
  echo -e  "--------Job submitted using pipeline version $VERSION--------"  >> $EXPECTED_DONE_FILES
  echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  if [ -z "$SAMPLE_SIZE" ]; then
    SAMPLE_SIZE=4
  fi
  if [ -z "${THREADS_ARRAY[$STEP]}" ]; then
      THREADS=$((SAMPLE_SIZE*4)) 
  else
      THREADS=${THREADS_ARRAY[$STEP]}
  fi
  if [ -z "${MEM_ARRAY[$STEP]}" ]; then
      MEM=$((SAMPLE_SIZE*4))g 
  else 
      MEM=${MEM_ARRAY[$STEP]}
  fi
  if [ -z "${WALLTIME_ARRAY[$STEP]}" ]; then
      WALLTIME=$((SAMPLE_SIZE*80)) 
  else
      WALLTIME=${WALLTIME_ARRAY[$STEP]}
  fi
  export THREADS=$THREADS
  export WALLTIME=$WALLTIME
  export MEM=$MEM
  if [  -d $OUTPUT_DIR/step3/objs3 ]; then 
    rm -rf  $OUTPUT_DIR/step3/objs3 ; mkdir -p $OUTPUT_DIR/step3/objs3 &
  else
    mkdir -p $OUTPUT_DIR/step3/objs3 &
  fi 
  if [  -d $OUTPUT_DIR/step3/fig3 ]; then  
    rm -rf  $OUTPUT_DIR/step3/figs3 ; mkdir -p $OUTPUT_DIR/step3/figs3  &
    else 
    mkdir -p $OUTPUT_DIR/step3/figs3 
  fi
  if [  -d $OUTPUT_DIR/step3/info3 ]; then    
    rm -rf  $OUTPUT_DIR/step3/info3 ; mkdir -p $OUTPUT_DIR/step3/info3 &  
    else 
    mkdir -p $OUTPUT_DIR/step3/info3   
  fi
fi 


if [[ $QUEUE =~ sbatch ]] && [[  ${MODE0[@]}  =~  3 ]]  &&  [[  ${MODE0[@]} =~ 2 ]] ; then
  # echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  # echo "Done using pipeline" $VERSION >> $EXPECTED_DONE_FILES
  # echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  echo "STEP 3 submitted following step 2"  >> $EXPECTED_DONE_FILES
  step_3="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem=${MEM} \
    --time=${WALLTIME} \
    --job-name $STEP \
    $DEPEND_step_2 \
    --export OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},R_VERSION=${R_VERSION},QUEUE=${QUEUE} \
    --output $JOB_OUTPUT_DIR/logs/%x.o%j \
    $PIPELINE_HOME/scrna/scripts/step3/pipeline_step3.sh"
  step_3=$($step_3 | grep -oP "\d+")
  echo "[Q] STEP 3         : $step_3 " >> $EXPECTED_DONE_FILES
  echo "THREADS: ${THREADS} " >> $EXPECTED_DONE_FILES
  echo "WALLTIME: ${WALLTIME} " >> $EXPECTED_DONE_FILES
  echo "MEM: ${MEM} " >> $EXPECTED_DONE_FILES
  DEPEND_step_3="--dependency=afterok:$step_3"; echo $DEPEND_step_3
  echo_general="STEP 3: Job Number:$step_3"; echo -e $echo_general
    # echo -e "\n------Parameters used to run this step-----" >> $EXPECTED_DONE_FILES
  # cat  $OUTPUT_DIR/job_info/parameters/step3_par.txt >> $EXPECTED_DONE_FILES
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES
  echo "The Output is under ${OUTPUT_DIR}/step3" >> $EXPECTED_DONE_FILES

elif [[ $QUEUE =~ sbatch ]] && [[  ${MODE0[@]}  =~  3  ]]  &&  [[  ${MODE0[@]} != 2 ]]; then
  # echo "just step 3 at" $TIMESTAMP >> $EXPECTED_DONE_FILES
  # echo "Done using pipeline" $VERSION >> $EXPECTED_DONE_FILES
  # echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  # echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  echo "STEP 3 submitted (not dependent on other job) "  >> $EXPECTED_DONE_FILES
  step_3="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem=${MEM} \
    --time=${WALLTIME} \
    --job-name $STEP \
    --export OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},R_VERSION=${R_VERSION},QUEUE=${QUEUE} \
    --output $JOB_OUTPUT_DIR/logs/%x.o%j \
    $PIPELINE_HOME/scrna/scripts/step3/pipeline_step3.sh"
  step_3=$($step_3 | grep -oP "\d+")
  echo "[Q] STEP 3         : $step_3 " >> $EXPECTED_DONE_FILES 
  echo "THREADS: ${THREADS} " >> $EXPECTED_DONE_FILES
  echo "WALLTIME: ${WALLTIME} " >> $EXPECTED_DONE_FILES
  echo "MEM: ${MEM} " >> $EXPECTED_DONE_FILES
  # DEPEND_step_3="--Nodependency on other JOB,\n --JOB Number:$step_3"; echo -e $DEPEND_step_3
  DEPEND_step_3="--dependency=afterok:$step_3"
  echo_general="STEP 3: Job Number:$step_3"; echo -e $echo_general
  # echo -e "\n------Parameters used to run this step-----" >> $EXPECTED_DONE_FILES
  # cat  $OUTPUT_DIR/job_info/parameters/step3_par.txt >> $EXPECTED_DONE_FILES
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES  
  echo "The Output is under ${OUTPUT_DIR}/step3" >> $EXPECTED_DONE_FILES
fi 

# ===============================================
# STEP 4: 
# ===============================================
#
# SCRNA_METHOD=HTO
# echo ${SCRNA_METHOD}
# echo [[  ${SCRNA_METHOD} =~ SCRNA ]] 
# if [[  ${MODE0[@]}  =~  4 ]]  &&  [[  ${SCRNA_METHOD} =~ SCRNA ]] ; then
#  echo -e 'NOTE: Standard scRNA-seq does not have step 4, please go to the next step'
# fi


STEP=step_4
# export STEP4_ANT_LAB=$OUTPUT_DIR/job_info/parameters/step4_antibody_label.txt
if [[ $QUEUE =~ sbatch ]] && [[  ${MODE0[@]}  =~  4 ]]  ; then
  echo -e "\n\n-----------------------------------------------------------" >> $EXPECTED_DONE_FILES
  echo -e  "--------Job submitted using pipeline version $VERSION--------"  >> $EXPECTED_DONE_FILES
  echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  if [ -z "$SAMPLE_SIZE" ]; then
    SAMPLE_SIZE=4
  fi
  if [ -z "${THREADS_ARRAY[$STEP]}" ]; then
      THREADS=$((SAMPLE_SIZE*4)) 
  else
      THREADS=${THREADS_ARRAY[$STEP]}
  fi
  if [ -z "${MEM_ARRAY[$STEP]}" ]; then
      MEM=$((SAMPLE_SIZE*4))g 
  else 
      MEM=${MEM_ARRAY[$STEP]}
  fi
  if [ -z "${WALLTIME_ARRAY[$STEP]}" ]; then
      WALLTIME=$((SAMPLE_SIZE*80)) 
  else
      WALLTIME=${WALLTIME_ARRAY[$STEP]}
  fi
  export THREADS=$THREADS
  export WALLTIME=$WALLTIME
  export MEM=$MEM

  if [ ! -d $OUTPUT_DIR/step4 ]; then 
  mkdir -p $OUTPUT_DIR/step4    
  fi
  if [ ! -d $OUTPUT_DIR/step4/objs4 ]; then 
  mkdir -p $OUTPUT_DIR/step4/objs4
  fi
  if [ ! -d $OUTPUT_DIR/step4/figs4 ]; then 
  mkdir -p $OUTPUT_DIR/step4/figs4
  fi  
  if [ ! -d $OUTPUT_DIR/step4/info4 ]; then 
  mkdir -p $OUTPUT_DIR/step4/info4
  fi    
fi

if [[ $QUEUE =~ sbatch ]] && [[  ${MODE0[@]}  =~  4 ]]  &&  [[  ${MODE0[@]} =~ 3 ]] && [[  ${MSD}  =~  F  ]] ; then
  if [  -d $OUTPUT_DIR/step4/objs4 ]; then 
    rm -rf  $OUTPUT_DIR/step4/objs4 ; mkdir -p $OUTPUT_DIR/step4/objs4 &
  else
    mkdir -p $OUTPUT_DIR/step4/objs4 &
  fi 
  # echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  # echo "Done using pipeline" $VERSION >> $EXPECTED_DONE_FILES
  # echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  echo "STEP 4 submitted following step 3"  >> $EXPECTED_DONE_FILES
  step_4="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem=${MEM} \
    --time=${WALLTIME} \
    --job-name $STEP \
    $DEPEND_step_3 \
    --export=OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},R_VERSION=${R_VERSION},STEP4_ANT_LAB=${STEP4_ANT_LAB},SCRNA_METHOD=${SCRNA_METHOD},QUEUE=${QUEUE} \
    --output $JOB_OUTPUT_DIR/logs/%x.o%j \
    $PIPELINE_HOME/scrna/scripts/step4/pipeline_step4.sh"
  step_4=$($step_4 | grep -oP "\d+")
  echo "[Q] STEP 4         : $step_4 " >> $EXPECTED_DONE_FILES
  echo "THREADS: ${THREADS} " >> $EXPECTED_DONE_FILES
  echo "WALLTIME: ${WALLTIME} " >> $EXPECTED_DONE_FILES
  echo "MEM: ${MEM} " >> $EXPECTED_DONE_FILES  
  # DEPEND_step_4="--dependency=afterok:$step_4"; echo -e $DEPEND_step_4
  DEPEND_step_4="--dependency=afterok:$step_4"
  echo_general="STEP 4: Job Number:$step_4"; echo -e $echo_general
  # echo -e "\n------Parameters used to run this step-----" >> $EXPECTED_DONE_FILES
  # cat  $OUTPUT_DIR/job_info/parameters/step4_par.txt >> $EXPECTED_DONE_FILES
  # echo -e "\nAnti Label used:" >> $EXPECTED_DONE_FILES
  # # cat  $OUTPUT_DIR/job_info/parameters/step4_antibody_label.txt >> $EXPECTED_DONE_FILES
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES
  echo "The Output is under ${OUTPUT_DIR}/step4" >> $EXPECTED_DONE_FILES
elif [[ $QUEUE =~ sbatch ]] && [[  ${MODE0[@]}  =~  4  ]]  &&  [[  ${MODE0[@]} != 3 ]] && [[  ${MSD}  =~  F  ]] ; then
  if [  -d $OUTPUT_DIR/step4/objs4 ]; then 
    rm -rf  $OUTPUT_DIR/step4/objs4 ; mkdir -p $OUTPUT_DIR/step4/objs4 &
  else
    mkdir -p $OUTPUT_DIR/step4/objs4 &
  fi 
  # echo "just step 4 at" $TIMESTAMP >> $EXPECTED_DONE_FILES
  # echo "Done using pipeline" $VERSION >> $EXPECTED_DONE_FILES
  # echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  # echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  echo "STEP 4 submitted (not dependent on other job) "  >> $EXPECTED_DONE_FILES
  step_4="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem=${MEM} \
    --time=${WALLTIME} \
    --job-name $STEP \
    --export=OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},R_VERSION=${R_VERSION},STEP4_ANT_LAB=${STEP4_ANT_LAB},SCRNA_METHOD=${SCRNA_METHOD},QUEUE=${QUEUE} \
    --output $JOB_OUTPUT_DIR/logs/%x.o%j \
    $PIPELINE_HOME/scrna/scripts/step4/pipeline_step4.sh"
  step_4=$($step_4 | grep -oP "\d+")
  echo "[Q] STEP 4         : $step_4 " >> $EXPECTED_DONE_FILES 
  echo "THREADS: ${THREADS} " >> $EXPECTED_DONE_FILES
  echo "WALLTIME: ${WALLTIME} " >> $EXPECTED_DONE_FILES
  echo "MEM: ${MEM} " >> $EXPECTED_DONE_FILES  
  # DEPEND_step_4="--Nodependency on other JOB, \n --JOB Number:$step_4"; echo -e $DEPEND_step_4
  DEPEND_step_4="--dependency=afterok:$step_4"
  echo_general="STEP 4: Job Number:$step_4"; echo -e $echo_general
  # echo -e "\n------Parameters used to run this step-----" >> $EXPECTED_DONE_FILES
  # cat  $OUTPUT_DIR/job_info/parameters/step4_par.txt >> $EXPECTED_DONE_FILES
  # echo -e "\nAnti Label used:" >> $EXPECTED_DONE_FILES
  # # cat  $OUTPUT_DIR/job_info/parameters/step4_antibody_label.txt >> $EXPECTED_DONE_FILES
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES
  echo "The Output is under ${OUTPUT_DIR}/step4" >> $EXPECTED_DONE_FILES
fi 

if [[ $QUEUE =~ sbatch ]] && [[  ${MODE0[@]}  =~  4 ]] && [[  ${MSD}  =~  T  ]] &&   [[  ${SCRNA_METHOD} =~ SCRNA ]] ; then
  echo -e 'NOTE: MSD  is not for Standard scRNA-seq'
fi

if [[ $QUEUE =~ sbatch ]] && [[  ${MODE0[@]}  =~  4 ]] && [[  ${MSD}  =~  T  ]] &&  [[  ${SCRNA_METHOD} =~ HTO ]] ; then
  # echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  # echo "Done using pipeline" $VERSION >> $EXPECTED_DONE_FILES
  # echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  echo "STEP 4-MSD "  >> $EXPECTED_DONE_FILES
  step_4="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem-per-cpu=6g \
    --time=${WALLTIME} \
    --job-name ${STEP}_MSD \
    --export=OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},R_VERSION=${R_VERSION},STEP4_ANT_LAB=${STEP4_ANT_LAB},QUEUE=${QUEUE} \
    --output $JOB_OUTPUT_DIR/logs/%x.o%j \
    $PIPELINE_HOME/scrna/scripts/step4/pipeline_step4_msd.sh"
  step_4=$($step_4 | grep -oP "\d+")
  echo "[Q] STEP 4         : $step_4 " >> $EXPECTED_DONE_FILES 
  echo "THREADS: ${THREADS} " >> $EXPECTED_DONE_FILES
  echo "WALLTIME: ${WALLTIME} " >> $EXPECTED_DONE_FILES
  echo "MEM: ${MEM} " >> $EXPECTED_DONE_FILES  
  # DEPEND_step_4="--Nodependency on other JOB. \n --JOB Number:$step_4"; echo -e $DEPEND_step_4
  DEPEND_step_4="--dependency=afterok:$step_4"
  echo_general="STEP 4: Job Number:$step_4"; echo -e $echo_general
  echo -e "\n------Parameters used to run this step-----" >> $EXPECTED_DONE_FILES
  # cat  $OUTPUT_DIR/job_info/parameters/step4_par.txt >> $EXPECTED_DONE_FILES
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES  
    echo "The Output is under ${OUTPUT_DIR}/step4" >> $EXPECTED_DONE_FILES
fi 

# ===============================================
# STEP 5: 
# ===============================================
#
STEP=step_5
if [[ $QUEUE =~ sbatch ]] && [[  ${MODE0[@]}  =~  5 ]]  ; then
  echo -e "\n\n-----------------------------------------------------------" >> $EXPECTED_DONE_FILES
  echo -e  "--------Job submitted using pipeline version $VERSION--------"  >> $EXPECTED_DONE_FILES
  echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  if [ -z "$SAMPLE_SIZE" ]; then
    SAMPLE_SIZE=4
  fi
  if [ -z "${THREADS_ARRAY[$STEP]}" ]; then
      THREADS=$((SAMPLE_SIZE*4)) 
  else
      THREADS=${THREADS_ARRAY[$STEP]}
  fi
  if [ -z "${MEM_ARRAY[$STEP]}" ]; then
      MEM=$((SAMPLE_SIZE*4))g 
  else 
      MEM=${MEM_ARRAY[$STEP]}
  fi
  if [ -z "${WALLTIME_ARRAY[$STEP]}" ]; then
      WALLTIME=$((SAMPLE_SIZE*80)) 
  else
      WALLTIME=${WALLTIME_ARRAY[$STEP]}
  fi
  export THREADS=$THREADS
  export WALLTIME=$WALLTIME
  export MEM=$MEM
  if [  -d $OUTPUT_DIR/step5/objs5 ]; then 
    rm -rf  $OUTPUT_DIR/step5/objs5 ; mkdir -p $OUTPUT_DIR/step5/objs5 &
  else
    mkdir -p $OUTPUT_DIR/step5/objs5 &
  fi 
  if [  -d $OUTPUT_DIR/step5/fig5 ]; then  
    rm -rf  $OUTPUT_DIR/step5/figs5 ; mkdir -p $OUTPUT_DIR/step5/figs5  &
    else 
    mkdir -p $OUTPUT_DIR/step5/figs5 
  fi
  if [  -d $OUTPUT_DIR/step5/info5 ]; then    
    rm -rf  $OUTPUT_DIR/step5/info5 ; mkdir -p $OUTPUT_DIR/step5/info5 &  
    else 
    mkdir -p $OUTPUT_DIR/step5/info5   
  fi
fi 
if [[ $QUEUE =~ sbatch ]] && [[  ${MODE0[@]}  =~  5 ]]  &&  [[  ${MODE0[@]} =~ 4 ]] ; then
  # echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  # echo "Done using pipeline" $VERSION >> $EXPECTED_DONE_FILES
  # echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  echo "STEP 5 submitted following step 4"  >> $EXPECTED_DONE_FILES
  step_5="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem=${MEM} \
    --time=${WALLTIME} \
    --job-name $STEP \
    $DEPEND_step_4 \
    --export=OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},R_VERSION=${R_VERSION},SCRNA_METHOD=${SCRNA_METHOD},QUEUE=${QUEUE} \
    --output $JOB_OUTPUT_DIR/logs/%x.o%j \
    $PIPELINE_HOME/scrna/scripts/step5/pipeline_step5.sh"
  step_5=$($step_5 | grep -oP "\d+")
  echo "[Q] STEP 5         : $step_5 " >> $EXPECTED_DONE_FILES
  echo "THREADS: ${THREADS} " >> $EXPECTED_DONE_FILES
  echo "WALLTIME: ${WALLTIME} " >> $EXPECTED_DONE_FILES
  echo "MEM: ${MEM} " >> $EXPECTED_DONE_FILES  
  # DEPEND_step_5="--dependency=afterok:$step_5"; echo $DEPEND_step_5
  DEPEND_step_5="--dependency=afterok:$step_5"
  echo_general="STEP 5: Job Number:$step_5"; echo -e $echo_general
  # echo -e "\n------Parameters used to run this step-----" >> $EXPECTED_DONE_FILES
  # cat  $OUTPUT_DIR/job_info/parameters/step5_par.txt >> $EXPECTED_DONE_FILES
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES  
   echo "The Output is under ${OUTPUT_DIR}/step5" >> $EXPECTED_DONE_FILES 
elif [[ $QUEUE =~ sbatch ]] && [[  ${MODE0[@]}  =~  5  ]]  &&  [[  ${MODE0[@]} != 4 ]]; then
  # echo "just step 5 at" $TIMESTAMP >> $EXPECTED_DONE_FILES
  # echo "Done using pipeline" $VERSION >> $EXPECTED_DONE_FILES
  # echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  # echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  echo "STEP 5 submitted (not dependent on other job) "  >> $EXPECTED_DONE_FILES
  step_5="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem=${MEM} \
    --time=${WALLTIME} \
    --job-name $STEP \
    --export=OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},R_VERSION=${R_VERSION},SCRNA_METHOD=${SCRNA_METHOD},QUEUE=${QUEUE} \
    --output $JOB_OUTPUT_DIR/logs/%x.o%j \
    $PIPELINE_HOME/scrna/scripts/step5/pipeline_step5.sh"
  step_5=$($step_5 | grep -oP "\d+")
  echo "[Q] STEP 5         : $step_5 " >> $EXPECTED_DONE_FILES 
  echo "THREADS: ${THREADS} " >> $EXPECTED_DONE_FILES
  echo "WALLTIME: ${WALLTIME} " >> $EXPECTED_DONE_FILES
  echo "MEM: ${MEM} " >> $EXPECTED_DONE_FILES  
  # DEPEND_step_5="--Nodependency on other JOB.\n --JOB Number:$step_5"; echo $DEPEND_step_5
  DEPEND_step_5="--dependency=afterok:$step_5"
  echo_general="STEP 5: Job Number:$step_5"; echo -e $echo_general
  # echo -e "\n------Parameters used to run this step-----" >> $EXPECTED_DONE_FILES
  # cat  $OUTPUT_DIR/job_info/parameters/step5_par.txt >> $EXPECTED_DONE_FILES
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES  
  echo "The Output is under ${OUTPUT_DIR}/step5" >> $EXPECTED_DONE_FILES
fi 


# ===============================================
# STEP 6: 
# ===============================================
#
STEP=step_6
if [[ $QUEUE =~ sbatch ]] && [[  ${MODE0[@]}  =~  6 ]]  ; then
  echo -e "\n\n-----------------------------------------------------------" >> $EXPECTED_DONE_FILES
  echo -e  "--------Job submitted using pipeline version $VERSION--------"  >> $EXPECTED_DONE_FILES
  echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  if [ -z "$SAMPLE_SIZE" ]; then
    SAMPLE_SIZE=4
  fi
  if [ -z "${THREADS_ARRAY[$STEP]}" ]; then
      THREADS=$((SAMPLE_SIZE*4)) 
  else
      THREADS=${THREADS_ARRAY[$STEP]}
  fi
  if [ -z "${MEM_ARRAY[$STEP]}" ]; then
      MEM=$((SAMPLE_SIZE*4))g 
  else 
      MEM=${MEM_ARRAY[$STEP]}
  fi
  if [ -z "${WALLTIME_ARRAY[$STEP]}" ]; then
      WALLTIME=$((SAMPLE_SIZE*80)) 
  else
      WALLTIME=${WALLTIME_ARRAY[$STEP]}
  fi
  export THREADS=$THREADS
  export WALLTIME=$WALLTIME
  export MEM=$MEM
  if [  -d $OUTPUT_DIR/step6/objs6 ]; then 
    rm -rf  $OUTPUT_DIR/step6/objs6 ; mkdir -p $OUTPUT_DIR/step6/objs6 &
  else
    mkdir -p $OUTPUT_DIR/step6/objs6 &
  fi 
  if [  -d $OUTPUT_DIR/step6/fig6 ]; then  
    rm -rf  $OUTPUT_DIR/step6/figs6 ; mkdir -p $OUTPUT_DIR/step6/figs6  &
    else 
    mkdir -p $OUTPUT_DIR/step6/figs6 
  fi
  if [  -d $OUTPUT_DIR/step6/info6 ]; then    
    rm -rf  $OUTPUT_DIR/step6/info6 ; mkdir -p $OUTPUT_DIR/step6/info6 &  
    else 
    mkdir -p $OUTPUT_DIR/step6/info6   
  fi
fi 
if [[ $QUEUE =~ sbatch ]] && [[  ${MODE0[@]}  =~  6 ]]  &&  [[  ${MODE0[@]} =~ 5 ]] ; then
  # echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  # echo "Done using pipeline" $VERSION >> $EXPECTED_DONE_FILES
  # echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  echo "STEP 6 submitted following step 5"  >> $EXPECTED_DONE_FILES
  step_6="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem=${MEM} \
    --time=${WALLTIME} \
    --job-name $STEP \
    $DEPEND_step_5 \
    --export=OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},R_VERSION=${R_VERSION},SCRNA_METHOD=${SCRNA_METHOD},QUEUE=${QUEUE} \
    --output $JOB_OUTPUT_DIR/logs/%x.o%j \
    $PIPELINE_HOME/scrna/scripts/step6/pipeline_step6.sh"
  step_6=$($step_6 | grep -oP "\d+")
  echo "[Q] STEP 6         : $step_6 " >> $EXPECTED_DONE_FILES
  echo "THREADS: ${THREADS} " >> $EXPECTED_DONE_FILES
  echo "WALLTIME: ${WALLTIME} " >> $EXPECTED_DONE_FILES
  echo "MEM: ${MEM} " >> $EXPECTED_DONE_FILES  
  DEPEND_step_6="--dependency=afterok:$step_6"
  echo_general="STEP 6: Job Number:$step_6"; echo -e $echo_general
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES
  echo "The Output is under ${OUTPUT_DIR}/step6" >> $EXPECTED_DONE_FILES
elif [[ $QUEUE =~ sbatch ]] && [[  ${MODE0[@]}  =~  6  ]]  &&  [[  ${MODE0[@]} != 5 ]]; then
  # echo "just step 6 at" $TIMESTAMP >> $EXPECTED_DONE_FILES
  # echo "Done using pipeline" $VERSION >> $EXPECTED_DONE_FILES
  # echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  # echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  echo "STEP 6 submitted (not dependent on other job) "  >> $EXPECTED_DONE_FILES
  step_6="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem=${MEM} \
    --time=${WALLTIME} \
    --job-name $STEP \
    --export=OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},R_VERSION=${R_VERSION},SCRNA_METHOD=${SCRNA_METHOD},QUEUE=${QUEUE} \
    --output $JOB_OUTPUT_DIR/logs/%x.o%j \
    $PIPELINE_HOME/scrna/scripts/step6/pipeline_step6.sh"
  step_6=$($step_6 | grep -oP "\d+")
  echo "[Q] STEP 6         : $step_6 " >> $EXPECTED_DONE_FILES 
  echo "THREADS: ${THREADS} " >> $EXPECTED_DONE_FILES
  echo "WALLTIME: ${WALLTIME} " >> $EXPECTED_DONE_FILES
  echo "MEM: ${MEM} " >> $EXPECTED_DONE_FILES  
  DEPEND_step_6="--dependency=afterok:$step_6"
  echo_general="STEP 6: Job Number:$step_6"; echo -e $echo_general
  # echo -e "\n------Parameters used to run this step-----" >> $EXPECTED_DONE_FILES
  # cat  $OUTPUT_DIR/job_info/parameters/step6_par.txt >> $EXPECTED_DONE_FILES
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES    
    echo "The Output is under ${OUTPUT_DIR}/step6" >> $EXPECTED_DONE_FILES
fi 



# ===============================================
# STEP 7: 
# ===============================================
#
STEP=step_7
if [[ $QUEUE =~ sbatch ]] && [[  ${MODE0[@]}  =~  7 ]]; then
  echo -e "\n\n-----------------------------------------------------------" >> $EXPECTED_DONE_FILES
  echo -e  "--------Job submitted using pipeline version $VERSION--------"  >> $EXPECTED_DONE_FILES
  echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  if [ ! -d $OUTPUT_DIR/step7 ]; then 
  mkdir -p $OUTPUT_DIR/step7    
  fi
  if [ ! -d $OUTPUT_DIR/step7/objs7 ]; then 
  mkdir -p $OUTPUT_DIR/step7/objs7
  fi
  if [ ! -d $OUTPUT_DIR/step7/figs7 ]; then 
  mkdir -p $OUTPUT_DIR/step7/figs7
  fi  
  if [ ! -d $OUTPUT_DIR/step7/info7 ]; then 
  mkdir -p $OUTPUT_DIR/step7/info7
  fi    
fi

STEP=step_7_markergsea
if [[ $QUEUE =~ sbatch ]] && [[  ${MODE0[@]}  =~  7  ]] && [[   ${STEP7markergsea}  =~  T ]]; then
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
      WALLTIME=$((SAMPLE_SIZE*80)) 
  else
      WALLTIME=${WALLTIME_ARRAY[$STEP]}
  fi
  export THREADS=$THREADS
  export WALLTIME=$WALLTIME
  export MEM=$MEM
  echo "STEP 7 submitted (not dependent on other job) "  >> $EXPECTED_DONE_FILES
  step_7="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem=${MEM} \
    --time=${WALLTIME} \
    --job-name ${STEP} \
    --export OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},R_VERSION=${R_VERSION},SCRNA_METHOD=${SCRNA_METHOD},QUEUE=${QUEUE} \
    --output $JOB_OUTPUT_DIR/logs/%x.o%j \
    $PIPELINE_HOME/scrna/scripts/step7/pipeline_step7_markergsea.sh"
  step_7=$($step_7 | grep -oP "\d+")
  echo "[Q] STEP 7 markergsea         : $step_7 " >> $EXPECTED_DONE_FILES 
  echo "THREADS: ${THREADS} " >> $EXPECTED_DONE_FILES
  echo "WALLTIME: ${WALLTIME} " >> $EXPECTED_DONE_FILES
  echo "MEM: ${MEM} " >> $EXPECTED_DONE_FILES  
  # DEPEND_step_7="--Nodependency on other JOB. \n --JOB Number:$step_7"; echo -e $DEPEND_step_7
  DEPEND_step_7="--dependency=afterok:$step_7"
  echo_general="STEP 7: Job Number:$step_7"; echo -e $echo_general
  # echo -e "\n------Parameters used to run this step-----" >> $EXPECTED_DONE_FILES
  # cat  $OUTPUT_DIR/job_info/parameters/step7_par.txt >> $EXPECTED_DONE_FILES
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES  
  echo "The Output is under ${OUTPUT_DIR}/step7" >> $EXPECTED_DONE_FILES
fi 

STEP=step_7_knownmarkers
if [[ $QUEUE =~ sbatch ]] && [[  ${MODE0[@]}  =~  7  ]] &&  [[   ${STEP7knownmarkers}  =~  T ]]; then
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
      WALLTIME=$((SAMPLE_SIZE*100)) 
  else
      WALLTIME=${WALLTIME_ARRAY[$STEP]}
  fi
  export THREADS=$THREADS
  export WALLTIME=$WALLTIME
  export MEM=$MEM
  echo "STEP 7 submitted (not dependent on other job) "  >> $EXPECTED_DONE_FILES
  step_7="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem=${MEM} \
    --time=${WALLTIME} \
    --job-name ${STEP} \
    --export OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},R_VERSION=${R_VERSION},SCRNA_METHOD=${SCRNA_METHOD},QUEUE=${QUEUE} \
    --output $JOB_OUTPUT_DIR/logs/%x.o%j \
    $PIPELINE_HOME/scrna/scripts/step7/pipeline_step7_knownmarkers.sh"
  step_7=$($step_7 | grep -oP "\d+")
  echo "[Q] STEP 7 knownmarkers         : $step_7 " >> $EXPECTED_DONE_FILES 
  echo "THREADS: ${THREADS} " >> $EXPECTED_DONE_FILES
  echo "WALLTIME: ${WALLTIME} " >> $EXPECTED_DONE_FILES
  echo "MEM: ${MEM} " >> $EXPECTED_DONE_FILES  
  DEPEND_step_7="--dependency=afterok:$step_7"
  echo_general="STEP 7: Job Number:$step_7"; echo -e $echo_general
  # echo -e "\n------Parameters used to run this step-----" >> $EXPECTED_DONE_FILES
  # cat  $OUTPUT_DIR/job_info/parameters/step7_par.txt >> $EXPECTED_DONE_FILES
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES  
  echo "The Output is under ${OUTPUT_DIR}/step7" >> $EXPECTED_DONE_FILES
fi 


STEP=step_7_referenceannotation
if [[ $QUEUE =~ sbatch ]] && [[  ${MODE0[@]}  =~  7  ]] &&  [[   ${STEP7referenceannotation}  =~  T ]]; then
  if [ -z "$SAMPLE_SIZE" ]; then
    SAMPLE_SIZE=4
  fi
  if [ -z "${THREADS_ARRAY[$STEP]}" ]; then
      THREADS=$((SAMPLE_SIZE*4)) 
  else
      THREADS=${THREADS_ARRAY[$STEP]}
  fi
  if [ -z "${MEM_ARRAY[$STEP]}" ]; then
      MEM=$((SAMPLE_SIZE*30))g 
  else 
      MEM=${MEM_ARRAY[$STEP]}
  fi
  if [ -z "${WALLTIME_ARRAY[$STEP]}" ]; then
      WALLTIME=$((SAMPLE_SIZE*100)) 
  else
      WALLTIME=${WALLTIME_ARRAY[$STEP]}
  fi
  export THREADS=$THREADS
  export WALLTIME=$WALLTIME
  export MEM=$MEM
  echo "STEP 7 submitted (not dependent on other job) "  >> $EXPECTED_DONE_FILES
  step_7="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem=${MEM} \
    --time=${WALLTIME} \
    --job-name ${STEP}_referenceannotation \
    --export OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},R_VERSION=${R_VERSION},SCRNA_METHOD=${SCRNA_METHOD},QUEUE=${QUEUE} \
    --output $JOB_OUTPUT_DIR/logs/%x.o%j \
    $PIPELINE_HOME/scrna/scripts/step7/pipeline_step7_referenceannotation.sh"
  step_7=$($step_7 | grep -oP "\d+")
  echo "[Q] STEP 7 referenceannotation         : $step_7 " >> $EXPECTED_DONE_FILES 
  echo "THREADS: ${THREADS} " >> $EXPECTED_DONE_FILES
  echo "WALLTIME: ${WALLTIME} " >> $EXPECTED_DONE_FILES
  echo "MEM: ${MEM} " >> $EXPECTED_DONE_FILES  
  # DEPEND_step_7="--Nodependency on other JOB. \n --JOB Number:$step_7"; echo -e $DEPEND_step_7
    DEPEND_step_7="--dependency=afterok:$step_7"
  echo_general="STEP 7: Job Number:$step_7"; echo -e $echo_general
    # echo -e "\n------Parameters used to run this step-----" >> $EXPECTED_DONE_FILES
  # cat  $OUTPUT_DIR/job_info/parameters/step7_par.txt >> $EXPECTED_DONE_FILES
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES
  echo "The Output is under ${OUTPUT_DIR}/step7" >> $EXPECTED_DONE_FILES
fi 

STEP=step_7_annotate
if [[ $QUEUE =~ sbatch ]] && [[  ${MODE0[@]}  =~  7  ]] &&  [[   ${STEP7annotate}  =~  T ]]; then
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
      WALLTIME=$((SAMPLE_SIZE*100)) 
  else
      WALLTIME=${WALLTIME_ARRAY[$STEP]}
  fi
  export THREADS=$THREADS
  export WALLTIME=$WALLTIME
  export MEM=$MEM
  echo "STEP 7 submitted (not dependent on other job) "  >> $EXPECTED_DONE_FILES
  step_7="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem=${MEM} \
    --time=${WALLTIME} \
    --job-name ${STEP} \
    --export OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},R_VERSION=${R_VERSION},SCRNA_METHOD=${SCRNA_METHOD},QUEUE=${QUEUE} \
    --output $JOB_OUTPUT_DIR/logs/%x.o%j \
    $PIPELINE_HOME/scrna/scripts/step7/pipeline_step7_annotate.sh"
  step_7=$($step_7 | grep -oP "\d+")
  echo "[Q] STEP 7 annotate         : $step_7 " >> $EXPECTED_DONE_FILES 
  echo "THREADS: ${THREADS} " >> $EXPECTED_DONE_FILES
  echo "WALLTIME: ${WALLTIME} " >> $EXPECTED_DONE_FILES
  echo "MEM: ${MEM} " >> $EXPECTED_DONE_FILES  
  # DEPEND_step_7="--Nodependency on other JOB. \n --JOB Number:$step_7"; echo -e $DEPEND_step_7
    DEPEND_step_7="--dependency=afterok:$step_7"
  echo_general="STEP 7: Job Number:$step_7"; echo -e $echo_general
    # echo -e "\n------Parameters used to run this step-----" >> $EXPECTED_DONE_FILES
  # cat  $OUTPUT_DIR/job_info/parameters/step7_par.txt >> $EXPECTED_DONE_FILES
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES
  echo "The Output is under ${OUTPUT_DIR}/step7" >> $EXPECTED_DONE_FILES
fi 


# ===============================================
# STEP 8: 
# ===============================================
#
STEP=step_8
# THRm=`wc -l < ${OUTPUT_DIR}/job_info/parameters/step8_contrast_genotype.txt`
# THRi=`wc -l < ${OUTPUT_DIR}/job_info/parameters/step8_contrast_celltype.txt`

if [[ $QUEUE =~ sbatch ]] && [[  ${MODE0[@]}  =~  8 ]]; then
  echo -e "\n\n-----------------------------------------------------------" >> $EXPECTED_DONE_FILES
  echo -e  "--------Job submitted using pipeline version $VERSION--------"  >> $EXPECTED_DONE_FILES
  echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  if [ ! -d $OUTPUT_DIR/step8 ]; then 
    mkdir -p $OUTPUT_DIR/step8    
  fi
  # if [ ! -d $OUTPUT_DIR/step8/figs8 ]; then 
  #   mkdir -p $OUTPUT_DIR/step8/figs8
  # fi  
  if [ ! -d $OUTPUT_DIR/step8/objs8 ]; then 
    mkdir -p $OUTPUT_DIR/step8/objs8
  fi
  if [ ! -d $OUTPUT_DIR/step8/info8 ]; then 
    mkdir -p $OUTPUT_DIR/step8/info8
  fi      
fi

STEP=step_8_addmeta

if [[ $QUEUE =~ sbatch ]] && [[  ${MODE0[@]}  =~  8  ]] &&  [[   ${STEP8addmeta}  =~  T ]]; then
  if [ -z "$SAMPLE_SIZE" ]; then
    SAMPLE_SIZE=4
  fi
  if [ -z "${THREADS_ARRAY[$STEP]}" ]; then
      THREADS=$((SAMPLE_SIZE*8)) 
  else
      THREADS=${THREADS_ARRAY[$STEP]}
  fi
  if [ -z "${MEM_ARRAY[$STEP]}" ]; then
      MEM=$((SAMPLE_SIZE*20))g 
  else 
      MEM=${MEM_ARRAY[$STEP]}
  fi
  if [ -z "${WALLTIME_ARRAY[$STEP]}" ]; then
      WALLTIME=$((SAMPLE_SIZE*100)) 
  else
      WALLTIME=${WALLTIME_ARRAY[$STEP]}
  fi
  export THREADS=$THREADS
  export WALLTIME=$WALLTIME
  export MEM=$MEM
  echo "STEP 8 submitted (not dependent on other job) "  >> $EXPECTED_DONE_FILES
  step_8="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem=${MEM} \
    --time=${WALLTIME} \
    --job-name ${STEP} \
    --export OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},R_VERSION=${R_VERSION},SCRNA_METHOD=${SCRNA_METHOD},QUEUE=${QUEUE} \
    --output $JOB_OUTPUT_DIR/logs/%x.o%j \
    $PIPELINE_HOME/scrna/scripts/step8/pipeline_step8_addmeta.sh"
  step_8=$($step_8 | grep -oP "\d+")
  echo "[Q] STEP 8         : $step_8 " >> $EXPECTED_DONE_FILES 
  echo "THREADS: ${THREADS} " >> $EXPECTED_DONE_FILES
  echo "WALLTIME: ${WALLTIME} " >> $EXPECTED_DONE_FILES
  echo "MEM: ${MEM} " >> $EXPECTED_DONE_FILES  
  # DEPEND_step_8_dgelist="--dependency=afterok:$step_8"; echo $DEPEND_step_8_dgelist
    DEPEND_step_8="--dependency=afterok:$step_8"
  echo_general="STEP 8: Job Number:$step_8"; echo -e $echo_general
  # echo -e "\n------Parameters used to run this step-----" >> $EXPECTED_DONE_FILES
  # cat  $OUTPUT_DIR/job_info/parameters/step8_par.txt >> $EXPECTED_DONE_FILES
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES    
    echo "The Output is under ${OUTPUT_DIR}/step8" >> $EXPECTED_DONE_FILES
fi 

STEP=step_8_rundge

if [[ $QUEUE =~ sbatch ]] && [[  ${MODE0[@]}  =~  8  ]] &&  [[  ${STEP8rundge}  =~  T ]]; then
  if [ -z "$SAMPLE_SIZE" ]; then
    SAMPLE_SIZE=4
  fi
  if [ -z "${THREADS_ARRAY[$STEP]}" ]; then
      THREADS=$((SAMPLE_SIZE*8)) 
  else
      THREADS=${THREADS_ARRAY[$STEP]}
  fi
  if [ -z "${MEM_ARRAY[$STEP]}" ]; then
      MEM=$((SAMPLE_SIZE*40))g 
  else 
      MEM=${MEM_ARRAY[$STEP]}
  fi
  if [ -z "${WALLTIME_ARRAY[$STEP]}" ]; then
      WALLTIME=$((SAMPLE_SIZE*100)) 
  else
      WALLTIME=${WALLTIME_ARRAY[$STEP]}
  fi
  export THREADS=$THREADS
  export WALLTIME=$WALLTIME
  export MEM=$MEM
  echo "STEP 8 submitted (not dependent on other job) "  >> $EXPECTED_DONE_FILES
  step_8="$QUEUE -A $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem=${MEM} \
    --time=${WALLTIME} \
    --job-name ${STEP} \
    --export OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},R_VERSION=${R_VERSION},SCRNA_METHOD=${SCRNA_METHOD},QUEUE=${QUEUE} \
    --output $JOB_OUTPUT_DIR/logs/%x.o%j \
    $PIPELINE_HOME/scrna/scripts/step8/pipeline_step8_rundge.sh"
  step_8=$($step_8 | grep -oP "\d+")
  echo "[Q] STEP 8         : $step_8 " >> $EXPECTED_DONE_FILES 
  echo "THREADS: ${THREADS} " >> $EXPECTED_DONE_FILES
  echo "WALLTIME: ${WALLTIME} " >> $EXPECTED_DONE_FILES
  echo "MEM: ${MEM} " >> $EXPECTED_DONE_FILES  
  # DEPEND_step_8_dgelist="--dependency=afterok:$step_8"; echo $DEPEND_step_8_dgelist
    DEPEND_step_8="--dependency=afterok:$step_8"
  echo_general="STEP 8: Job Number:$step_8"; echo -e $echo_general
  # echo -e "\n------Parameters used to run this step-----" >> $EXPECTED_DONE_FILES
  # cat  $OUTPUT_DIR/job_info/parameters/step8_par.txt >> $EXPECTED_DONE_FILES
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES    
    echo "The Output is under ${OUTPUT_DIR}/step8" >> $EXPECTED_DONE_FILES
fi 


exit 0
