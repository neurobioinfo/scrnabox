#!/bin/bash

# The pipeline is done as part of Project Dark Genome 
# Copyright belong MNI BIOINFO CORE (https://github.com/neurobioinfo)
# The pipeline is written by Saeid Amiri (saeid.amiri@mcgill.ca)


source $OUTPUT_DIR/job_info/configs/scrnabox_config.ini
source $OUTPUT_DIR/job_info/.tmp/temp_config.ini
source $PIPELINE_HOME/tools/utils.sh

# if [ $MODE == 'ALL' ]; then  MODE00=`echo {2..8}`; MODE0=`eval echo $MODE00`; fi 
# if [[ "$MODE" == *"-"* ]]; then
#   MODE00=`echo $MODE | sed  "s/-/../g"  | awk '{print "{"$0"}"}'`
#   MODE0=`eval echo $MODE00`
# else 
#   MODE0=$MODE
# fi

MODE0=$MODE

# ===============================================
# ===============================================

# if  [ ! -f $JOB_OUTPUT_DIR/.tmp/sample.list ]; then
#     search_dir=${OUTPUT_DIR}/step1
#     for entry in "$search_dir"/*
#     do
#       echo $(basename  "$entry") >> ${OUTPUT_DIR}/job_info/.tmp/sample.list
#     done
# fi


# if  [ ! -f $JOB_OUTPUT_DIR/.tmp/sample_dir.list ]; then
#     search_dir=${OUTPUT_DIR}/step1
#     for entry in "$search_dir"/*
#     do
#       echo "$entry" >> ${OUTPUT_DIR}/job_info/.tmp/sample_dir.list
#     done
# fi 

# SAMPLE_SIZE=`wc -l < ${OUTPUT_DIR}/job_info/.tmp/sample.list`


export QUEUE=${QUEUE}

# ===============================================
# STEP 1: 
# ===============================================
#
STEP=step_1

if [[ $QUEUE =~ bash ]] && [[  ${MODE0[@]}  =~  1  ]]; then
  # export ACCOUNT=$ACCOUNT
  if [ ! -d $OUTPUT_DIR/step1 ]; then 
    mkdir -p $OUTPUT_DIR/step1 
  fi
  step1_par_auto=`grep 'par_automated_library_prep=' ${OUTPUT_DIR}/job_info/parameters/step1_par.txt`
  step1_par_auto1=`echo ${step1_par_auto//[[:blank:]]/} | tr 'A-Z' 'a-z'`
  if [[ "${step1_par_auto1}" == "par_automated_library_prep=\"yes\"" ]]; then
        # module load r/$R_VERSION 
        Rscript ${PIPELINE_HOME}/hto/scripts/step1/hto_step1_auto.R $OUTPUT_DIR  $R_LIB_PATH 
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
  bash ${PIPELINE_HOME}/hto/scripts/step1/create_cellranger_hto_local.sh $OUTPUT_DIR /job_info/.tmp/hto_cellranger_local.sh $OUTPUT_DIR/job_info/.tmp/step1_par.txt

  while read item
  do
    # cp ${PIPELINE_HOME}/scrna/scripts/step1/slurm.template $item
    # sed -i $item/slurm.template -e 's/account0/'$ACCOUNT'/'
    cp $OUTPUT_DIR/job_info/.tmp/hto_cellranger_local.sh $item
  done < ${OUTPUT_DIR}/job_info/.tmp/sample_dir.list

  export OUTPUT_DIR=$OUTPUT_DIR
  while read item
  do
      cd ${item}; bash  ${item}/hto_cellranger_local.sh -r ouput_folder &
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
if [[ $QUEUE =~ bash ]] && [[  ${MODE0[@]}  =~  2 ]]  ; then
  echo -e "\n\n-----------------------------------------------------------" >> $EXPECTED_DONE_FILES
  echo -e  "--------Job submitted using pipeline version $VERSION--------"  >> $EXPECTED_DONE_FILES
  echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
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
# fi 
# if [[ $QUEUE =~ bash ]] && [[  ${MODE0[@]}  =~  2  ]]; then
  # echo "just step 2 at" $TIMESTAMP >> $EXPECTED_DONE_FILES
  # echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  # echo "Done using pipeline" $VERSION >> $EXPECTED_DONE_FILES
  # echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  echo "STEP 2 submitted "  >> $EXPECTED_DONE_FILES
  $QUEUE $PIPELINE_HOME/hto/scripts/step2/pipeline_step2.sh --export=OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},R_VERSION=${R_VERSION},SCRNA_METHOD=${SCRNA_METHOD}  &> $JOB_OUTPUT_DIR/logs/${STEP}.$(date +%FT%H.%M.%S) &
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES
  echo "The output is under ${OUTPUT_DIR}/step2" >> $EXPECTED_DONE_FILES  
fi 



# ===============================================
# STEP 3: 
# ===============================================
#
STEP=step_3
if [[ $QUEUE =~ bash ]] && [[  ${MODE0[@]}  =~  3 ]]  ; then
  echo -e "\n\n-----------------------------------------------------------" >> $EXPECTED_DONE_FILES
  echo -e  "--------Job submitted using pipeline version $VERSION--------"  >> $EXPECTED_DONE_FILES
  echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
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
  # echo $TIMESTAMP  >> $EXPECTED_DONE_FILES
  # echo "Done using pipeline" $VERSION >> $EXPECTED_DONE_FILES
  # echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
  echo "STEP 3 submitted"  >> $EXPECTED_DONE_FILES
  $QUEUE $PIPELINE_HOME/hto/scripts/step3/pipeline_step3.sh --export=OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},R_VERSION=${R_VERSION} &> $JOB_OUTPUT_DIR/logs/${STEP}.$(date +%FT%H.%M.%S) 

  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES
  echo "The Output is under ${OUTPUT_DIR}/step3" >> $EXPECTED_DONE_FILES
fi

# ===============================================
# STEP 4: 
# ===============================================
#


STEP=step_4

if [[ $QUEUE =~ bash ]] &&  [[  ${MODE0[@]}  =~  4 ]]  ; then
  echo -e "\n\n-----------------------------------------------------------" >> $EXPECTED_DONE_FILES
  echo -e  "--------Job submitted using pipeline version $VERSION--------"  >> $EXPECTED_DONE_FILES
  echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
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
  if [  -d $OUTPUT_DIR/step4/objs4 ]; then 
    rm -rf  $OUTPUT_DIR/step4/objs4 ; mkdir -p $OUTPUT_DIR/step4/objs4 &
  else
    mkdir -p $OUTPUT_DIR/step4/objs4 &
  fi 
  echo "STEP 4 submitted "  >> $EXPECTED_DONE_FILES
  $QUEUE $PIPELINE_HOME/hto/scripts/step4/pipeline_step4.sh --export=OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},R_VERSION=${R_VERSION},SCRNA_METHOD=${SCRNA_METHOD} &> $JOB_OUTPUT_DIR/logs/${STEP}.$(date +%FT%H.%M.%S) 
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES
  echo "The Output is under ${OUTPUT_DIR}/step4" >> $EXPECTED_DONE_FILES
fi 


# ===============================================
# STEP 5: 
# ===============================================
#
STEP=step_5
if [[ $QUEUE =~ bash ]] && [[  ${MODE0[@]}  =~  5 ]]  ; then
  echo -e "\n\n-----------------------------------------------------------" >> $EXPECTED_DONE_FILES
  echo -e  "--------Job submitted using pipeline version $VERSION--------"  >> $EXPECTED_DONE_FILES
  echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
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
  echo "STEP 4 submitted "  >> $EXPECTED_DONE_FILES
  $QUEUE $PIPELINE_HOME/hto/scripts/step5/pipeline_step5.sh  --export=OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},R_VERSION=${R_VERSION},SCRNA_METHOD=${SCRNA_METHOD} &> $JOB_OUTPUT_DIR/logs/${STEP}.$(date +%FT%H.%M.%S) 
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES  
  echo "The Output is under ${OUTPUT_DIR}/step5" >> $EXPECTED_DONE_FILES
fi 


# ===============================================
# STEP 6: 
# ===============================================
#
STEP=step_6
if [[ $QUEUE =~ bash ]] && [[  ${MODE0[@]}  =~  6 ]]  ; then
  echo -e "\n\n-----------------------------------------------------------" >> $EXPECTED_DONE_FILES
  echo -e  "--------Job submitted using pipeline version $VERSION--------"  >> $EXPECTED_DONE_FILES
  echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
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
  $QUEUE $PIPELINE_HOME/hto/scripts/step6/pipeline_step6.sh --export=OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},R_VERSION=${R_VERSION},SCRNA_METHOD=${SCRNA_METHOD} &> $JOB_OUTPUT_DIR/logs/${STEP}.$(date +%FT%H.%M.%S) 
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES    
    echo "The Output is under ${OUTPUT_DIR}/step6" >> $EXPECTED_DONE_FILES
fi 



# ===============================================
# STEP 7: 
# ===============================================
#
STEP=step_7
if [[ $QUEUE =~ bash ]] && [[  ${MODE0[@]}  =~  7 ]]; then
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
if [[ $QUEUE =~ bash ]] && [[  ${MODE0[@]}  =~  7  ]] && [[   ${STEP7markergsea}  =~  T ]]; then
  echo "STEP 7 submitted "  >> $EXPECTED_DONE_FILES
  $QUEUE $PIPELINE_HOME/hto/scripts/step7/pipeline_step7_markergsea.sh --export=OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},R_VERSION=${R_VERSION},SCRNA_METHOD=${SCRNA_METHOD} &> $JOB_OUTPUT_DIR/logs/${STEP}.$(date +%FT%H.%M.%S) 
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES  
  echo "The Output is under ${OUTPUT_DIR}/step7" >> $EXPECTED_DONE_FILES
fi 

STEP=step_7_knownmarkers
if [[ $QUEUE =~ bash ]] && [[  ${MODE0[@]}  =~  7  ]] &&  [[   ${STEP7knownmarkers}  =~  T ]]; then
  echo "STEP 7 submitted "  >> $EXPECTED_DONE_FILES
  $QUEUE $PIPELINE_HOME/hto/scripts/step7/pipeline_step7_knownmarkers.sh --export=OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},R_VERSION=${R_VERSION},SCRNA_METHOD=${SCRNA_METHOD}  &> $JOB_OUTPUT_DIR/logs/${STEP}.$(date +%FT%H.%M.%S) 
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES  
  echo "The Output is under ${OUTPUT_DIR}/step7" >> $EXPECTED_DONE_FILES
fi 


STEP=step_7_referenceannotation
if [[ $QUEUE =~ bash ]] && [[  ${MODE0[@]}  =~  7  ]] &&  [[   ${STEP7referenceannotation}  =~  T ]]; then
  echo "STEP 7 submitted "  >> $EXPECTED_DONE_FILES
    $QUEUE $PIPELINE_HOME/hto/scripts/step7/pipeline_step7_referenceannotation.sh --export=OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},R_VERSION=${R_VERSION},SCRNA_METHOD=${SCRNA_METHOD}  &> $JOB_OUTPUT_DIR/logs/${STEP}.$(date +%FT%H.%M.%S) 
  echo "-------------------------------------------" >> $EXPECTED_DONE_FILES
  echo "The Output is under ${OUTPUT_DIR}/step7" >> $EXPECTED_DONE_FILES
fi 

STEP=step_7_annotate
if [[ $QUEUE =~ bash ]] && [[  ${MODE0[@]}  =~  7  ]] &&  [[   ${STEP7annotate}  =~  T ]]; then
  echo "STEP 7 submitted "  >> $EXPECTED_DONE_FILES
  $QUEUE $PIPELINE_HOME/hto/scripts/step7/pipeline_step7_annotate.sh --export=OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},R_VERSION=${R_VERSION},SCRNA_METHOD=${SCRNA_METHOD}  &> $JOB_OUTPUT_DIR/logs/${STEP}.$(date +%FT%H.%M.%S) 
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

if [[ $QUEUE =~ bash ]] && [[  ${MODE0[@]}  =~  8 ]]; then
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

if [[ $QUEUE =~ bash ]] && [[  ${MODE0[@]}  =~  8  ]] &&  [[   ${STEP8addmeta}  =~  T ]]; then
  echo "STEP 8 submitted (not dependent on other job) "  >> $EXPECTED_DONE_FILES
    $QUEUE $PIPELINE_HOME/hto/scripts/step8/pipeline_step8_addmeta.sh  --export=OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},R_VERSION=${R_VERSION},SCRNA_METHOD=${SCRNA_METHOD}  &> $JOB_OUTPUT_DIR/logs/${STEP}.$(date +%FT%H.%M.%S) 
    echo "-------------------------------------------" >> $EXPECTED_DONE_FILES    
    echo "The Output is under ${OUTPUT_DIR}/step8" >> $EXPECTED_DONE_FILES
fi 

STEP=step_8_rundge

if [[ $QUEUE =~ bash ]] && [[  ${MODE0[@]}  =~  8  ]] &&  [[  ${STEP8rundge}  =~  T ]]; then
    echo "STEP 8 submitted (not dependent on other job) "  >> $EXPECTED_DONE_FILES
    $QUEUE $PIPELINE_HOME/hto/scripts/step8/pipeline_step8_rundge.sh --export=OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},R_VERSION=${R_VERSION},SCRNA_METHOD=${SCRNA_METHOD} &> $JOB_OUTPUT_DIR/logs/${STEP}.$(date +%FT%H.%M.%S) 
    echo "-------------------------------------------" >> $EXPECTED_DONE_FILES    
    echo "The Output is under ${OUTPUT_DIR}/step8" >> $EXPECTED_DONE_FILES
fi 


exit 0

