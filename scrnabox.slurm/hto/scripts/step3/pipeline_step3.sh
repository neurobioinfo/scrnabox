#!/bin/bash

umask 002

source $PIPELINE_HOME/tools/utils.sh

if [[ $QUEUE =~ bash ]]; then
   call_parameter $1
fi

#----------------------------------------------------------------#
#                                                     #
# INITIALIZE VARIABLES                                #
#                                                     #
#----------------------------------------------------------------#

echo "-------------------------------------------"
echo "* step3 submitted at `date +%FT%H.%M.%S`"
echo "-------------------------------------------"
echo "* Slurm JOB ID:            $SLURM_JOB_ID"
echo "* PIPELINE_HOME:           $PIPELINE_HOME"
echo "* DIR:                     $OUTPUT_DIR"
echo "* DIR:                     $R_LIB_PATH"
echo "-------------------------------------------"
echo "------Parameters used in this step---------"
cat  $OUTPUT_DIR/job_info/parameters/step3_par.txt
echo "-------------------------------------------"
echo -e "------Output of Run------------------------\n\n"
#----------------------------------------------------------------#
# START PIPELINE                                      #
#----------------------------------------------------------------#


if [[ $QUEUE =~ sbatch ]]; then
   module load r/$R_VERSION 
fi

Rscript ${PIPELINE_HOME}/hto/scripts/step3/hto_step3.R $OUTPUT_DIR  $R_LIB_PATH 

