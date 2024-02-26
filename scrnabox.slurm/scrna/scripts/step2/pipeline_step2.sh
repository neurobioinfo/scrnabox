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
echo "* step2 submitted at `date +%FT%H.%M.%S`"
echo "-------------------------------------------"
echo "* PIPELINE_HOME:        $PIPELINE_HOME"
echo "* OUTPUT DIR:                  $OUTPUT_DIR"
echo "* R LIB PATH:                  $R_LIB_PATH"
echo "* scRNA method:          $SCRNA_METHOD"
echo "-------------------------------------------"
echo "------Parameters used in this step---------"
cat  $OUTPUT_DIR/job_info/parameters/step2_par.txt
echo "-------------------------------------------"
echo -e "------Output of Run------------------------\n\n"

#----------------------------------------------------------------#
# START PIPELINE                                      #
#----------------------------------------------------------------#

if [[ $QUEUE =~ sbatch ]]; then
module load r/$R_VERSION 
fi

Rscript ${PIPELINE_HOME}/scrna/scripts/step2/scrna_step2.R  $OUTPUT_DIR  $R_LIB_PATH 

