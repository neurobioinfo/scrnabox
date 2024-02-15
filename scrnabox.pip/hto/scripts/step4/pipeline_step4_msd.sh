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
echo "* step4 submitted at `date +%FT%H.%M.%S`"
echo "-------------------------------------------"
echo "* PIPELINE_HOME:        $PIPELINE_HOME"
echo "* DIR:                  $OUTPUT_DIR"
echo "* R_LIB_PATH:           $R_LIB_PATH"
# echo "* STEP4_ANT_LAB:            $STEP4_ANT_LAB"
echo "-------------------------------------------"
echo "------Parameters used in this step---------"
cat  $OUTPUT_DIR/job_info/parameters/step4_par.txt
# echo "Anti Label used:"
# cat  $OUTPUT_DIR/job_info/parameters/step4_antibody_label.txt
echo "-------------------------------------------"
echo -e "------Output of Run------------------------\n\n"

#----------------------------------------------------------------#
# START PIPELINE                                      #
#----------------------------------------------------------------#

if [[ $QUEUE =~ sbatch ]]; then
module load r/$R_VERSION 
fi

Rscript ${PIPELINE_HOME}/hto/scripts/step4/hto_step4_msd.R $OUTPUT_DIR  $R_LIB_PATH

