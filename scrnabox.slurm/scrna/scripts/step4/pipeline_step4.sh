#!/bin/bash

umask 002

source $PIPELINE_HOME/tools/utils.sh

if [[ $QUEUE =~ bash ]]; then
    call_parameter $1
fi

# echo timestamp $(date +%s)

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
echo "* STEP4_ANT_LAB:        $STEP4_ANT_LAB"
echo "* SCRNA_METHOD:         ${SCRNA_METHOD}"
echo "-------------------------------------------"
echo "------Parameters used in this step---------"
cat  $OUTPUT_DIR/job_info/parameters/step4_par.txt
echo "-------------------------------------------"
echo -e "------Output of Run------------------------\n\n"
#----------------------------------------------------------------#
# START PIPELINE                                      #
#----------------------------------------------------------------#

if [[ $QUEUE =~ sbatch ]]; then
module load r/$R_VERSION 
fi
# module load r/$R_VERSION
Rscript ${PIPELINE_HOME}/scrna/scripts/step4/scrna_step4.R $OUTPUT_DIR  $R_LIB_PATH  $STEP4_ANT_LAB




