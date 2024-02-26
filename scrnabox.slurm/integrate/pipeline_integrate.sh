#!/bin/bash

umask 002
# echo timestamp $(date +%s)
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
echo "* Integrate step submitted at `date +%FT%H.%M.%S`"
echo "-------------------------------------------"
echo "* PIPELINE_HOME:        $PIPELINE_HOME"
echo "* DIR:                  $OUTPUT_DIR"
echo "* SEULIST:                  ${SEULIST}"
echo "-------------------------------------------"
echo "-------------------------------------------"
echo "------Parameters used to run this step-----"
cat  $OUTPUT_DIR/job_info/parameters/stepint_par.txt
echo "-------------------------------------------"
#----------------------------------------------------------------#
# START PIPELINE                                      #
#----------------------------------------------------------------#

if [[ $QUEUE =~ sbatch ]]; then
module load r/$R_VERSION 
fi

Rscript ${PIPELINE_HOME}/integrate/integrate.R $OUTPUT_DIR  $R_LIB_PATH  ${SEULIST}

