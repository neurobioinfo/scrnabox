#!/bin/bash

umask 002
# echo timestamp $(date +%s)

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

module load r/$R_VERSION 
Rscript ${PIPELINE_HOME}/soft/Rcheck/rcheck.R $OUTPUT_DIR $R_LIB_PATH 

