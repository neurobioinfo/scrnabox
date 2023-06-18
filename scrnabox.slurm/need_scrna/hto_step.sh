
STEP=step_2

echo -e "\n\n-----------------------------------------------------------" >> $EXPECTED_DONE_FILES
echo -e  "--------Job submitted using pipeline version $VERSION--------"  >> $EXPECTED_DONE_FILES
echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
echo "STEP 2 submitted"  >> $EXPECTED_DONE_FILES
QUEUE aa0.sh \
--export OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},R_VERSION=${R_VERSION},SCRNA_METHOD=${SCRNA_METHOD} \
--output $JOB_OUTPUT_DIR/logs/%x.o%j \
---bash0 $PIPELINE_HOME/hto/scripts/step2/pipeline.step2.qsub"
echo "-------------------------------------------" >> $EXPECTED_DONE_FILES
echo "The Output is under ${OUTPUT_DIR}/step2" >> $EXPECTED_DONE_FILES


