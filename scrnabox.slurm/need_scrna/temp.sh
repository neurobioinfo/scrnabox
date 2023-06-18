function bash_step () {
    echo "$1" # arguments are accessible through $1, $2,...
}

$1=2  $2=hto $2=$EXPECTED_DONE_FILES   $3=$VERSION 4=$OUTPUT_DIR
bash_step $1 $2 
echo -e "\n\n-----------------------------------------------------------" >> $EXPECTED_DONE_FILES
echo -e  "--------Job submitted using pipeline version $VERSION--------"  >> $EXPECTED_DONE_FILES
echo "-----------------------------------------------------------"  >> $EXPECTED_DONE_FILES
echo "STEP 2"  >> $EXPECTED_DONE_FILES
if [  -d $OUTPUT_DIR/step$1/objs$1 ]; then 
rm -rf  $OUTPUT_DIR/step$1/objs$1 ; mkdir -p $OUTPUT_DIR/step$1/objs$1 &
else
mkdir -p $OUTPUT_DIR/step$1/objs$1 &
fi 
if [  -d $OUTPUT_DIR/step$1/figs$1 ]; then  
rm -rf  $OUTPUT_DIR/step$1/figs$1 ; mkdir -p $OUTPUT_DIR/step$1/figs$1  &
else 
mkdir -p $OUTPUT_DIR/step$1/figs$1 
fi
if [  -d $OUTPUT_DIR/step$1/info$1 ]; then    
rm -rf  $OUTPUT_DIR/step$1/info$1 ; mkdir -p $OUTPUT_DIR/step$1/info$1 &  
else 
mkdir -p $OUTPUT_DIR/step$1/info$1   
fi

--export OUTPUT_DIR=${OUTPUT_DIR},PIPELINE_HOME=${PIPELINE_HOME},R_LIB_PATH=${R_LIB_PATH},R_VERSION=${R_VERSION},SCRNA_METHOD=${SCRNA_METHOD} \
$PIPELINE_HOME/$2/scripts/step$1/pipeline_step$1.qsub"
--output $JOB_OUTPUT_DIR/logs/%x.o%j \
echo "-------------------------------------------" >> $EXPECTED_DONE_FILES
echo "The Output is under ${OUTPUT_DIR}/step$1" >> $EXPECTED_DONE_FILES

