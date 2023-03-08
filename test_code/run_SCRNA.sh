screen -S run_scrnabox_SCRNA
# mkdir -p /lustre03/project/6070393/COMMON/Dark_Genome/analysis_DarkGenome2weeks_test
export SCRNABOX_HOME=/lustre03/project/6070393/COMMON/Dark_Genome/samamiri/pipeline/scrnaboxbeta2.slurm
export SCRNABOX_PWD=/lustre03/project/6070393/COMMON/Dark_Genome/analysis_DarkGenome2weeks_SCRNA

sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 0

sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 1

sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 2

sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 3 \
--nFRNA 300 \
--nCRNA 21000 \
--pmt 15

sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 4 

sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 5

sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 6

sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--marker T


sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--fta T

sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--enrich T
