screen -S run_scrnabox_HTO
# mkdir -p /lustre03/project/6070393/COMMON/Dark_Genome/analysis_DarkGenome2weeks_HTO
export SCRNABOX_HOME=/lustre03/project/6070393/COMMON/Dark_Genome/samamiri/pipeline/scrnabox010.slurm
export SCRNABOX_PWD=/lustre03/project/6070393/COMMON/Dark_Genome/analysis_DarkGenome2weeks_HTO

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
--nFRNAL 300 \
--nCRNAU 6500 \
--pmtU 25


sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 4 \
--msd T 

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
--steps 8 \
--dgelist T


sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--genotype T


sh $SCRNABOX_HOME/launch_pipeline.scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--celltype T

