screen -S run_scrnabox_HTO
# mkdir -p /lustre03/project/6070393/COMMON/Dark_Genome/analysis_DarkGenome2weeks_HTO
export SCRNABOX_HOME=/lustre03/project/6070393/COMMON/Dark_Genome/samamiri/pipeline/scrnabox.slurm
export SCRNABOX_PWD=/lustre03/project/6070393/COMMON/Dark_Genome/analysis_DarkGenome2weeks_HTO

bash $SCRNABOX_HOME/launch_scrnabox.sh -h

bash $SCRNABOX_HOME/launch_scrnabox.sh\
-d ${SCRNABOX_PWD} \
--steps 0 \
--method HTO

bash $SCRNABOX_HOME/launch_scrnabox.sh\
-d ${SCRNABOX_PWD} \
--steps 1

bash $SCRNABOX_HOME/launch_scrnabox.sh\
-d ${SCRNABOX_PWD} \
--steps 2


bash $SCRNABOX_HOME/launch_scrnabox.sh\
-d ${SCRNABOX_PWD} \
--steps 3 \
--nFeature_RNA_L 1000 \
--nCount_RNA_U 65000 \
--mitochondria_percent_U 25


bash $SCRNABOX_HOME/launch_scrnabox.sh\
-d ${SCRNABOX_PWD} \
--steps 4 \
--msd T 

bash $SCRNABOX_HOME/launch_scrnabox.sh\
-d ${SCRNABOX_PWD} \
--steps 4 

bash $SCRNABOX_HOME/launch_scrnabox.sh\
-d ${SCRNABOX_PWD} \
--steps 5 

bash $SCRNABOX_HOME/launch_scrnabox.sh\
-d ${SCRNABOX_PWD} \
--steps 6 

bash $SCRNABOX_HOME/launch_scrnabox.sh\
-d ${SCRNABOX_PWD} \
--steps 7 \
--marker T

bash $SCRNABOX_HOME/launch_scrnabox.sh\
-d ${SCRNABOX_PWD} \
--steps 7 \
--fta T

bash $SCRNABOX_HOME/launch_scrnabox.sh\
-d ${SCRNABOX_PWD} \
--steps 7 \
--fta enrich

bash $SCRNABOX_HOME/launch_scrnabox.sh\
-d ${SCRNABOX_PWD} \
--steps 8 \
--dgelist T


bash $SCRNABOX_HOME/launch_scrnabox.sh\
-d ${SCRNABOX_PWD} \
--steps 8 \
--genotype T


bash $SCRNABOX_HOME/launch_scrnabox.sh\
-d ${SCRNABOX_PWD} \
--steps 8 \
--celltype T


