cd /lustre04/scratch/fiorini9/scrna_pipeline/standard_v1.39/pipeline

export SCRNABOX_HOME=/lustre04/scratch/fiorini9/scrna_pipeline/standard_v1.39/scrnabox0/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=/lustre04/scratch/fiorini9/scrna_pipeline/standard_v1.39/pipeline

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 1
