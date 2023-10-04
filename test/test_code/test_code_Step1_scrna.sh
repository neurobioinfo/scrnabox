cd /lustre04/scratch/fiorini9/scrna_pipeline/test_code/pipeline_standard

export SCRNABOX_HOME=/lustre04/scratch/fiorini9/scrna_pipeline/test_code/scrnabox0/scrnabox/scrnabox.slurm
export SCRNABOX_PWD=/lustre04/scratch/fiorini9/scrna_pipeline/test_code/pipeline_standard

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 1
