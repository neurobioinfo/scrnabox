screen -ls 
screen -S scrnabox014hto

salloc -A def-tdurcan --time=0-1 -c 1 --mem-per-cpu=4g

sacct  -S 2020-10-01  --format="JobID,JobName%30,time,start,end,state"

JOBID 

The configuration files do not exist, the pipeline will create them during execution.


The pipeline, scrnabox.slurm, has the capability to automatically generate the required configuration files if they do not already exist. If the necessary configuration files are not present, the pipeline will create them during execution. This feature ensures that the pipeline can seamlessly proceed even if the configuration files are missing.


cd /lustre03/project/6070393/COMMON/Dark_Genome/analysis_DarkGenome2weeks_HTO
export SCRNABOX_HOME=/home/samamiri/NB043_dge/comm/Dark_Genome/samamiri/pipeline/scrnabox014.slurm
export SCRNABOX_PWD=/lustre03/project/6070393/COMMON/Dark_Genome/analysis_DarkGenome2weeks_HTO


bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 0 \
--method HTO
y

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 1

echo ${pwd_dir}
echo $(basename $(dirname ${pwd_dir}))
echo $OUTPUT_DIR/job_info/logs/$(basename $(dirname ${pwd_dir}))_step1.log


aa=$(pwd)
echo $aa
echo $(dirname $aa)
pwd_dir=$(pwd)
$OUTPUT_DIR/job_info/logs/$(basename $(dirname pwd))_step1.log


bash $SCRNABOX_HOME/launch_scrnabox.sh -h 

seu<-readRDS("/home/samamiri/NB043_dge/comm/Dark_Genome/analysis_DarkGenome2weeks_HTO/step2/objs2/seu1.rds")

sink(paste(output_dir,'/step2/info2/nCount_RNA.txt', sep=""))
seu$nCount_RNA %>% summary
sink()

sink(paste(output_dir,'/step2/info2/nCount_RNA.txt', sep=""))
seu$nCount_RNA %>% summary
sink()

writeLines(seu$nCount_RNA %>% summary,'/lustre03/project/6070393/COMMON/Dark_Genome/analysis_DarkGenome2weeks_HTO/nCount_RNA.txt')
writeLines(capture.output(sessionInfo()), paste(output_dir,'/step2/info2/sessionInfo.txt', sep=""))
seu$nCount_RNA %>% summary
seu$nFeature_RNA %>% summary
Seurat$pt.mito %>% summarys


bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 2

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 2-3 \
--nFeature_RNA_L 1000 \
--nCount_RNA_U 65000 \
--mitochondria_percent_U 25

bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 3 \
--nFeature_RNA_L 1000 \
--nCount_RNA_U 65000 \
--mitochondria_percent_U 25


sh $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 4 \
--msd T 


sh $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 4

sh $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 5 

setwd("/lustre03/project/6070393/COMMON/Dark_Genome/analysis_DarkGenome2weeks_HTO")

seu1<-readRDS("/lustre03/project/6070393/COMMON/Dark_Genome/analysis_DarkGenome2weeks_HTO/step2/objs2/seu1.rds")
DefaultAssay(seu1) <- "HTO"


source(paste(output_dir,'/job_info/parameters/step4_par.txt',sep=""))
library('tidyverse')
dd<-read.csv(paste(output_dir,'/job_info/parameters/step4_antibody_label.txt',sep=''), header=FALSE, sep="")

B0251-TotalSeqB,B0252-TotalSeqB,B0253-TotalSeqB,B0254-TotalSeqB,B0255-TotalSeqB,B0256-TotalSeqB,Doublet,Negative
AIW002,SNCA-A53T,GBA-KO,Parkin-KO,PINK1-KO,SNCA-KO,Doublet,Negative

B0251-TotalSeqB,B0252-TotalSeqB,B0253-TotalSeqB,B0254-TotalSeqB,B0255-TotalSeqB,B0256-TotalSeqB,Doublet,Negative
AIW002,SNCA-A53T,GBA-KO,Parkin-KO,PINK1-KO,SNCA-KO,Doublet,Negative

output_dir=setwd("output_dir")
output_dir<-"/lustre03/project/6070393/COMMON/Dark_Genome/analysis_DarkGenome2weeks_HTO"
dd<-read.csv(paste(output_dir,'/job_info/parameters/step4_antibody_label.txt',sep=''), header=FALSE, sep="")

old_antibody_label=c("B0251-TotalSeqB","B0252-TotalSeqB","B0253-TotalSeqB","B0254-TotalSeqB","B0255-TotalSeqB","B0256-TotalSeqB","Doublet","Negative")
new_antibody_label=c("AIW002","SNCA-A53T","GBA-KO","Parkin-KO","PINK1-KO","SNCA-KO","Doublet","Negative")
sh $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 6 

length(new_antibody_label)

sh $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--marker T


bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 7 \
--fta T

step8_clus_label<-c('DAneurons','DAneurons','other','other','neurons','ependymal','neurons','neurons','neurons','DAneurons','DAneurons','progenitors','Glia','other','progenitors','progenitors','neuroblast','other','neurons','other','other','progenitors','other','DAneurons','neurons','DAneurons','progenitors','neurons','progenitors')

clus_label<-read.csv(paste(output_dir,'/job_info/parameters/step8_clus_label.txt',sep=''), header=FALSE, sep="")
cluster.ids<-str_split(clus_label[1,],",")[[1]]


. /home/samamiri/.bash_profile
backup_zip  /home/samamiri/NB043_dge/comm/Dark_Genome/samamiri/pipeline/scrnabox014.slurm


awk '/Pipestance/'  /home/samamiri/NB043_dge/comm/Dark_Genome/analysis_DarkGenome2weeks_HTO/job_info/logs/LaunchSample1_step1.log
tail -3 /home/samamiri/NB043_dge/comm/Dark_Genome/analysis_DarkGenome2weeks_HTO/job_info/logs/LaunchSample1_step1.log



sh $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--dgelist T

sh $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--celltype T


sh $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 8 \
--genotype T

