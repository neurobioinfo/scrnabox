# Load R
module load r/4.2.1
# Folder for R packages 
R_PATH=~/scratch/R/x86_64-pc-linux-gnu-library/4.2
mkdir -p $R_PATH
# Install package
Rscript ./scrnabox.slurm/soft/R/install_packages_scrnabox.R $R_PATH
