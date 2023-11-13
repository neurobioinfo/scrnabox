#!/usr/bin/env Rscript

####################
# Chech the R libraries
####################

## load parameters
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
r_list=args[3]
summary_r_check=args[4]
## load library
.libPaths(r_lib_path)
# list_R_packages.ini
    write(paste0("START:"),file=summary_r_check,append=TRUE)
    write(paste0("---"),file=summary_r_check,append=TRUE)
    write(paste0(" "),file=summary_r_check,append=TRUE)
r_list1=read.table(r_list, header = FALSE, sep = "")
i0=0
for (pkg in basename(r_list1$V1)) {
    if (suppressWarnings(suppressMessages(!library(pkg, character.only = TRUE, logical.return = TRUE, quietly = TRUE)))){
        write(paste0("Package <<",pkg,">> is not installed"),file=summary_r_check,append=TRUE)
        i0<-i0+1
    }
}
    write(paste0("---"),file=summary_r_check,append=TRUE)
    write(paste0("END"),file=summary_r_check,append=TRUE)
    if (i0==0) {
     write(paste0("Great news, all the necessary libraries are installed."),file=summary_r_check,append=TRUE)
    } else if (i0>0) {
    write(paste0("Please install the libraries listed above."),file=summary_r_check,append=TRUE)
    }
