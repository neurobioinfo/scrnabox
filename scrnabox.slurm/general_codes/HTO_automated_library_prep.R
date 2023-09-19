##########################################
# Cell Hashtag Analysis Track library prep (Step 1)
########################################## 

## load parameters
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]
pipeline_home=args[3]

## load libraries
.libPaths(r_lib_path)

## load parameters
source(paste(output_dir,'/job_info/parameters/step1_par.txt',sep=""))

## automated library prep 
if (tolower(par_automated_library_prep)=='yes') {
    dir.create(paste(output_dir, "/samples_info", sep = ""))
    new_dir <- paste(output_dir, "/samples_info/", sep = "")

for(i in length(par_seq_run_names)) {
    #create directory
    dir.create(paste(new_dir,par_seq_run_names[i], sep = ""))

    #write the library.csv file 
    fastqs <- c(par_fastq_directory, par_fastq_directory)
    sample <- c(par_RNA_run_names[i], par_HTO_run_names[i])
    library_type <- c("Gene Expression","Antibody Capture") 
    df <- data.frame(fastqs, sample,library_type)
    write.table(df, file = paste(new_dir,par_seq_run_names[i],"/library.csv", sep=""),sep=",", row.names = FALSE, quote = FALSE) 

    #write feature_ref.csv files 
    length <- length(id)
    pattern <- rep.int(pattern, length)
    read <- rep.int(read, length)
    feature_type <- "Antibody Capture"
    feature_type <- rep.int(feature_type, length)
    df <- data.frame(id,name,read, pattern, sequence, feature_type)
    write.table(df, file = paste(new_dir,par_seq_run_names[i],"/feature_ref.csv", sep=""),sep=",", row.names = FALSE, quote = FALSE) 
}   
}

if (tolower(par_automated_library_prep)=='no') {
print("Skipping automated library prep for CellRanger")  
}



