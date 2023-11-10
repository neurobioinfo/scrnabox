# Step 1: FASTQ to gene expression matrix
- [Library preparation](#library-preparation)
- [Running Step 1](#running-step-1)
 - - - -
In Step 1, gene expression matrices are generated from FASTQ files using the CellRanger _counts_ pipeline. Prior to running CellRanger, `library.csv` files must be prepared to define the FASTQ files for each sample. In addition, `feature_ref.csv` files must be prepared for the HTO track to define the HTOs used in the experiment. ScRNAbox provides an option for automating library preparation but the correct information must still be defined in the parameters. Alternatively, users may manually prepapre the libraries. For a tutorial on manual library preparation, please see the [tutorial](library_prep.md).<br />
 - - - -
## Library preparation
#### library.csv
The `library.csv` file defines the necessary information of the FASTQ files for the experiment, including the gene expression and antibody assays. The structure of the `library.csv` file should be: <br />
```
fastqs,sample,library_type
path/to/fastqs/directory/,SampleNameGEX,Gene Expression
path/to/fastqs/directory/,SampleNameHTO,Antibody Capture
```
- The `fastqs` column defines the path to the directory that contains the FASTQ files for the experiment. <br /> 
- The `sample` column defines the sample name of the corresponding FASTQ file. Please note that FASTQ files must be named according to standard CellRanger nomenclature (e.g. CTRL1_S1_L001_R1_001.fastq). For more information please visit CellRanger's [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/fastq-input). <br />
- The `library_type` column defines the assay type. For the standard analysis track, this will always be "Gene expression". For the HTO analysis track, each sequencing run should have a "Gene Expression" and "Antibody Capture" assay. For more information, please visit CellRanger's [documentation]("https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis") <br />

For more information regarding the preparation of the `library.csv`, please see CellRanger's [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis).

#### feature_ref.csv
The `feature_ref.csv` file defines the necessary information for processing HTOs that will eventually be used to demultiplex the pooled samples. For example, if there are four samples pooled together with four unique HTOs, the structure of the `feature_ref.csv` file should be:
```
id,name,read,pattern,sequence,feature_type
Hash1,B0251_TotalSeqB,R2,5PNNNNNNNNNN(BC),GTCAACTCTTTAGCG,Antibody Capture
Hash2,B0252_TotalSeqB,R2,5PNNNNNNNNNN(BC),TGATGGCCTATTGGG,Antibody Capture
Hash3,B0253_TotalSeqB,R2,5PNNNNNNNNNN(BC),TTCCGCCTCTCTTTG,Antibody Capture
Hash4,B0254_TotalSeqB,R2,5PNNNNNNNNNN(BC),AGTAAGTTCAGCGTA,Antibody Capture
```
- The `id` column defines the barcode ID which will be used to track the feature counts. <br /> 
- The `name` column defines the arbitrary name for the barcode identifier. <br /> 
- The `read` column defines which RNA sequencing read contains the barcode sequence. This value Will be either R1 or R2.<br /> 
- The `pattern` column defines the pattern of the barcode identifiers. For more information please visit the 10X Genomics [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis#pattern)<br /> 
- The `sequence` column defines nucleotide sequence associated with the barcode identifier.<br /> 
- The `feature_type` column defines the type of feature used for sample identification. Please ensure that the feature_type in the `feature_ref.csv` file matches a library_type in the `library.csv` file.  <br /> 

For more information regarding the preparation of the `feature_ref.csv`, please see CellRanger's [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis).
 - - - -
 
## Running Step 1
 
The following parameters are adjustable for Step 1 of the **standard track** (`~/working_directory/job_info/parameters/step1_par.txt`):

|Parameter|Default|Description|
|:--|:--|:--|
|par_automated_library_prep|No| Whether or not to perform automated library prep. Alternatively, you may set this parameter to "no" and manually prepare the libraries.|
|par_fastq_directory|NULL|Path to directory containing the FASTQ files. This directory should only contain FASTQ files for the experiment.|
|par_sample_names|NULL|The sample names used to name the FASTQ files according to CellRanger nomeclature|
|par_rename_samples|Yes| Whether or not you want to rename your samples. These names will be used to identify cells in the Seurat objects|
|par_new_sample_names|NULL| New sample names. Make sure they are defined in the same order as 'par_sample_names'|
|par_paired_end_seq|Yes| Whether or not paired-end sequencing was performed|
|par_ref_dir_grch|NULL|Path to reference genome for FASTQ alignment. 10X Genomics reference genomes are available for download. For more information see the 10X Genomics [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ct).|
|par_r1_length|NULL|Minimum number of bases to retain for R1 sequence of gene expression|
|par_r2_length|NULL|Minimum number of bases to retain for R2 sequence of gene expression|
|par_mempercode|30|For clusters whose job managers do not support memory requests, it is possible to request memory in the form of cores. This option will scale up the number of threads requested via the __MRO_THREADS__ variable according to how much memory a stage requires when given to the ratio of memory on your nodes.|
|par_include_introns|No|Whether or not to include intronic reads in the gene expression matrix|
|par_no_target_umi_filter|No| Whether or not to tirn of CellRanger's target UMI filtering subpipeline|
|par_expect_cells|NULL| Expected number of cells. By default, CellRanger's auto-estimate algorithm will be used|
|par_force_cells|NULL| Force the CellRanger count ipeline to use N cells.|
|par_no_bam|No| Whether or not to skip the bam file generation in the CellRanger pipeline.|

The following parameters are adjustable for Step 1 of the **HTO track** (`~/working_directory/job_info/parameters/step1_par.txt`):

|Parameter|Default|Description|
|:--|:--|:--|
|par_automated_library_prep|Yes|Whether or not to perform automated library prep. Alternatively, you may set this parameter to "no" and manually prepare the libraries.|
|par_fastq_directory|NULL|Path to directory containing the FASTQ files. This directory should only contain FASTQ files for the experiment.|
|par_RNA_run_names|NULL|The names of the sequencing runs for the RNA assay|
|par_HTO_run_names|NULL|The names of the sequencing runs for the HTO assay|
|par_seq_run_names|NULL|The user-selected name for the sequencing run.  These names will be used to identify cells in the Seurat objects|
|par_paired_end_seq|Yes|Whether or not paired-end sequencing was performed|
|id|NULL|Barcode ID which will be used to track the feature counts|
|name|NULL|The user-selected name for the barcode identifier|
|read|R2|Which RNA sequencing read contains the barcode sequence. This value Will be either R1 or R2.|
|pattern|NULL|The pattern of the barcode identifiers|
|sequence|NULL|The nucleotide sequence associated with the barcode identifier|
|par_ref_dir_grch|NULL|Path to reference genome for FASTQ alignment. 10X Genomics reference genomes are available for download. For more information see the 10X Genomics [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ct).|
|par_r1_length|NULL|Minimum number of bases to retain for R1 sequence of gene expression|
|par_r2_length|NULL|Minimum number of bases to retain for R2 sequence of gene expression|
|par_mempercode|30|For clusters whose job managers do not support memory requests, it is possible to request memory in the form of cores. This option will scale up the number of threads requested via the __MRO_THREADS__ variable according to how much memory a stage requires when given to the ratio of memory on your nodes.|
|par_include_introns|No|Whether or not to include intronic reads in the gene expression matrix|
|par_no_target_umi_filter|No| Whether or not to tirn of CellRanger's target UMI filtering subpipeline|
|par_expect_cells|NULL| Expected number of cells. By default, CellRanger's auto-estimate algorithm will be used|
|par_force_cells|NULL| Force the CellRanger count ipeline to use N cells.|
|par_no_bam|No| Whether or not to skip the bam file generation in the CellRanger pipeline.|

Given that CellRanger runs a user interface, it is recommended to run Step 1 in a **'screen'** which will allow the the task to keep running if the connection is broken. To run Step 1, use the following command:
```
screen -S run_scrnabox
bash $SCRNABOX_HOME/launch_scrnabox.sh \
-d ${SCRNABOX_PWD} \
--steps 1
```

The resulting output files are deposited into `~/working_directory/step1`. The expression filtered matrix, features, and barcode files outputed by CellRanger are located in `~/working_directory/step1/sample/ouput_folder/outs/filtered_feature_bc_matrix`.

**Note:** If CellRanger was successfull it will display _Pipestance completed successfully!_ If this message is not displayed, you must re-run the Step 1. If CellRanger fails a second time, please contact the developers.

**Note:** If you do not have access to FASTQ files for your experiment, you may intiate the pipeline at which ever Analytical Step takes your data object as input. In the case where FASTQ files are not available, users do not have to create the `samples_info` folder. For more information see [Processed Data](PROC.md).  



