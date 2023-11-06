## Manual CellRanger Library Preparation
- [Introduction](#introduction)
- [Standard analysis track](#standard-analysis-track)
    - [library.csv](#librarycsv)
- [HTO analysis track](#hto-analysis-track)
    - [library.csv](#librarycsv)
    - [feature_ref.csv](#feature-ref-csv)


 - - - -

## Introduction
Step 1 (FASTQ to gene expression matrix) of the scRNAbox pipeline leverages the CellRanger _counts_ pipeline to generate gene expression matrices from FASTQ files. In order to run the CellRanger _counts_ pipeline, libraries must be generated to define the information of the FASTQ files. For the [Standard analysis track](#standard-analysis-track), a `library.csv` file must be generated for each sample. For the [HTO analysis track](#hto-analysis-track), a `library.csv` and `feature_ref.csv` file must be generated for each sequencing run. Athough scRNAbox provides an option to automate the library preparation process, it is important that users understand the information that is required for these files. 

In this tutorial we demonstrate how users can manually prepare the `library.csv` and `feature_ref.csv` files. The information presented in this tutorial can also be used to inform the parameters that must be defined by the user for automated library preparation.

 - - - -
## Standard analysis track
#### library.csv
The `library.csv` file defines the necessary information of the FASTQ files for the experiment, including the gene expression and antibody assays. The structure of the `library.csv` file should be: <br />
```
fastqs,sample,library_type
path/to/fastqs/directory/,SampleName,Gene Expression
```
- The `fastqs` column defines the path to the directory that contains the FASTQ files for the experiment. <br /> 
- The `sample` column defines the sample name of the corresponding FASTQ file. Please note that FASTQ files must be named according to standard CellRanger nomenclature (e.g. CTRL1_S1_L001_R1_001.fastq). For more information, please visit CellRanger's [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/fastq-input). <br />
- The `library_type` column defines the assay type. For the Standard Analysis Track, the `library_type` will be "Gene Expression". <br />

#### Example: manual library preparation
For the Standard Analysis Track, if the experiment comprises four samples (two case and two controls), the following steps should be taken for manual library preparation: <br />

1) Navigate to the working directory and create a `samples_info` folder: <br />
```
cd ~/working_directory
mkdir samples_info
```
2) Navigate to the `samples_info` folder and create a folder for each sample run: <br />
```
cd samples_info
mkdir case1
mkdir case2
mkdir control1
mkdir control2
```
3) Navigate to the folder for each sample and create the `library.csv` file. <br />

After performing steps 1-3 above, the structure of the samples_info folder (`~working_directory/samples_info`) for an experiment with four samples should be:
```
working_directory
└── samples_info
    ├── case1
    │   └──library.csv
    ├── case2
    │   └── library.csv
    ├── control1
    │   └── library.csv
    └── control2
        └── library.csv
```


 - - - -

## HTO analysis track
#### library.csv
The `library.csv` file defines the necessary information of the FASTQ files for the experiment, including the gene expression and antibody assays. The structure of the `library.csv` file should be: <br />
```
fastqs,sample,library_type
path/to/fastqs/directory/,SampleNameGEX,Gene Expression
path/to/fastqs/directory/,SampleNameHTO,Antibody Capture
```
- The `fastqs` column defines the path to the directory that contains the FASTQ files for the experiment. <br /> 
- The `sample` column defines the sample name of the corresponding FASTQ file. Please note that FASTQ files must be named according to standard CellRanger nomenclature (e.g. CTRL1_S1_L001_R1_001.fastq). For more information please visit CellRanger's [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/fastq-input). <br />
- The `library_type` column defines the assay type. For the Cell Hashtag Analysis track, each sequencing run should have a "Gene Expression" and "Antibody Capture" assay. For more information, please visit CellRanger's [documentation]("https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis") <br />

#### feature_ref.csv
The `feature_ref.csv` file defines the necessary information for processing the sample-specific barcodes that will eventually be used to demultiplex the pooled samples. For example, if there are four samples pooled together with four unique barcode identifiers, the structure of the `feature_ref.csv` file should be:
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

#### Example: manual library preparation
For the Cell Hashtag Analysis Track, if the experiment comprises four sequencing, the following steps should be taken for manual library preparation: <br />

1) Navigate to the working directory and create a `samples_info` folder: <br />
```
cd ~/working_directory
mkdir samples_info
```
2) Navigate to the `samples_info` folder and create a folder for each sequencing run: <br />
```
cd samples_info
mkdir run1
mkdir run2
mkdir run3
mkdir run4
```
3) Navigate to the folder for each sequencing and create the `library.csv` file. <br />

4) Navigate to the folder for each sequencing and create the `feature_ref.csv` file. <br />

After performing steps 1-4 above, the structure of the samples_info folder (`~working_directory/samples_info`) for an experiment with four sequencing runs should be:
```
working_directory
├── samples_info
    ├── run1
    │    ├── library.csv
    │    └── feature_ref.csv
    ├── run2
    │    ├── library.csv
    │    └── feature_ref.csv
    ├── run3
    │    ├── library.csv
    │    └── feature_ref.csv   
    └── run4
        ├── library.csv
        └── feature_ref.csv   
```

