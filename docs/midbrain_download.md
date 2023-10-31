## Downloading the midbrain dataset
The scRNAseq data produced by [Smajic et al.](https://academic.oup.com/brain/article/145/3/964/6469020) is publicly available in the Gene Expression Omnibus with accession code [GSE157783](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157783). To download the data, we must first install SRAtoolkit (if this is not already installed on your High-Performance Computing (HPC) system). We will create a directory for our raw data and download SRAtoolkit with the following code:

```
mkdir data_download
cd data_download
wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -vxzf sratoolkit.tar.gz
export PATH=$PATH:$PWD/sratoolkit.3.0.5-ubuntu64/bin
```

For more information regarding the SRAtoolkit, please visit the [documentation](https://github.com/ncbi/sra-tools/wiki).

The Sequence Read Archive (SRA) run identifiers for each of the 11 samples in the midbrain dataset are:

|Sample|SRR|
|:--|:--|
|PD1|SRR12621862|
|PD2|SRR12621863|
|PD3|SRR12621864|
|PD4|SRR12621865|
|PD5|SRR12621866|
|CTRL1|SRR12621867|
|CTRL2|SRR12621868|
|CTRL3|SRR12621869|
|CTRL4|SRR12621870|
|CTRL5|SRR12621871|
|CTRL6|SRR12621872|

**Note**: If you simply want to test scRNAbox's Standard scRNAseq Analysis Track, it may be best to only incorportate a subset of samples in a test run, as using all 11 samples will take substantially longer. In this case, we suggest including at least three PD sample and three control to facilitate differential gene expression (DGE) contrasts in Step 8. 

To download the FASTQ files for all 11 samples, run the following code. Please note that this may take a very long time. 

```
export PATH=$PATH:$PWD/sratoolkit.3.0.5-ubuntu64/bin
module load StdEnv/2020 gcc/9.3.0
module load sra-toolkit/3.0.0 

#PD1
prefetch SRR12621862 
fasterq-dump SRR12621862 

#PD2
prefetch SRR12621863 
fasterq-dump SRR12621863  

#PD3
prefetch SRR12621864 
fasterq-dump SRR12621864  

#PD4
prefetch SRR12621865 
fasterq-dump SRR12621865  

#PD5
prefetch SRR12621866 
fasterq-dump SRR12621866 

#CTRL1
prefetch SRR12621867 
fasterq-dump SRR12621867 

#CTRL2
prefetch SRR12621868 
fasterq-dump SRR12621868  

#CTRL3
prefetch SRR12621869 
fasterq-dump SRR12621869  

#CTRL4
prefetch SRR12621870 
fasterq-dump SRR12621870  

#CTRL5
prefetch SRR12621871 
fasterq-dump SRR12621871  

#CTRL6
prefetch SRR12621872
fasterq-dump SRR12621872 
```
If the FASTQ files for all 11 samples have been downloaded properly, the `data_download` folder should contain the following:
```
data_download
├── SRR12621862
│   └── SRR12621862.sra
├── SRR12621862_1.fastq
├── SRR12621862_2.fastq
├── SRR12621863
│   └── SRR12621863.sra
├── SRR12621863_1.fastq
├── SRR12621863_2.fastq
├── SRR12621864
│   └── SRR12621864.sra
├── SRR12621864_1.fastq
├── SRR12621864_2.fastq
├── SRR12621865
│   └── SRR12621865.sra
├── SRR12621865_1.fastq
├── SRR12621865_2.fastq
├── SRR12621866
│   └── SRR12621866.sra
├── SRR12621866_1.fastq
├── SRR12621866_2.fastq
├── SRR12621867
│   └── SRR12621867.sra
├── SRR12621867_1.fastq
├── SRR12621867_2.fastq
├── SRR12621868
│   └── SRR12621868.sra
├── SRR12621868_1.fastq
├── SRR12621868_2.fastq
├── SRR12621869
│   └── SRR12621869.sra
├── SRR12621869_1.fastq
├── SRR12621869_2.fastq
├── SRR12621870
│   └── SRR12621870.sra
├── SRR12621870_1.fastq
├── SRR12621870_2.fastq
├── SRR12621871
│   └── SRR12621871.sra
├── SRR12621871_1.fastq
├── SRR12621871_2.fastq
├── SRR12621872
│   └── SRR12621872.sra
├── SRR12621872_1.fastq
└── SRR12621872_2.fastq
```


Next, we will rename the FASTQ files according to the CellRanger nomenclature and transfer the FASTQ files to a folder named `fastqs`. For more information regarding the nomeclature required by the CellRanger _counts_ pipeline, please visit CellRanger's [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/fastq-input).

**Note**: The `fastqs` folder should only contain FASTQ files for the experiment.

```
mkdir fastqs

#PD1
cp ~/data_download/SRR12621862_1.fastq ~/fastqs/PD1_S1_L001_R1_001.fastq
cp ~/data_download/SRR12621862_2.fastq ~/fastqs/PD1_S1_L001_R2_001.fastq

#PD2
cp ~/data_download/SRR12621863_1.fastq ~/fastqs/PD2_S1_L001_R1_001.fastq
cp ~/data_download/SRR12621863_2.fastq ~/fastqs/PD2_S1_L001_R2_001.fastq

#PD3
cp ~/data_download/SRR12621864_1.fastq ~/fastqs/PD3_S1_L001_R1_001.fastq
cp ~/data_download/SRR12621864_2.fastq ~/fastqs/PD3_S1_L001_R2_001.fastq

#PD4
cp ~/data_download/SRR12621865_1.fastq ~/fastqs/PD4_S1_L001_R1_001.fastq
cp ~/data_download/SRR12621865_2.fastq ~/fastqs/PD4_S1_L001_R2_001.fastq

#PD5
cp ~/data_download/SRR12621866_1.fastq ~/fastqs/PD5_S1_L001_R1_001.fastq
cp ~/data_download/SRR12621866_2.fastq ~/fastqs/PD5_S1_L001_R2_001.fastq

#Ctrl1
cp ~/data_download/SRR12621867_1.fastq ~/fastqs/CTRL1_S1_L001_R1_001.fastq
cp ~/data_download/SRR12621867_2.fastq ~/fastqs/CTRL1_S1_L001_R2_001.fastq

#Ctrl2
cp ~/data_download/SRR12621868_1.fastq ~/fastqs/CTRL2_S1_L001_R1_001.fastq
cp ~/data_download/SRR12621868_2.fastq ~/fastqs/CTRL2_S1_L001_R2_001.fastq

#Ctrl3
cp ~/data_download/SRR12621869_1.fastq ~/fastqs/CTRL3_S1_L001_R1_001.fastq
cp ~/data_download/SRR12621869_2.fastq ~/fastqs/CTRL3_S1_L001_R2_001.fastq

#Ctrl4
cp ~/data_download/SRR12621870_1.fastq ~/fastqs/CTRL4_S1_L001_R1_001.fastq
cp ~/data_download/SRR12621870_2.fastq ~/fastqs/CTRL4_S1_L001_R2_001.fastq

#Ctrl5
cp ~/data_download/SRR12621871_1.fastq ~/fastqs/CTRL5_S1_L001_R1_001.fastq
cp ~/data_download/SRR12621871_2.fastq ~/fastqs/CTRL5_S1_L001_R2_001.fastq

#Ctrl6
cp ~/data_download/SRR12621872_1.fastq ~/fastqs/CTRL6_S1_L001_R1_001.fastq
cp ~/data_download/SRR12621872_2.fastq ~/fastqs/CTRL6_S1_L001_R2_001.fastq
```
If the above steps were conducted properly, the `fastqs` folder should contain the following files:
```
├── CTRL1_S1_L001_R1_001.fastq
├── CTRL1_S1_L001_R2_001.fastq
├── CTRL2_S1_L001_R1_001.fastq
├── CTRL2_S1_L001_R2_001.fastq
├── CTRL3_S1_L001_R1_001.fastq
├── CTRL3_S1_L001_R2_001.fastq
├── CTRL4_S1_L001_R1_001.fastq
├── CTRL4_S1_L001_R2_001.fastq
├── CTRL5_S1_L001_R1_001.fastq
├── CTRL5_S1_L001_R2_001.fastq
├── CTRL6_S1_L001_R1_001.fastq
├── CTRL6_S1_L001_R2_001.fastq
├── PD1_S1_L001_R1_001.fastq
├── PD1_S1_L001_R2_001.fastq
├── PD2_S1_L001_R1_001.fastq
├── PD2_S1_L001_R2_001.fastq
├── PD3_S1_L001_R1_001.fastq
├── PD3_S1_L001_R2_001.fastq
├── PD4_S1_L001_R1_001.fastq
├── PD4_S1_L001_R2_001.fastq
├── PD5_S1_L001_R1_001.fastq
└── PD5_S1_L001_R2_001.fastq
```
