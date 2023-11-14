## Downloading the PBMC dataset
The scRNAseq data produced by [Stoeckius et al.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1603-1) is publicly available in the Gene Expression Omnibus with accession code [GSE108313](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108313). To download the data, we must first install SRAtoolkit (if this is not already installed on your High-Performance Computing (HPC) system). We will create a directory for our raw data and download SRAtoolkit with the following code:

```
mkdir data_download
cd data_download
wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -vxzf sratoolkit.tar.gz
export PATH=$PATH:$PWD/sratoolkit.3.0.5-ubuntu64/bin
```

For more information regarding the SRAtoolkit, please visit the [documentation](https://github.com/ncbi/sra-tools/wiki).

The Sequence Read Archive (SRA) run identifiers for the RNA and antibody assays are:

|Assay|SRR|
|:--|:--|
|RNA|SRR8281306|
|Antibody|SRR8281307|

To download the FASTQ files for the RNA and antibody assays, run the following code. Please note that this may take a very long time. 

```
export PATH=$PATH:$PWD/sratoolkit.3.0.5-ubuntu64/bin
module load StdEnv/2020 gcc/9.3.0
module load sra-toolkit/3.0.0 

#RNA
prefetch SRR8281306 --max-size 100GB
fasterq-dump SRR8281306

#Antibody
prefetch SRR8281307 --max-size 100GB
fasterq-dump SRR8281307 

```
If the FASTQ files for the RNA and antibody assays have been downloaded properly, the `data_download` folder should contain the following:
```
data_download
├── SRR8281306
│   └── SRR8281306.sra
├── SRR8281306_1.fastq
├── SRR8281306_2.fastq
├── SRR8281307
│   └── SRR8281307.sra
├── SRR8281307_1.fastq
└── SRR8281307_2.fastq
```


Next, we will rename the FASTQ files according to the CellRanger nomenclature and transfer the FASTQ files to a folder named `fastqs`. For more information regarding the nomeclature required by the CellRanger _counts_ pipeline, please visit CellRanger's [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/fastq-input).

**Note**: The `fastqs` folder should only contain FASTQ files for the experiment.

```
mkdir fastqs

# RNA assay
cp ~/data_download/SRR8281306_1.fastq  ~/fastqs/run1GEX_S1_L001_R1_001.fastq
cp ~/data_download/SRR8281306_2.fastq  ~/fastqs/run1GEX_S1_L001_R2_001.fastq

# HTO assay
cp ~/data_download/SRR8281307_1.fastq ~/fastqs/run1HTO_S1_L001_R1_001.fastq
cp ~/data_download/SRR8281307_2.fastq ~/fastqs/run1HTO_S1_L001_R2_001.fastq

```
If the above steps were conducted properly, the `fastqs` folder should contain the following files:
```
fastqs
├── run1GEX_S1_L001_R1_001.fastq
├── run1GEX_S1_L001_R2_001.fastq
├── run1HTO_S1_L001_R1_001.fastq
└── run1HTO_S1_L001_R2_001.fastq
```


