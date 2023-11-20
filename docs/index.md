# Welcome to scRNAbox's documentation!
ScRNAbox is a single-cell RNA sequencing (scRNAseq) pipeline specifically designed for analyzing data under a High-Performance Computing (HPC) systems using the [Slurm Workload Manager](https://slurm.schedmd.com/). The scRNAbox pipeline incorporates eight analytical steps into a comprehensive scRNAseq analysis that provides the foundation for further investigations. The eight analytical steps are outlined below. 

 <p align="center">
 <img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/bd671431-4a28-4a3a-967a-28f04580a5f8" width="550" height="100">
 </p>

The scRNAbox pipeline provides two distinct, yet highly comparable analysis tracks:

1. **Standard analysis track**
2. **HTO analysis track**

<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/029c737c-c696-4cd3-8a00-3c434cf4e731" height="300"> 

The **standard analysis track** is designed for experiments where each sample is captured and sequenced separately, while the **HTO analysis track** is designed for multiplexed experiments where samples are tagged with sample-specific oligonucleotide tagged Hashtag antibodies (HTO), pooled, and sequenced together. The two tracks share all of the same processes except for Steps 1 and 4. In Step 1, the FASTQ processing for HTO data requires different configurations than the standard analysis. In Step 4, the HTO data is demultiplexed and a different doublet removal method is used. 

For a comprehenseive description of each step, please see the **Pipeline** section of the scRNAbox documentation or see our [pre-print manuscript](https://www.biorxiv.org/content/10.1101/2023.11.13.566851v1). <br/>

For a tutorial that leverages the datasets used as the application cases in our pre-print manuscript, please see [scRNAbox analysis of the midbrain dataset](Dataset1.md).

 - - - -

## Contents
- Pipeline:
    - [Installation](installation.md)
    - [Step 0: Set up](Step0.md)
    - [Step 1: FASTQ to expression matrix](Step1.md)
    - [Step 2: Seurat object and ambient RNA](Step2.md)
    - [Step 3: Quality control and filtering](Step3.md)
    - [Step 4: Doublet removal (standard)](Step4.md)
    - [Step 4: Demultiplexing (HTO)](Step4HTO.md)
    - [Step 5: Integration](Step5.md)
    - [Step 6: Clustering](Step6.md)
    - [Step 7: Cluster annotation](Step7.md)
    - [Step 8: Differential gene expression](Step8.md)
- Documentation:    
    - [Job configurations](config.md)
    - [Execution parameters](reference.md) 
    - [Outputs](outputs.md) 
- Tutorial:
    - [Downloading the Midbrain dataset](midbrain_download.md)
    - [SCRNA analysis track: Midbrain dataset](Dataset1.md)
    - [Downloading the PBMC dataset](pbmc_download.md)
    - [HTO analysis track: PBMC dataset](Dataset2.md)
    - [Analysis of DGE outputs](DEG.md)
    - [Manual CellRanger library preparation](library_prep.md)      
           
- About:
    - [License](LICENSE.md)
    - [Contributing](contributing.md)
    - [Acknowledgement](Acknowledgement.md)