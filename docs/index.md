# Welcome to scRNAbox's documentation!
ScRNAbox is a single-cell RNA sequencing (scRNAseq) pipeline specifically designed for analyzing data under a High-Performance Computing (HPC) systems using the [Slurm Workload Manager](https://slurm.schedmd.com/). The scRNAbox pipeline incorporates eight analytical steps into a comprehensive scRNAseq analysis that provides the foundation for further investigations. The eight analytical steps are outlined below. 

 <p align="center">
 <img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/3221e078-68d4-4fde-8a75-9d36691c1cf5" width="550" height="100">
 </p>

The scRNAbox pipeline provides two distinct, yet highly comparable analysis tracks:

1. **Standard analysis track**
2. **HTO analysis track**

The **standard analysis track** is designed for experiments where each sample is captured and sequenced separately, while the **HTO analysis track** is designed for multiplexed experiments where samples are tagged with sample-specific oligonucleotide tagged Hashtag antibodies (HTO), pooled, and sequenced together. The **HTO analysis track** is distinguished by an additional sample demultiplexing step that assigns cells to their sample-of-origin via the sample-specific HTOs. 

<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/3a6df83e-e104-45d2-9b04-fe246642c6a8" height="300"> 

For a comprehenseive description of each step, please see the **Pipeline** section of the scRNAbox documentation or see our [pre-print manuscript](). <br/>

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
- Additional information:    
    - [Job configurations](config.md)
    - [Execution parameters](reference.md) 
    - [Outputs](outputs.md) 
- Tutorial:
    - [Downloading the midbrain dataset](midbrain_download.md)
    - [scRNAbox analysis of midbrain dataset](Dataset1.md)
    - [Analysis of DEGs](DEG.md)
    - [Manual CellRanger library preparation](library_prep.md)      
           
- About:
    - [License](LICENSE.md)
    - [Changelog](changelog.md)
    - [Contributing](contributing.md)
    - [Acknowledgement](Acknowledgement.md)