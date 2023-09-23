# Welcome to scRNAbox's documentation!
ScRNAbox is a single-cell RNA sequencing (scRNAseq) pipeline specifically designed for analyzing data under a High-Performance Computing (HPC) systems using the [Slurm Workload Manager](https://slurm.schedmd.com/). The scRNAbox pipeline incorporates nine Analytical Steps into a comprehensive scRNAseq analysis and provides the foundation for further investigations. The nine Analytical Steps are outlined below. 

<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/eccddd8e-4ea2-4c1e-9427-8ba40e6418ba" width="550" height="100">

The scRNAbox pipeline provides two distinct, yet highly comparable Analysis Tracks:

1. **Standard scRNAseq**
2. **Cell Hashtag scRNAseq**

The **Standard Analysis Track** is designed for experiments where each sample is captured and sequenced separately, while the **Cell Hashtag Analysis Track** is designed for multiplexed experiments, whereby samples are tagged with sample-specific barcodes, pooled, and sequenced together. The Cell Hashtag Analysis Track is distinguished by an additional sample demultiplexing Step that assigns cells to their sample-of-origin via the sample-specific barcodes. 

<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/3a6df83e-e104-45d2-9b04-fe246642c6a8" height="300"> 

For a comprehenseive description of each Analytical Step, please see [Standard Analysis Track](SCRNA.md) and [Cell Hashtag Analysis Track](HTO.md). <br/>

For a tutorial that leverages the datasets used as the application cases in our pre-print manuscript, please see [Standard Analysis: Midbrain dataset](Dataset1.md) and [Cell Hashtag Analysis: PBMC dataset](Dataset2.md).

 - - - -

## Contents
- [Installation](installation.md)
- Overview:
    - [Standard Analysis Track](SCRNA.md)
    - [Cell Hashtag Analysis Track](HTO.md)
    - [Execution parameters](reference.md)
    - [Outputs](outputs.md)
- Tutorial   
    - [Standard Analysis Track: Midbrain dataset](Dataset1.md)
    - [Cell Hashtag Analysis Track: PBMC dataset](Dataset2.md)
    - [Processed Data](PROC.md)
- [FAQ](FAQ.md)
