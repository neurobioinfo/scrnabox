# Welcome to scRNAbox's documentation!
ScRNAbox is a single-cell RNA sequencing (scRNAseq) pipeline specifically designed for analyzing data under a High-Performance Computing (HPC) systems using the [Slurm Workload Manager](https://slurm.schedmd.com/). ScRNAbox provides two distinct, yet highly comparable Analysis Tracks:

1. **Standard scRNAseq**
2. **Cell Hashtag scRNAseq**

The Standard Analysis Track is designed for experiments where each sample is captured and sequenced separately, while the Cell Hashtag Analysis Track is designed for multiplexed experiments, whereby samples are tagged with sample-specific barcodes, pooled, and sequenced together. The Cell Hashtag Analysis Track is distinguished by an additional sample demultiplexing Step that assigns cells to their sample-of-origin via the sample-specific barcodes. 

<img src="https://github.com/neurobioinfo/scrnabox/assets/110110777/3a6df83e-e104-45d2-9b04-fe246642c6a8" height="300"> 

For instructions on how to run each Analytical Step of the [Standard scRNAseq](SCRNA.md) and [Cell Hashtag scRNAseq](HTO.md) Analysis Track please see the respective tutorials. For a demonstration that leverages the datasets used as the application cases in the manuscript please see [Dataset1: Smajic et al.](Dataset1.md) and [Datset2: Stoeckius et al.](Dataset2.md) for the Standard scRNAseq and Cell Hashtag scRNAseq Analysis Track, respectively.

## Contents
- [Installation](installation.md)
- [Tutorial:]()
    - [Standard scRNAseq](SCRNA.md)
    - [Cell Hashtag scRNAseq](HTO.md)
    - [Processed Data](PROC.md)
    - [Dataset1: Smajic et al.](Dataset1.md)
    - [Datset2: Stoeckius et al.](Dataset2.md)
- [FAQ](FAQ.md)
- [Reference](reference.md)
