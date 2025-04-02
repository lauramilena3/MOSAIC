**MOSAIC** is a modularized Snakemake workflow designed for reproducible viral metagenomics. It integrates several bioinformatics tools for processing and analyzing viral sequencing data, enabling efficient, scalable, and reproducible workflows for metagenomic studies.

## Requirements

- **Conda**: A package manager for managing environments and dependencies
- **Mamba**: A fast, drop-in replacement for Conda, also used for managing environments and dependencies. It provides faster package resolution and environment creation compared to Conda.

You can learn how to install it here:

[Installing Mamba/Anaconda](https://www.notion.so/Installing-Mamba-Anaconda-15f1465d41b780c7be0bfa9b13fbf605?pvs=21)

## Installation

### 1. Clone the Repository

Create a dedicated folder for your bioinformatics applications to keep them organized. Navigate to this folder before cloning the repository.

```bash
mkdir ~/apps
cd ~/apps
```

Clone the MOSAIC repository to your local machine:

```bash
git clone https://github.com/lauramilena3/MOSAIC
cd MOSAIC/mosaic
```

### 2. Create Conda Environment

Create a conda environment with the required dependencies by running:

```bash
mamba env create -n Mosaic -f Mosaic.yaml
```

Activate the environment:

```bash
conda activate Mosaic
```
