# Data Processing Pipeline
## For Advanced Users Only

**Note:** Workshop participants do NOT need to run this pipeline. All necessary processed files are already provided in `workshop_4/codes/0.2.data_analysis/`. This documentation is provided for reference and reproducibility.

## Dataset Information
**Source:** NCBI BioProject [PRJNA880162](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA880162)

**Study:** Mother-infant microbiome and mycobiome analysis from milk and gut samples

## Setup Instructions

### 1. Download Raw Data
Download samples from NCBI and place them in:
```
workshop_4/data/PRJNA880162/input_data/
```

### 2. Redistribute raw files in 16S and ITS folders
Run the following scripts:
```
workshop_4/codes/0.1.data_processing/Separate 16S and ITS samples.ipynb
```

### 3. Download Reference Database
Follow the instructions in [Workshop 1](https://github.com/DAWGPSU/DAWG_workshops_2025-2026/tree/main/workshop_1) to download the UNITE fungal ITS database.

Place the database files in:
```
workshop_4/data/database/
```

### 4. Install snakemake
Follow the instructions in [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

## Running the Pipeline

The processing pipeline uses Snakemake on a SLURM cluster system. The pipeline must be run twice: once for 16S (bacteria) and once for ITS (fungi).

### Process 16S rRNA Data (Bacteria)
```bash
snakemake --snakefile Snakemodule \
          --configfile config_16S.yaml \
          --cluster-config cluster_config.yaml \
          --cluster "sbatch -o {cluster.log} --partition=open --account=open \
                     --mem={cluster.mem} --time={cluster.time} \
                     --nodes={cluster.nodes} --ntasks={cluster.ntasks}" \
          --use-conda \
          --conda-frontend conda \
          --conda-prefix ../../../conda_envs \
          -p \
          --jobs 60 \
          --cores 10 \
          --latency-wait 60
```

### Process ITS Data (Fungi)
```bash
snakemake --snakefile Snakemodule \
          --configfile config_ITS.yaml \
          --cluster-config cluster_config.yaml \
          --cluster "sbatch -o {cluster.log} --partition=open --account=open \
                     --mem={cluster.mem} --time={cluster.time} \
                     --nodes={cluster.nodes} --ntasks={cluster.ntasks}" \
          --use-conda \
          --conda-frontend conda \
          --conda-prefix ../../../conda_envs \
          -p \
          --jobs 60 \
          --cores 10 \
          --latency-wait 60
```

## Pipeline Outputs
- Phyloseq objects for bacteria (16S) and fungi (ITS)
- ASV tables with taxonomic assignments
- Quality control reports

## System Requirements
- SLURM cluster access
- Conda/Mamba package manager
- Sufficient storage for raw sequencing data and intermediate files

## Pipeline Steps
The Snakemake pipeline follows the DADA2 workflow from Workshop 1:
1. Quality filtering
2. Denoising with DADA2
3. Chimera removal
4. Taxonomic assignment
5. Phyloseq object creation

## Questions?
For pipeline troubleshooting, refer to Workshop 1 documentation or contact the workshop organizers.
