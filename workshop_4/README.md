# Workshop 4: Integrated Analysis of Microbiome and Mycobiome Data
**Date:** January 29th, 2026

## Overview
This workshop demonstrates integrated analysis techniques for microbiome and mycobiome data, focusing on multivariate statistical approaches and phylogenetic methods.

**Main Analysis Script:** `workshop_4/codes/0.2.data_analysis/1_Main_analysis_PCoA_Procrustes_CCA.ipynb`

## Workshop Content
- **Principal Coordinates Analysis (PCoA)** based on Bray-Curtis distance
- **Phylogenetic tree construction** using 16S rRNA sequences with ghost-tree approach for primary taxonomic annotation
- **UniFrac-based PCoA** incorporating phylogenetic relationships
- **Procrustes analysis** for comparing PCoA ordinations between microbiome and mycobiome datasets
- **Canonical Correlation Analysis (CCA)** for identifying correlated patterns across data types

## Dataset
We use dataset **PRJNA880162** from a mother-infant study examining microbiome and mycobiome composition in milk and gut samples.

**Publication:** [https://doi.org/10.3389/fmicb.2022.1050574](https://doi.org/10.3389/fmicb.2022.1050574)

### Data Processing
**Step 1.** Raw sequencing reads were processed following the DADA2 pipeline described in [Workshop 1](https://github.com/DAWGPSU/DAWG_workshops_2025-2026/tree/main/workshop_1). The complete processing pipeline and instructions are available in `workshop_4/codes/0.1.data_processing`.

**Step 2.** Generation of the phylogenetic trees for microbiome (`MUSCLE` + `fasttree`) and mycobiome data (`ghost-tree`). The corresponding scripts: `workshop_4/codes/0.2.data_analysis/0.1-0.3`

**Processing outputs:**
- Two phyloseq objects (one for bacteria, one for fungi)
- Two phyloseq objects with integrated phylogenetic trees

**Note:** All files necessary for the main analysis script are provided in `workshop_4/codes/0.2.data_analysis/`, so participants **do not need to run** any of the processing steps.

## Installation

### Option 1: Using Conda (Recommended)
```bash
conda env create -f workshop_4/codes/0.2.data_analysis/conda_recipes/r_libraries.yaml
conda activate sp26_workshop
```

### Option 2: Using R
Run the following code in your R console:

```r
# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install("phyloseq")
BiocManager::install("GUniFrac")

# Install CRAN packages
install.packages(c("vegan", "ape", "ggplot2", "ggpubr", 
                   "PMA", "reshape2", "tidyr", "phytools", "compositions"))

# Note: 'parallel' is part of base R

# Verify installations
cat("\n=== Verifying installations ===\n")
libraries <- c("phyloseq", "vegan", "ape", "ggplot2", "ggpubr", 
               "PMA", "parallel", "reshape2", "tidyr", "GUniFrac", "phytools", "compositions")

for (lib in libraries) {
  if (require(lib, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("✓ %s loaded successfully\n", lib))
  } else {
    cat(sprintf("✗ %s failed to load\n", lib))
  }
}

cat("\nInstallation complete!\n")
```

## Getting Started
1. Clone this repository
2. Install required libraries (see Installation section)
3. Navigate to `workshop_4/codes/0.2.data_analysis/`
4. Open `1_Main_analysis_PCoA_Procrustes_CCA.ipynb` in RStudio
5. Run the R Markdown document interactively or knit to generate a report

## Repository Structure
```
workshop_4/
├── README.md
├── codes/
│   ├── 0.1.data_processing/         # Raw data processing pipeline            
│   └── 0.2.data_analysis/
│       ├── 0.1. Building phylogenetic trees. Data pre-processing. R code.ipynb
│       ├── 0.2. Building phylogenetic trees. Ghost-tree usage. Terminal.ipynb
│       ├── 0.3. Building phylogenetic trees. Data post-processing. R code.ipynb
│       ├── 1. Main_analysis_PCoA_Procrustes_CCA.ipynb    # Main workshop script
│       ├── conda_recipes/
│       │   ├── ghost-tree.yaml
│       │   └── r_libraries.yaml
│       └── scripts/
│           └── convert_taxonomy_to_fasta.py
└── data/
    ├── PRJNA880162/
    │   └── processed/phyloseq/       # Processed phyloseq objects
    └── data_analysis/                # Analysis outputs and intermediate files
```

**Key directories:**
- `codes/0.1.data_processing/` - Optional preprocessing pipeline (not needed for workshop)
- `codes/0.2.data_analysis/` - Main workshop scripts and all necessary analysis files
- `data/` - Input data and analysis outputs (created during analysis)

## Questions?
Contact the workshop organizers or open an issue in this repository.
