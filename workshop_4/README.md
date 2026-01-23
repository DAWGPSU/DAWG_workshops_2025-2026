# Workshop 4. January 29th, 2026.
# Integrated analysis of microbiome and mycobiome data
Main workshop script file: `workshop_4/codes/0.2.data_analysis/1_Main_analysis_PCoA_Procrustes_CCA`.
Workshop content:
- PCoA based on Bray-Curtis distance
- Building phylogenetic tree based on 16S rRNA and using primary taxanomic annotation (ghost-tree) (workshop_4/codes/0.2.data_analysis/1_Main_analysis_PCoA_Procrustes_CCA)
- PCoA based on Unifrac metric
- Comparing PCoAs of two datatypes using Procrustes
- Canonical correlation analysis.

The dataset (PRJNA880162) used in the workshop consists of the microbiome and mycobiome sequences extracted from milk and gut samples of the mother-infant study (https://doi.org/10.3389/fmicb.2022.1050574).
The raw reads were processed according the instructions provided in the [workshop_1](https://github.com/DAWGPSU/DAWG_workshops_2025-2026/tree/main/workshop_1). The processing pipelines and instructions could be found in `workshop_4/codes/0.1.data_processing`.
The data processing pipeline generates two phyloseq files containing microbial and mycobial data. The phyloseq files were then processed using scripts 0.1-0.3, that add phylogenetic trees to the phyloseq objects. 
All the files, neccesary for the work of the main analysis script are provided in this repository (`workshop_4/codes/0.2.data_analysis/1_Main_analysis_PCoA_Procrustes_CCA`).  

### Installation
In order to run the main script file (1_Main_analysis_PCoA_Procrustes_CCA), please install all neccesary libraries in advance, either through conda envrironment: `conda env create -f "workshop_4/0.2.data_analysis/conda_recipes/r_libraries.yaml"` or using R installation code:
```
# Please run the following code to install all neccesary libraries, 
# unless you have already installed them using conda recipe: conda_recipes/r_libraries.yaml.

# R Library Installation Script
# Install required packages for microbiome analysis

# Install BiocManager if not already installed (needed for Bioconductor packages)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install("phyloseq")
BiocManager::install("GUniFrac")

# Install CRAN packages
install.packages("vegan")
install.packages("ape")
install.packages("ggplot2")
install.packages("ggpubr")
install.packages("PMA")
install.packages("reshape2")
install.packages("tidyr")
install.packages("phytools")

# Note: 'parallel' is part of base R and doesn't need to be installed

# Verify installations by loading libraries
cat("\n=== Verifying installations ===\n")
libraries <- c("phyloseq", "vegan", "ape", "ggplot2", "ggpubr", 
               "PMA", "parallel", "reshape2", "tidyr", "GUniFrac", "phytools")

for (lib in libraries) {
  if (require(lib, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("✓ %s loaded successfully\n", lib))
  } else {
    cat(sprintf("✗ %s failed to load\n", lib))
  }
}

cat("\nInstallation complete!\n")
```

