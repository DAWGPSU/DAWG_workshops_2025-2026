"""
Snakemake module for DADA2 amplicon sequence variant (ASV) analysis.

This module handles:
- Learning error rates from sequencing data
- Dereplication and denoising with DADA2
- Merging paired-end reads (if applicable)
- Removing chimeric sequences
- Taxonomic assignment
- Creating phyloseq objects at multiple taxonomic levels
"""

import os
import yaml
import pandas as pd

# Configuration variables
conda_dir = config["conda_dir"]
scripts_dir = config["scripts_dir"]

wildcard_constraints:
    study_name = '|'.join(config["study_names"]),
    agglom = 'asv|genus|family'


##############################################
#                FUNCTIONS                   #
##############################################


def get_seqtab_input(wildcards):
    """
    Dynamically determine input for sequence table based on read type.
    
    This function identifies whether the dataset consists of single-end or 
    paired-end sequencing reads:
    - For paired-end: uses merged reads from merge_pairs rule
    - For single-end: uses denoised reads directly from dada rule
    
    Args:
        wildcards: Snakemake wildcards containing study_name
        
    Returns:
        list: Path to appropriate RDS file (mergers.rds or dada_R1.rds)
    """
    # Read strand configuration from checkpoint
    strands_config = checkpoints.dada2_get_strands.get(**wildcards).output[0]
    
    with open(strands_config, 'r') as f:
        strand = f.read().strip()
    
    # If paired-end (R2 present), use merged reads
    if strand == 'R2':
        return expand(
            config["study_dir"] + "processed/dada2/{study_name}/mergers.rds",
            study_name=[wildcards.study_name]
        )
    # If single-end, use R1 dada results directly
    else:
        return expand(
            config["study_dir"] + "processed/dada2/{study_name}/dada_R1.rds",
            study_name=[wildcards.study_name]
        )


##############################################
#                  RULES                     #
##############################################


rule dada2_all:
    """
    Master rule to generate phyloseq objects at multiple taxonomic levels.
    
    Produces phyloseq objects agglomerated at:
    - ASV level (no agglomeration)
    - Genus level
    - Family level
    """
    input:
        phyloseq = expand(
            config["study_dir"] + "processed/phyloseq/{study_name}/phyloseq.rds",
            study_name=config['study_names']
        )


# rule download_microbial_database:
#     """
#     Download reference database from remote URL.
    
#     This rule:
#         - Downloads taxonomic reference databases (e.g., Silva, GTDB)
#         - Saves to centralized database directory
#         - Uses wget for reliable downloading with retry capability
        
#     Database types include:
#         - Silva taxonomic training sets
#         - Species-level assignment databases
        
#     The database URL is specified in config['database_url'] with {database}
#     placeholder that gets filled with the wildcard value.
    
#     Example config:
#         database_url: "https://zenodo.org/record/XXX/files/{database}"
#     """
#     output:
#         config["study_dir"] + "database/%s/{database}"%config['database_name']
#     params:
#         url = lambda w: config['database_url'].format(database=w.database)
#     shell:
#         """
#         wget -O {output} {params.url}
#         """


checkpoint dada2_get_strands:
    """
    Determine if sequencing data is single-end or paired-end.
    
    This checkpoint:
        - Reads the sample table from host removal step
        - Checks for presence of R2 (reverse reads) column
        - Writes 'R2' if paired-end, 'R1' if single-end
        - Enables conditional workflow (merge vs direct to seqtab)
    """
    input:
        config["study_dir"] + "processed/trimmed/{study_name}.csv"
    output:
        config["service_config"] + "dada2/{study_name}.strands.txt"
    run:
        # Read sample table
        sample_table = pd.read_csv(input[0], sep='\t')
        
        # Check if R2 column exists (paired-end indicator)
        with open(output[0], 'w') as f:
            if 'R2' in sample_table.columns.values.tolist():
                f.write('R2')
            else:
                f.write('R1')


rule learn_errors:
    """
    Learn error rates from sequencing data using DADA2.
    
    DADA2's error learning:
        - Models substitution error rates
        - Uses parametric error model
        - Learns from a subset of reads
        - Essential for accurate denoising
        
    Runs separately for forward (R1) and reverse (R2) reads.
    """
    input:
        samples_table = config["study_dir"] + "processed/trimmed/{study_name}.csv"
    output:
        config["study_dir"] + "processed/dada2/{study_name}/err_R{strand}.rds"
    params:
        strand = 'R{strand}'
    conda:
        f"{conda_dir}/dada2.yaml"
    shell:
        f"""
        Rscript {scripts_dir}/dada2/learn_errors.r {{input}} {{output}} {{params}} 
        """


rule derep:
    """
    Dereplicate sequences to reduce computational burden.
    
    Dereplication:
        - Combines identical sequences
        - Retains quality information
        - Significantly reduces dataset size
        - Maintains abundance information
        
    Runs separately for forward (R1) and reverse (R2) reads.
    """
    input:
        samples_table = config["study_dir"] + "processed/trimmed/{study_name}.csv"
    output:
        config["study_dir"] + "processed/dada2/{study_name}/derep_R{strand}.rds"
    params:
        strand = 'R{strand}'
    conda:
        f"{conda_dir}/dada2.yaml"
    shell:
        f"""
        Rscript {scripts_dir}/dada2/derep.r {{input}} {{output}} {{params}}
        """


rule dada:
    """
    Denoise sequences using DADA2 algorithm.
    
    DADA2 denoising:
        - Infers true biological sequences
        - Corrects sequencing errors
        - Resolves sequences differing by single nucleotide
        - Uses learned error rates
        - Outputs amplicon sequence variants (ASVs)
        
    Pooling strategy controlled by config['pool']:
        - FALSE: process samples independently
        - TRUE: pool all samples for rare variant detection
        - pseudo: compromise between speed and sensitivity
    """
    input:
        err_file = config["study_dir"] + "processed/dada2/{study_name}/err_R{strand}.rds",
        derep_file = config["study_dir"] + "processed/dada2/{study_name}/derep_R{strand}.rds"
    output:
        config["study_dir"] + "processed/dada2/{study_name}/dada_R{strand}.rds"
    params:
        params_file = config["params_dir"] + "dada2/{study_name}.yaml"
    conda:
        f"{conda_dir}/dada2.yaml"
    shell:
        f"""
        Rscript {scripts_dir}/dada2/dada.r {{input}}  {{output}} {{params}} 
        """


rule merge_pairs:
    """
    Merge forward and reverse reads for paired-end data.
    
    Merging:
        - Aligns forward and reverse reads
        - Removes pairs with mismatches in overlap region
        - Reconstructs full amplicon sequence
        
    Note:
        Only runs for paired-end data; skipped for single-end via get_seqtab_input()
    """
    input:
        dada_R1 = config['study_dir'] + 'processed/dada2/{study_name}/dada_R1.rds',
        dada_R2 = config['study_dir'] + 'processed/dada2/{study_name}/dada_R2.rds',
        derep_R1 = config['study_dir'] + 'processed/dada2/{study_name}/derep_R1.rds',
        derep_R2 = config['study_dir'] + 'processed/dada2/{study_name}/derep_R2.rds'
    output:
        config["study_dir"] + "processed/dada2/{study_name}/mergers.rds"
    conda:
        f"{conda_dir}/dada2.yaml"
    shell:
        f"""
        Rscript {scripts_dir}/dada2/merge_pairs.r {{input}}  {{output}}
        """


rule make_seqtab:
    """
    Create sequence table (sample x ASV matrix).
    
    The sequence table:
        - Rows: samples
        - Columns: unique ASV sequences
        - Values: abundance counts
        - Input for chimera removal
        
    Input source determined dynamically:
        - Paired-end: merged sequences
        - Single-end: denoised R1 sequences
    """
    input:
        get_seqtab_input
    output:
        config["study_dir"] + "processed/dada2/{study_name}/seqtab.rds"
    conda:
        f"{conda_dir}/dada2.yaml"
    shell:
        f"""
        Rscript {scripts_dir}/dada2/make_seqtab.r {{input}} {{output}}
        """


rule seqtab_nochim:
    """
    Remove chimeric sequences from sequence table.
    
    Chimera detection:
        - Identifies sequences formed from two+ parent sequences
        - Checks if sequence can be reconstructed from more abundant sequences
        - Typically removes 5-20% of sequences but <5% of reads
        
    Essential for accurate diversity estimates and taxonomic assignment.
    """
    input:
        config["study_dir"] + "processed/dada2/{study_name}/seqtab.rds"
    output:
        config["study_dir"] + "processed/dada2/{study_name}/seqtab_nochim.rds"
    conda:
        f"{conda_dir}/dada2.yaml"
    shell:
        f"""
        Rscript {scripts_dir}/dada2/seqtab_nochim.r {{input}} {{output}}
        """


rule assign_taxa:
    """
    Assign taxonomy to ASVs using reference databases.
    
    Taxonomic assignment:
        - Uses naive Bayesian classifier
        - Assigns taxonomy from Kingdom to Species (if possible)
        - Silva database for genus-level assignment
        - Species-level database for exact matches
        
    Databases:
        - First database: Silva reference (general taxonomy)
        - Second database: Species-level reference for exact matching
    """
    input:
        seqtab_nochim = config["study_dir"] + "processed/dada2/{study_name}/seqtab_nochim.rds",
        databases = expand(
            config["study_dir"] + "database/{database_name}/{database}",
            database=config['databases'], database_name=config['database_name']
        )
    output:
        config["study_dir"] + "processed/dada2/{study_name}/taxa.rds"
    params:
        db = lambda w, input: input.databases
    conda:
        f"{conda_dir}/dada2.yaml"
    shell:
        f"""
        Rscript {scripts_dir}/dada2/assign_taxa.r \
            {{input.seqtab_nochim}} \
            {{params}}\
            {{output}}
        """


rule create_phyloseq:
    """
    Convert DADA2 outputs to phyloseq object and export tables.
    
    This rule:
        - Creates phyloseq object from seqtab and taxonomy
        - Exports ASV-level phyloseq object
        - Exports OTU table (ASV abundances)
        - Exports taxonomy table
        
    Output:
        - phyloseq.rds: R object for downstream analysis
        - otu_table.csv: abundance matrix
        - tax_table.csv: taxonomic assignments
    """
    input:
        seqtab_nochim = config["study_dir"] + "processed/dada2/{study_name}/seqtab_nochim.rds",
        taxa = config["study_dir"] + "processed/dada2/{study_name}/taxa.rds",
        metadata = config['metadata']
    output:
        ps = config["study_dir"] + "processed/phyloseq/{study_name}/phyloseq.rds"
    conda:
        f"{conda_dir}/dada2.yaml"
    shell:
        f"""
        Rscript {scripts_dir}/dada2/create_phyloseq.r {{input}} {{output}}
        """

