"""
Snakemake module for filtering and trimming sequencing reads with DADA2.

This module handles:
- Organizing samples by batch
- Creating per-batch parameter files
- Running DADA2 filterAndTrim on single-end and paired-end reads
- Generating sample tables with trimming statistics
"""

import os
import yaml
import pandas as pd
import numpy as np

# Configuration variables
conda_dir = config["conda_dir"]
scripts_dir = config["scripts_dir"]

wildcard_constraints:
    study_name = '|'.join(config["study_names"])


##############################################
#                FUNCTIONS                   #
##############################################


def get_trimming_logs(wildcards):
    """
    Generate list of trimming log files for all samples in a study.
    
    Args:
        wildcards: Snakemake wildcards containing study_name
        
    Returns:
        list: Paths to log files for all samples (using last strand only)
        
    Note:
        For paired-end data, only requests the _2.log file (which implies _1 completed)
    """
    config_batches = checkpoints.filter_and_trim_samples.get(**wildcards).output[0]
    with open(config_batches, "r") as file:
        config_batches_data = yaml.safe_load(file)
    
    # Get the last strand (2 for paired-end, 1 for single-end)
    strands = [config_batches_data['strands'].split(";")[-1]]
    
    # Collect all sample names across batches
    sample_names = []
    for batch, samples in config_batches_data['samples'].items():
        sample_names.extend(samples)
    
    return expand(
        config["study_dir"] + "processed/trimmed/filter_and_trim/{study_name}/{sample}/{sample}_{strand}.log",
        study_name=[wildcards.study_name],
        sample=sample_names,
        strand=strands
    )


def form_trimmomatic_params(study_name):
    """
    Construct Trimmomatic parameter string from YAML configuration.
    
    Args:
        study_name: Name of the study
        
    Returns:
        str: Formatted parameter string for Trimmomatic
    """
    param_file = rules.trimmomatic.input.params_file.format(study_name=study_name)
    
    with open(param_file) as f:
        params = yaml.safe_load(f)
    
    trimmomatic_params = []
    for param, values in params.items():
        param_str = [param] + [str(values[key]) for key in values['params_order']]
        param_str = ':'.join(param_str)
        trimmomatic_params.append(param_str)
    
    return ' '.join(trimmomatic_params)


##############################################
#                  RULES                     #
##############################################

rule dada_filter_and_trim_all:
    """
    Master rule to filter and trim all samples across all studies.
    
    Ensures all studies have completed trimming and generated sample tables.
    """
    input:
        config_expansion = expand(
            config["study_dir"] + "processed/trimmed/{study_name}.csv",
            study_name=config["study_names"]
        )


checkpoint filter_and_trim_samples:
    """
    Create YAML configuration mapping samples to batches and strand information.
    
    This checkpoint:
        - Reads batch information from batches.csv
        - Extracts strand configuration (single vs paired-end)
        - Groups samples by batch
        - Creates trigger file for dynamic batch processing
        
    Output structure:
        strands: "1" or "1;2"
        batches:
            batch1: [sample1, sample2, ...]
            batch2: [sample3, sample4, ...]
    """
    input:
        config["study_dir"] + "{study_name}.batch.csv"
    output:
        config["service_config"] + "filter_and_trim/{study_name}.yaml"
    run:
        # Read batch information
        data = pd.read_csv(input[0])
        
        # Extract strand information (same for all samples in study)
        strands = str(data.strands.values[0])
        
        # Create batch dictionary with samples key
        batch_dict = {'strands': strands, 'samples': {}}
        
        for batch, data_batch in data.groupby('batch'):
            batch_dict['samples'][str(batch)] = data_batch.run_accession.values.tolist()
        
        # Write configuration file
        with open(output[0], "w") as file:
            yaml.dump(batch_dict, file, default_flow_style=False)


rule dada_filter_and_trim_paired:
    """
    Filter and trim paired-end reads using DADA2.
    
    DADA2 filterAndTrim:
        - Removes low-quality bases
        - Trims adapters and primers
        - Filters reads by quality score
        - Removes reads below minimum length
        - Outputs filtered FASTQ files and statistics
        
    Note:
        Uses batch-specific parameters via symlinked parameter files
    """
    input:
        R1 = config["study_dir"] + "{study_name}/{sample}/{sample}_1.fastq.gz",
        R2 = config["study_dir"] + "{study_name}/{sample}/{sample}_2.fastq.gz",
        dada2_params = config["params_dir"] + "filter_and_trim/{study_name}.yaml"
    output:
        log_file = config["study_dir"] + "processed/trimmed/filter_and_trim/{study_name}/{sample}/{sample}_2.log"
    params:
        R1_out = config["study_dir"] + "processed/trimmed/filter_and_trim/{study_name}/{sample}/{sample}_1.fastq.gz",
        R2_out = config["study_dir"] + "processed/trimmed/filter_and_trim/{study_name}/{sample}/{sample}_2.fastq.gz"
    conda:
        f"{conda_dir}/dada2.yaml"
    shell:
        f"""
        Rscript {scripts_dir}/download_trim/dada2_filter_trim.r \
            {{input.R1}} {{input.R2}} {{input.dada2_params}} \
            {{params.R1_out}} {{params.R2_out}} {{output.log_file}}
        """


rule create_sample_table_trimmed:
    """
    Create summary table of trimming statistics for all samples.
    
    This rule:
        - Collects all trimming log files
        - Parses statistics (reads in, reads out, reads filtered)
        - Generates CSV table for quality assessment
        - Used for tracking data loss through filtering pipeline
    """
    input:
        trimming_logs = get_trimming_logs
    output:
        config["study_dir"] + "processed/trimmed/{study_name}.csv"
    params:
        reads_dir = config["study_dir"] + "processed/trimmed/filter_and_trim/{study_name}/"
    conda:
        f"{conda_dir}/dada2.yaml"
    shell:
        f"""
        Rscript {scripts_dir}/download_trim/create_sample_table.r \
            {{params.reads_dir}} {{output}}
        """
