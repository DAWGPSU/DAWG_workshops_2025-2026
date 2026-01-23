#!/usr/bin/env python3
"""
Convert Taxonomy to FASTA Format for Phylogenetic Tree Building

This script processes taxonomy and OTU tables from phyloseq objects:
1. Renames OTU IDs to simplified format (OTU1, OTU2, etc.)
2. Creates FASTA file from taxonomy sequences
3. Formats taxonomy for ghost-tree compatibility

Usage:
    python convert_taxonomy_to_fasta.py \
        input_otu.csv \
        input_taxonomy.csv \
        output_otu_renamed.csv \
        output_taxonomy_renamed.csv \
        output_taxonomy_ghost.csv \
        output_fasta.fasta

Author: Data Analysis Workshop
"""

import pandas as pd
import sys

# Parse command line arguments
if len(sys.argv) != 7:
    print("Error: Incorrect number of arguments")
    print("Usage: python convert_taxonomy_to_fasta.py input_otu input_tax output_otu output_tax output_ghost output_fasta")
    sys.exit(1)

input_otu_file = sys.argv[1]
input_tax_file = sys.argv[2]
output_otu_file = sys.argv[3]
output_tax_file = sys.argv[4]
output_tax_ghost_file = sys.argv[5]
output_fasta_file = sys.argv[6]

print("Loading input files...")
print(f"  OTU file: {input_otu_file}")
print(f"  Taxonomy file: {input_tax_file}")

# Load input files
df_tax = pd.read_csv(input_tax_file, sep=",", index_col=0)
df_otu = pd.read_csv(input_otu_file, sep=",", index_col=0)

# Fill missing taxonomy values with empty strings
df_tax = df_tax.fillna("")

# Extract sequences (these are the original OTU IDs/sequences)
sequences = df_tax.index.values

# Create simplified OTU names (OTU0, OTU1, OTU2, ...)
otu_names = ['OTU%d' % i for i, seq in enumerate(sequences)]

print(f"\nProcessing {len(sequences)} OTUs...")

# Create FASTA file with renamed OTU IDs
print(f"Creating FASTA file: {output_fasta_file}")
otu_rename = {}
with open(output_fasta_file, 'w') as f:
    for name, seq in zip(otu_names, sequences):
        otu_rename[seq] = name
        f.write(">%s\n%s\n" % (name, seq))

# Update taxonomy table with new OTU names
df_tax.index = otu_names
df_tax.to_csv(output_tax_file, sep="\t", index=True, header=True)
print(f"Created renamed taxonomy file: {output_tax_file}")

# Update OTU table with new OTU names
df_otu.columns = [otu_rename[val] for val in df_otu.columns.values]
df_otu.to_csv(output_otu_file, sep="\t", index=True, header=True)
print(f"Created renamed OTU file: {output_otu_file}")

# Format taxonomy for ghost-tree
# Ghost-tree expects specific prefixes for each taxonomic level
ghosttree_levels = ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
taxonomy_cols = df_tax.columns[0:7]  # First 7 columns are Domain to Species

def format_tax(row):
    """
    Format taxonomy string for ghost-tree compatibility.
    
    Args:
        row: DataFrame row containing taxonomy information
        
    Returns:
        Semicolon-separated taxonomy string with ghost-tree prefixes
    """
    cleaned = []
    for i in range(7):
        val = row[taxonomy_cols[i]]
        prefix = ghosttree_levels[i]
        
        if val != '':
            # Clean up the taxonomy value
            val = val.strip()
            
            # For species level, extract only the species epithet
            if i == 6:  # Species level
                parts = val.split()
                val = parts[-1] if parts else ""
            
            # Replace spaces with underscores
            val = val.replace(" ", "_")
            cleaned.append(val)
        else:
            # Use just the prefix for missing levels
            cleaned.append(prefix)
    
    # Join all levels with semicolons
    return ";".join(cleaned)

# Apply formatting function to all rows
print("Formatting taxonomy for ghost-tree...")
df_tax_ghost = df_tax.apply(format_tax, axis=1)

# Save ghost-tree formatted taxonomy
df_tax_ghost.to_csv(output_tax_ghost_file, sep="\t", index=True, header=False)
print(f"Created ghost-tree taxonomy file: {output_tax_ghost_file}")

print("\nConversion complete!")
print(f"  Generated {len(otu_names)} renamed OTUs")
print(f"  Output files:")
print(f"    - {output_otu_file}")
print(f"    - {output_tax_file}")
print(f"    - {output_tax_ghost_file}")
print(f"    - {output_fasta_file}")
