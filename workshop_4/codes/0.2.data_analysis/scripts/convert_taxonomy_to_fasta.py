import pandas as pd
import sys

input_otu_file = sys.argv[1]
input_tax_file = sys.argv[2]
output_otu_file = sys.argv[3]
output_tax_file = sys.argv[4]
output_tax_ghost_file = sys.argv[5]
output_fasta_file = sys.argv[6]

# Load input
df_tax = pd.read_csv(input_tax_file, sep=",", index_col=0)
df_otu = pd.read_csv(input_otu_file, sep=",", index_col=0)
df_tax = df_tax.fillna("")

sequences = df_tax.index.values
otu_names = ['OTU%d'%i for i, seq in enumerate(df_tax.index.values)]

otu_rename = {}
with open(output_fasta_file, 'w') as f:
    for name, seq in zip(otu_names, sequences):
        otu_rename[seq] = name
        f.write(">%s\n%s\n"%(name, seq))

df_tax.index = otu_names
df_tax.to_csv(output_tax_file, sep="\t", index=True, header=True)

df_otu.columns = [otu_rename[val] for val in df_otu.columns.values]
df_otu.to_csv(output_otu_file, sep="\t", index=True, header=True)

# Column order expected for ghost-tree
ghosttree_levels = ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
taxonomy_cols = df_tax.columns[0:7]  # first 7 columns are Domain to Species

# Remove _1 suffix and build taxonomy string
def format_tax(row):
    cleaned = []
    for i in range(7):
        val = row[taxonomy_cols[i]]
        prefix = ghosttree_levels[i]
        if val != '':
            val = val.strip()#.replace("_1", "")
            if i == 6:  # species level
                parts = val.split()
                val = parts[-1] if parts else ""
            val = val.replace(" ", "_")
            cleaned.append(val)
        else:
            cleaned.append(prefix)
    return ";".join(cleaned)
df_tax = df_tax.apply(format_tax, axis=1)

# Output file
df_tax.to_csv(output_tax_ghost_file, sep="\t", index=True, header=False)

