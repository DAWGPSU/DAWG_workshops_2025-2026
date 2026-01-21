##################################################################
# Step convert_rds: create a phyloseq object and incorporates metadata.
#    Creates a phyloseq object that contains:
#    read counts (otu_table), taxonomy (tax_table), and provided metadata (sample_data).
# Output: a phyloseq object save in a .rds format
##################################################################

# load input data
args <- commandArgs(TRUE)

# read and save the input
seqtab <- readRDS(args[[1]]) # reads sequence table
taxa <- readRDS(args[[2]])   # reads a taxonomy table
# metadata <-  read.table(args[[3]], sep='\t', header=1) # reads a metadata table
output_path <- args[[4]] # an output file path

metadata <- read.csv(args[[3]], row.names = 1, sep='\t')
metadata <- metadata[row.names(seqtab), ]

# load libraries
library("phyloseq")
# combine all the tables into a phyloseq object
ps <- phyloseq(otu_table(seqtab, taxa_are_rows = FALSE),
               tax_table(taxa), 
               sample_data(metadata))
# save in rds format
saveRDS(ps, output_path)   