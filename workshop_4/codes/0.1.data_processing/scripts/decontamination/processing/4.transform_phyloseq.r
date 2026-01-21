##################################################################
# Step transform_rds: Agglomerate and normalize phyloseq objects
#     Agglomerates ASV at the genus and family levels.
#     Normalizes counts to relative abundances if requested.
#     The resulting phyloseq objects are saved in the .rds format.
#     In addition, for compatibility with other programming languages, 
#     each part of the phyloseq object data (otu_table, tax_table, sample_data) 
#     is saved in the form of csv tables.
# Output: a phyloseq object in .rds format
#         an ASV table in .csv format
#         a taxonomy table in .csv format
#         a metadata table in .csv format
##################################################################

# load input data
args <- commandArgs(TRUE)

# read and save the input
ps <- readRDS(args[[1]])          # read the phyloseq object
output_path_ps <- args[[2]]       # an output path for the phyloseq object (rds)
output_path_otu <- args[[3]]      # an output path for the ASV table (csv)
output_path_tax <- args[[4]]      # an output path for the taxonomy table (csv)
output_path_meta <- args[[5]]     # an output path for the sample data table (csv)
agglom <- args[[6]]               # a parameter that indicates the level of ASV agglomeration
normalization <- args[[7]]        # a parameter for normalization type

# load libraries
library("phyloseq")

# agglomerate ASVs at the genus level, if requested
if (agglom == 'genus') {
    ps <- tax_glom(ps, "Genus", NArm = FALSE)
}

# agglomerate ASVs at the family level, if requested
if (agglom == 'family') {
    ps <- tax_glom(ps, "Family", NArm = FALSE)
}

# normalize to relative abundance if requested
if (normalization == 'normalized') {
    ps <- transform_sample_counts(ps, function(x) x / sum(x))
}

# save a phyloseq object to an rds format
saveRDS(ps, output_path_ps)

# save ASV table to a csv format
write.csv(otu_table(ps), output_path_otu)

# save taxonomy table to a csv format
write.csv(tax_table(ps), output_path_tax)

# save metadata table to a csv format
write.csv(data.frame(sample_data(ps)), output_path_meta)
