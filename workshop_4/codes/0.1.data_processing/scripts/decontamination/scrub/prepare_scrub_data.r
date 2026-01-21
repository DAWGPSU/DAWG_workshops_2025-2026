library(phyloseq)
args <- commandArgs(TRUE)

# LOAD PARAMS
ps_decontaminated_path <- readRDS(args[[1]])
ps_NC <- readRDS(args[[2]])

otu_table_path <- args[[3]]
sample_data_path <- args[[4]]

# Merge the two phyloseq objects
ps_merged <- merge_phyloseq(ps_decontaminated_path, ps_NC)
write.csv(otu_table(ps_merged), otu_table_path)
write.csv(data.frame(sample_data(ps_merged)), sample_data_path)
