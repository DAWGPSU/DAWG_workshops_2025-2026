library(phyloseq)

args <- commandArgs(TRUE)

# LOAD PARAMS
ps <- readRDS(args[[1]])
output_ps <- args[[2]]
threshold <- args[[3]]

ps_normalized <- ps
ps_normalized <- transform_sample_counts(ps_normalized, function(x){x / sum(x)})
ps_normalized_samples_only <- subset_samples(ps_normalized, (i.Diagnosis != "Negative Control"))
ps_normalized_nc_only <- subset_samples(ps_normalized, (i.Diagnosis == "Negative Control"))

# remove low represented based on 0.5.1 threshold=0.01
maxval <- as.vector(apply(otu_table(ps_normalized_samples_only), 2, max)) # max abundance for each bacteria
names(maxval) <- colnames(otu_table(ps_normalized_samples_only))
ps <- subset_taxa(ps, maxval > threshold)
ps

saveRDS(ps, output_ps)
