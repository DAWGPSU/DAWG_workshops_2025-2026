library(phyloseq)

args <- commandArgs(TRUE)

# LOAD PARAMS
ps_decontam <- readRDS(args[[1]])
ps_scrub <- readRDS(args[[2]])
output_ps <- args[[3]]
                                  
ps <- prune_taxa(taxa_names(ps_decontam), ps_scrub)

saveRDS(ps, output_ps)

