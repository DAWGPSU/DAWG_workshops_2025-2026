library(phyloseq)

args <- commandArgs(TRUE)

# LOAD PARAMS
ps <- readRDS(args[[1]])
output_ps <- args[[2]]

ps <- subset_taxa(ps, (Kingdom=="Bacteria")|(Kingdom=="Archaea"))
ps
ps <- subset_taxa(ps, (Family!="Mitochondria"))
ps
ps <- subset_samples(ps, (sample != "p72_S72_L001"))
ps
saveRDS(ps, output_ps)
