library(phyloseq)
args <- commandArgs(TRUE)

# LOAD PARAMS
otu_decontaminated <- read.table(args[[1]], row.names=1)
ps <- readRDS(args[[2]])
ps_decontaminated_path <- args[[3]]

ps_decontaminated <- phyloseq(otu_table(otu_decontaminated, taxa_are_rows=FALSE),
                              tax_table(ps),
                              sample_data(ps))

saveRDS(ps_decontaminated, ps_decontaminated_path)