library(phyloseq)
args <- commandArgs(TRUE)

# LOAD PARAMS
ps <- readRDS(args[[1]])
ps_NC_path <- args[[2]]
sample_data_NC_path <- args[[3]]
otu_table_NC_path <- args[[4]]

ps <- subset_samples(ps, (sample != 'p72_S72_L001')&(sample != ''))

ps_NC <- subset_samples(ps, (i.Diagnosis == "Negative Control"))
saveRDS(ps_NC, ps_NC_path)

write.csv(otu_table(ps), otu_table_NC_path)
write.csv(data.frame(sample_data(ps)), sample_data_NC_path)