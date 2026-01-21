library(phyloseq)
library(decontam)

args <- commandArgs(TRUE)

# LOAD PARAMS
ps <- readRDS(args[[1]])
output_ps <- args[[2]]
output_decontam <- args[[3]]
threshold <- as.numeric(args[[4]])

sample_data(ps)$is.neg <- sample_data(ps)$type == "Negative Control"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=threshold, normalize=FALSE)

table(contamdf.prev$contaminant)
contamdf.prev <- contamdf.prev[order(contamdf.prev$p),]
contamdf.prev['name'] = unlist(lapply(row.names(contamdf.prev), function(x){
                                    paste(tax_table(ps)[x, ], collapse=';')}))
taxa <- row.names(contamdf.prev[contamdf.prev$contaminant == FALSE, ])
ps_freq <- prune_taxa(taxa, ps)

saveRDS(ps_freq, output_ps)
write.csv(contamdf.prev, output_decontam)

