library(phyloseq)
library(decontam)

args <- commandArgs(TRUE)

# LOAD PARAMS
ps <- readRDS(args[[1]])
output_ps <- args[[2]]
output_decontam <- args[[3]]

contamdf.freq <- isContaminant(ps, method="frequency", conc='concentrations', normalize=FALSE)
table(contamdf.freq$contaminant)
contamdf.freq.sort <- contamdf.freq[order(contamdf.freq$p),]
contamdf.freq.sort['name'] = unlist(lapply(row.names(contamdf.freq.sort), function(x){
                                    paste(tax_table(ps)[x, ], collapse=';')}))
taxa <- row.names(contamdf.freq.sort[contamdf.freq.sort$contaminant == FALSE, ])
ps_freq <- prune_taxa(taxa, ps)

saveRDS(ps_freq, output_ps)
write.csv(contamdf.freq.sort, output_decontam)

