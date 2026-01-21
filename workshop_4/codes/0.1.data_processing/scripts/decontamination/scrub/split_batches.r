library(phyloseq)
args <- commandArgs(TRUE)

# LOAD PARAMS
metadata <- read.csv(args[[1]], row.names=1)
otu_data <- read.csv(args[[2]], row.names=1)
col_diagnosis <- args[[3]]
col_well <- args[[4]]
col_batch <- args[[5]]
metadata_batches_path <- args[[6]]
otu_data_batches_path <- args[[7]]

metadata_scrub <- metadata[, c(col_diagnosis, col_batch, col_well)]
metadata_scrub <- metadata[!is.na(metadata[col_batch]), ]
metadata_scrub['is_control'] <- unlist(lapply(metadata_scrub[, col_diagnosis], function(x){
                                    if (is.na(x)){
                                        return(FALSE)
                                    }
                                    if (x == 'Negative Control') {
                                        return(TRUE)
                                    } else {
                                        return(FALSE)
                                    }}))
batches <- unique(metadata_scrub[[col_batch]])
metadata_scrub <- metadata_scrub[, c('is_control', col_diagnosis, col_well, col_batch)]
colnames(metadata_scrub) <- c('is_control', 'sample_type', 'sample_well', col_batch)
metadata_scrub <- lapply(batches, function(batch){
    return(metadata_scrub[metadata_scrub[, col_batch] == batch,
                          c('is_control', 'sample_type', 'sample_well')])
})
names(metadata_scrub) <- batches

batches_names <- names(metadata_scrub) 
otu_data_scrub <- lapply(batches_names, function(batch){
    return(otu_data[row.names(metadata_scrub[[batch]]),])
})
names(otu_data_scrub) <- batches

saveRDS(metadata_scrub, metadata_batches_path)  
saveRDS(otu_data_scrub, otu_data_batches_path)  