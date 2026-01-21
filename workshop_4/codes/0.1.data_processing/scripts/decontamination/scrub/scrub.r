if (!requireNamespace("SCRuB", quietly = TRUE)) {
    devtools::install_github("shenhav-and-korem-labs/SCRuB", upgrade=FALSE)
    torch::install_torch()
}

library(SCRuB)
args <- commandArgs(TRUE)

metadata_scrub <- readRDS(args[[1]])
otu_data_scrub <- readRDS(args[[2]])
batches_names <- names(otu_data_scrub)
decontaminated_data_path <- args[[3]]

decontaminated_data_ <- lapply(batches_names, function(batch){
                        print(batch)
                        res <- SCRuB(otu_data_scrub[[batch]], metadata_scrub[[batch]])
                        return(res$decontaminated_samples)
})
decontaminated_data <- decontaminated_data_[[1]]
res <- lapply(c(2:length(decontaminated_data_)), function(i){
    decontaminated_data <<- rbind(decontaminated_data, decontaminated_data_[[i]])
})
write.table(decontaminated_data, decontaminated_data_path)
