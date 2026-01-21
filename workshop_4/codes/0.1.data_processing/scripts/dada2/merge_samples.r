library(yaml) 

args <- commandArgs(TRUE)
# LOAD PARAMS
params <- yaml.load_file(args[[1]])
samples = params$samples
merged_files = params$merged_files
output_file <- args[[2]]
 
merged_list <- list()
for(i in seq_along(merged_files)) {
     file <- merged_files[i]
     merged_list[[i]] <- readRDS(file)
}
if (length(merged_list) == 1) {
    merged_list <- merged_list[[1]]
  } else {
    if (is.null(names(merged_files))) {
      names(merged_list) <- samples
    } else {
      names(merged_list) <- samples
    }
  }
saveRDS(merged_list, output_file)