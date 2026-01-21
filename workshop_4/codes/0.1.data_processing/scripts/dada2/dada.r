##################################################################
# Step dada: run the core dada2 algorithm of ASV identification.
#     Utilizes a error model to identify ASVs on dereplicated data.
#     Runs forward and reverse strands separately.
#     Creates one object per strand for all samples in the dataset.
# Output: a dada-class object in .rds format
##################################################################

# load input data
library(yaml) 
library("dada2")
args <- commandArgs(TRUE)

# read and save the input
err <- readRDS(args[[1]])    # read the error model
derep <- readRDS(args[[2]])  # read the file with dereplicated reads
out <- c(args[[3]]) # a path for the output 
params <- yaml.load_file(args[[4]]) # read a pool parameter for the dada pipeline

# Unlsee the pool parameter is pseudo, transform to the bool format
if (params$pool != 'pseudo') {
    params$pool <- as.logical(params$pool)
}
params$derep <- derep
params$err <- err
params$multithread <- TRUE
params$verbose <- TRUE

# load libraries
# run a dada algorithm
pool <- do.call(dada, params)
# save the result to the rds format
saveRDS(pool, out)