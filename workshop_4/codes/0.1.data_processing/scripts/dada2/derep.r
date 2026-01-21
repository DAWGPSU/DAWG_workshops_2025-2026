##################################################################
# Step derep: dereplicate reads.
#    Dereplicates the reads for each sample in the dataset.
#    Creates one object per strand for all samples in the dataset.
# Output: a derep-class object in a .rds format.
##################################################################

# load libraries
library("dada2")

# load input data
args <- commandArgs(TRUE)
input_file <- args[[1]] # a fastq.gz file path or table file path
output_path <- args[[2]] # an output file path
strand <- args[[3]] # read the strand the algorithm will be applied to

# read and save the input
if (grepl(".fastq.gz", input_file, fixed = TRUE)) {
    files_in <- c(input_file)
    samples <- c(strand)
} else {
    reads <- read.table(file = args[[1]], sep = '\t', header = TRUE) # read a sample table
    files_in <- as.character(reads[[strand]]) # ensure non-factors format
    samples <- as.character(reads$sample) # ensure non-factors format
}
# dereplicate
derep <- derepFastq(files_in, verbose=TRUE)
# ensure the correct names
names(derep) <- samples
# save to the rds format
saveRDS(derep, output_path)