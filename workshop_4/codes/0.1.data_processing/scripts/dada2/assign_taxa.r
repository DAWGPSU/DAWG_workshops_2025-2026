##################################################################
# Step assign_taxa: assign taxonomy
#    Assigns taxonomy using SILVA and SILVA species databases.
# Output: a taxonomy table in a .rds format with columns: Kingdom, Phylum, Class, Order, Family, Genus, Species.
##################################################################

# load input data
args <- commandArgs(TRUE)

# read and save the input
seqtab <- readRDS(args[[1]]) # reads sequence table
db_loc <- c(args[[2]])       # a path to the SILVA database
if (length(args) == 4){
    db2_loc <- c(args[[3]])      # a path to the SILVA species database
    out <- c(args[[4]])          # a path to the output file
} else {
    out <- c(args[[3]])          # a path to the output file
}

# load libraries
library("dada2")
# assign taxonomy
taxa <- assignTaxonomy(seqtab, db_loc, multithread=TRUE, tryRC=TRUE)
if (length(args) == 4){
    # assign species
    taxa <- addSpecies(taxa, db2_loc)
}
# save in rds format
saveRDS(taxa, out)