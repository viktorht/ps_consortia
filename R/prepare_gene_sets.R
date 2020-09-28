# Create gene sets from eggnog mapper
rm(list=ls()) # Removes variables 
library(here)
library(data.table)
library(KEGGREST)
source("R/custom_functions.R")

bac <- 'S'
updatePathwayNameFile <- TRUE
eggnog.files <- list("S" = c("data/raw/eGGNOG_mapper/Bacillus.tsv"), 
                     "P" = c("data/raw/eGGNOG_mapper/job_MM_skhrncty_annotations_pseudomonas.tsv"))
eggnogTable <- readEggnogTable(eggnog.files[[1]])
geneSets <- makeKEGGGeneSets(eggnogTable = eggnogTable)


if (updatePathwayNameFile) { # This takes a few minutes 
  geneSetNames <- fetchKeggNames(names(geneSets)) # fetch the names of the pathways ways 
  write.table(geneSetNames, paste0('data/raw/pathway_names_', bac, '.txt')) # save to file 
}
geneSetNames <- read.table(paste0('data/raw/pathway_names_', bac, '.txt'))

geneSets[!(names(geneSets) %in% geneSetNames$kegg.ko)] # find pathways that were not found in kegg



# Create genes sets
# Add kegg names

bacteriaList <- c(quote(P), quote(S)) # quote is used to the make use of the bac variable to subset data.table