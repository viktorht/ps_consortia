## Main workflow: PSconsortia
## Constributors: viktorht@yahoo.dk
# Create gene sets from eggnog mapper
rm(list=ls()) # Removes variables 
library(here)
library(data.table)
library(KEGGREST)
source("R/custom_functions.R")

###### Input data ##########
bacteriaList <- c(quote(P), quote(S)) # quote is used to the make use of the bac variable to subset data.table
updatePathwayNameFile <- FALSE
kobasFiles <- list("S" = c('data/raw/kobasAnnotation/sqr9.txt'), 
                    "P" = c("data/raw/kobasAnnotation/pae.txt"))
keggOrgId <- list('S' = c('bamy'), # the org id that was used by kobas as reference 
                  'P' = c('pae'))
############################

for (bac in bacteriaList){
  kobasAnnot <- readKobasAnnotation(file = kobasFiles[[bac]]) # load the kobas file
  
  genesAndPathways <- keggLink(keggOrgId[[bac]], 'pathway') # get the mapping between kegg Gene ID and pathways
  keggGeneSets_long <- data.table('keggGeneID' = genesAndPathways, 'kegg.ko' = names(genesAndPathways)) # creates a data table of the above list
  geneSets_long <- merge.data.table(x = kobasAnnot, # merge query IDs with keggID to map query IDs to genes sets
                                    y = keggGeneSets_long,
                                    by.x = 'keggID',
                                    by.y = 'keggGeneID')
  
  # convert gene set long format to list format
  uniqueKos <- unique(geneSets_long$kegg.ko)
  geneSets <- sapply(uniqueKos, function(ko){geneSets_long[kegg.ko == ko, query]})
  
  saveRDS(geneSets, file = paste0('data/tidy/geneSets/geneSets_Kobas_', bac, '.rds'))
}

