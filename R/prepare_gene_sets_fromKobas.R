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
keggOrgId <- list('S' = c('bamy'),
                  'P' = c('pae'))
############################

for (bac in bacteriaList){
  kobasAnnot <- readKobasAnnotation(file = kobasFiles[[bac]])
  
  genesAndPathways <- keggLink(keggOrgId[[bac]], 'pathway')
  refGeneSets_long <- data.table('refGeneID' = genesAndPathways, 'kegg.ko' = names(genesAndPathways))
  geneSets_long <- merge.data.table(x = kobasAnnot,
                                    y = refGeneSets_long,
                                    by.x = 'keggID',
                                    by.y = 'refGeneID')
  
  # convert gene set long format to list format
  uniqueKos <- unique(geneSets_long$kegg.ko)
  geneSets <- sapply(uniqueKos, function(ko){geneSets_long[kegg.ko == ko, query]})
  
  saveRDS(geneSets, file = paste0('data/tidy/geneSets/geneSets_Kobas_', bac, '.rds'))
}

