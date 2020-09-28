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
############################
bac <- 'S'

kobasAnnot <- readKobasAnnotation(file = kobasFiles[[bac]])

genesAndPathways <- keggLink('bamy', 'pathway')
refGeneSets_long <- data.table('refGeneID' = genesAndPathway, 'kegg.ko' = names(genesAndPathways))
geneSets_long <- merge.data.table(x = kobasAnnot,
                                 y = refGeneSets_long,
                                 by.x = 'keggID',
                                 by.y = 'refGeneID')

# compare the gene sets from kobas to eggnog
geneSets_eggnog <- readRDS(file = 'data/tidy/geneSets/geneSets_S.rds')

# length of BCAA synthesis gene set 
length(geneSets_long[kegg.ko == 'path:bamy00290', query]) # from kobas
length(geneSets_eggnog$ko00290) # from eggnog
# The length is the same

sum(!geneSets_long[kegg.ko == 'path:bamy00290', query] %in% geneSets_eggnog$ko00290)
# The two method finds the exact same genes to be part of that pathway
