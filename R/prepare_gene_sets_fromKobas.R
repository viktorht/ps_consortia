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
bac <- 'P'

kobasAnnot <- readKobasAnnotation(file = kobasFiles[[bac]])

genesAndPathways <- keggLink('pae', 'pathway')
refGeneSets_long <- data.table('refGeneID' = genesAndPathways, 'kegg.ko' = names(genesAndPathways))
geneSets_long <- merge.data.table(x = kobasAnnot,
                                 y = refGeneSets_long,
                                 by.x = 'keggID',
                                 by.y = 'refGeneID')

# compare the gene sets from kobas to eggnog
geneSets_eggnog <- readRDS(file = 'data/tidy/geneSets/geneSets_P.rds')

# length of BCAA synthesis gene set 
length(geneSets_long[kegg.ko == 'path:pae00290', query]) # from kobas
length(geneSets_eggnog$ko00290) # from eggnog

geneSets_eggnog$ko00290[!(geneSets_eggnog$ko00290 %in% geneSets_long[kegg.ko == 'path:pae00290', query])]
# Three genes are not in the kobas gene set