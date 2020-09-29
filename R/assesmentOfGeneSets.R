## Main workflow: PSconsortia
## Constributors: viktorht@yahoo.dk
## Assesment of gene set methods

rm(list=ls()) # Removes variables 
library(here)
library(data.table)
library(KEGGREST)
source("R/custom_functions.R")

bac <- 'P'

# load files
geneSets_kobas <- readRDS(file = paste0('data/tidy/geneSets/geneSets_Kobas_', bac, '.rds'))
geneSets_eggnog <- readRDS(file = paste0('data/tidy/geneSets/geneSets_', bac, '.rds'))



# length of BCAA synthesis gene set 
length(geneSets_long[kegg.ko == 'path:pae00290', query]) # from kobas
length(geneSets_eggnog$ko00290) # from eggnog

geneSets_eggnog$ko00290[!(geneSets_eggnog$ko00290 %in% geneSets_long[kegg.ko == 'path:pae00290', query])]
# Three genes are not in the kobas gene set