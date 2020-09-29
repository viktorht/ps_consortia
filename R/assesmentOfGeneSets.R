## Main workflow: PSconsortia
## Constributors: viktorht@yahoo.dk
## Assesment of gene set methods

library(here)
library(data.table)
library(KEGGREST)
source("R/custom_functions.R")

bac <- 'P'

# load files
geneSets_kobas <- readRDS(file = paste0('data/tidy/geneSets/geneSets_Kobas_', bac, '.rds'))
geneSets_eggnog <- readRDS(file = paste0('data/tidy/geneSets/geneSets_', bac, '.rds'))
richfactor_kobas <- read.csv(file = paste0('data/tidy/geneSetEnrichmentAnalysis/richFactor_kobas_', bac, '.csv'))
richfactor_eggnog <- read.csv(file = paste0('data/tidy/geneSetEnrichmentAnalysis/richFactor_eggnog_', bac, '.csv'))

# number of gene set assocciations 
sum(richfactor_eggnog$size)
sum(richfactor_kobas$size)
# There are way fewer genes associated to gene sets in kobas (for P)


summary(richfactor_eggnog$size)
summary(richfactor_kobas$size)

boxplot(data.frame('eggnog' = richfactor_eggnog$size, 'kobas' = richfactor_kobas$size), outline = F)


# length of BCAA synthesis gene set 
length(geneSets_long[kegg.ko == 'path:pae00290', query]) # from kobas
length(geneSets_eggnog$ko00290) # from eggnog

geneSets_eggnog$ko00290[!(geneSets_eggnog$ko00290 %in% geneSets_long[kegg.ko == 'path:pae00290', query])]
# Three genes are not in the kobas gene set