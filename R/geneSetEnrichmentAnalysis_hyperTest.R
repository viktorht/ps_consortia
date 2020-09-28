## Main workflow: PSconsortia
## Constributors: viktorht@yahoo.dk

library(here)
library(data.table)
source('R/custom_functions.R')

use.data = 'PS'

if (use.data == "PS"){
  eggnog.files <- list("S" = c("data/raw/eGGNOGmapper/Bacillus.tsv"), 
                       "P" = c("data/raw/eGGNOGmapper/job_MM_skhrncty_annotations_pseudomonas.tsv"))
  dat.file <- list('P' = c('DESeq2Results_PSvsP-P_conditions_all_shrunkenLFC.csv'),
                   'S' = c('DESeq2Results_PSvsS-S_conditions_all_shrunkenLFC.csv'))
  bacteriaList <- c(quote(S), quote(P)) # quote is used to the make use of the bac variable to subset data.table
}

degCountsAll <- data.table(NULL) # for collecting data
for (bac in bacteriaList) {
  # load deseq2 results
  deseq2Results <- fread(input = paste0('data/tidy/DESeq2/', dat.file[[bac]]))
  setnames(deseq2Results, 'V1', 'geneID')
  
  # Make geneSets from eggnog
  geneSets <- readRDS(paste0('data/tidy/geneSets/geneSets_', bac, '.rds'))
  
  hyperTestResults <- GeneSetHyperTest(degResults = deseq2Results, 
                                       geneSet = geneSets, 
                                       lfcCut = 2, 
                                       pvalCut = 0.05)
  
  write.csv(x = hyperTestResults, # Saves data to file
            file = paste0('data/tidy/geneSetEnrichmentAnalysis/hyperTestResults_', bac, '.csv'),
            row.names = F)
}

