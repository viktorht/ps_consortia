library(here)
library(data.table)
source('R/custom_functions.R')

use.data = 'PS'
use.geneSets = 'eggnog'
if (use.data == "PS"){
  eggnog.files <- list("S" = c("data/raw/eGGNOGmapper/Bacillus.tsv"), 
                       "P" = c("data/raw/eGGNOGmapper/job_MM_skhrncty_annotations_pseudomonas.tsv"))
  dat.file <- list('P' = c('DESeq2Results_PSvsP-P_conditions_all_shrunkenLFC.csv'),
                   'S' = c('DESeq2Results_PSvsS-S_conditions_all_shrunkenLFC.csv'))
  bacteriaList <- c(quote(S), quote(P)) # quote is used to the make use of the bac variable to subset data.table
}


for (bac in bacteriaList) {
  geneSets <- loadGeneSets(folderPath = 'data/tidy/geneSets/', # selected the correct gene set file
                           bacteria = bac, 
                           use.geneSets = use.geneSets)
  
  # load deseq2 results
  deseq2Results <- fread(input = paste0('data/tidy/DESeq2/', dat.file[[bac]]))
  setnames(deseq2Results, 'V1', 'geneID')
  
  # Calculate rich factor 
  degCounts <- calcRichfactor(degseq2Results = deseq2Results, geneSets = geneSets, padjCutOff = 0.05, lfcCutOff = 2) 
  degCounts[, Bacteria := as.character(bac)] # adds the bacteria 

  
  write.csv(x = degCounts, # Saves data to file
            file = paste0('data/tidy/geneSetEnrichmentAnalysis/richFactor_', use.geneSets, '_', bac, '.csv'),
            row.names = F)
}
