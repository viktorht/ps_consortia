## Main workflow: PSconsortia
## Constributors: viktorht@yahoo.dk

## Main file to run all
library(here)

source('R/setup.R')

# Data processing
source("R/inferrence_of_degs.R") # use DESeq2 to identify DEGs
source('R/prepare_gene_sets.R') # create kegg based gene sets from eggnog annotation

# Analysis
source('R/geneSetEnrichmentAnalysis_hyperTest.R')
source('R/richFactor.R')

# plotting
