## Main workflow: PSconsortia
## Constributors: viktorht@yahoo.dk

# Install required packeges
if (!requireNamespace("BiocManager", quietly = TRUE)) # Bioconductor installer
  install.packages("BiocManager")

if (!requireNamespace("DESeq2", quietly = TRUE)){ # DESeq2 package
  BiocManager::install("DESeq2")
}

if (!requireNamespace("KEGGREST", quietly = TRUE)){ # KEGGREST package
  BiocManager::install("KEGGREST")
}

if (!requireNamespace("ggplot2", quietly = TRUE)){ # ggplot2 package
  install.packages('ggplot2')
}


if (!requireNamespace("purrr", quietly = TRUE)){ # purr package
  install.packages('purrr')
}

if (!requireNamespace("data.table", quietly = TRUE)){ # DESeq2 package
  install.packages('data.table')
}

# Script sets up directories
library(purrr)
library(here)

folder_names <- c("data/raw", 
                  "data/tidy", "data/tidy/deseq2", "data/tidy/geneSets", 'data/tidy/geneSetEnrichmentAnalysis',
                  "refs", 
                  "R", 
                  "analysis", 
                  "figures/deseq2", 
                  "man")
map(folder_names, dir.create) # Creates folder 