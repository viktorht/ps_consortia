## Main workflow: PSconsortia
## Constributors: viktorht@yahoo.dk

# Script sets up directories
library(purrr)
library(here)

folder_names <- c("data/raw", 
                  "data/tidy", "data/tidy/deseq2", "data/tidy/geneSets",
                  "refs", 
                  "R", 
                  "analysis", 
                  "figures/deseq2", 
                  "man")
map(folder_names, dir.create) # Creates folder 