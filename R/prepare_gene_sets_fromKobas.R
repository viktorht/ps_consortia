## Main workflow: PSconsortia
## Constributors: viktorht@yahoo.dk
# Create gene sets from eggnog mapper
rm(list=ls()) # Removes variables 
library(here)
library(data.table)
library(KEGGREST)
source("R/custom_functions.R")

