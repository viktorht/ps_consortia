##########
# Viktor Hesselberg-Thomsen 24-09-2020
##########

# Script sets up directories
library(purrr)
library(here)

folder_names <- c("data/raw", "data/tidy", "refs", "R", "analysis", "figures", "man")
#sapply(folder_names, dir.create) # base R way
map(folder_names, dir.create) #purrr-fect way