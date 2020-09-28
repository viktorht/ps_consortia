## Main workflow: PSconsortia
## Constributors: viktorht@yahoo.dk

## Plotting richfactor col plot 
rm(list = ls())
library(data.table)
library(here)
library(ggplot2)
bacteriaList <- c('P', 'S')
richfactorCut <- 0.3
padjCut <- 0.1
sizeCut <- 5

# load and aggregate data
dfAll <- data.frame() # df for data collection of all strains
for (bac in bacteriaList){
  pathwayNames <- read.table(paste0("data/tidy/geneSets/pathwayNames_", bac, ".txt")) # read kegg pathway names
  richfactor <- read.csv(file = paste0("data/tidy/geneSetEnrichmentAnalysis/richFactor_", bac, ".csv")) # read richfactor data
  hyperTestResults <- read.csv(file = paste0("data/tidy/geneSetEnrichmentAnalysis/hyperTestResults_", bac, ".csv")) # read results of hypergeometric test
  
  df <- merge.data.frame(x = pathwayNames, # merge data to one df
                         y = richfactor)
  df <- merge.data.frame(x = df,
                         y = hyperTestResults,
                         by = 'kegg.ko')
  dfAll <- rbind.data.frame(dfAll, df) # combine with other strains
  rm(df) # cleans
}

# format data for plotting
dfAll <- as.data.table(dfAll)
dfAll.long <- melt.data.table(data = dfAll[, .(Bacteria, kegg.ko, pathway.name, size, richfactor.up, richfactor.down, padj.up, padj.down)],
                                     measure.vars = patterns('richfactor', 'padj'), # use all columns that contain the two patterns as measurement variables
                                     id.vars = c('Bacteria', 'kegg.ko', 'pathway.name', 'size'), 
                                     value.name = c('richfactor', 'padj'), # naming the columns
                                     variable.name = 'direction') # naming the columns
dfAll.long[direction == 1, direction := 'up'] # renaming value of assigned variable during melting of dfAll
dfAll.long[direction == 2, direction := 'down']


mask <- dfAll.long[, kegg.ko] %in% dfAll.long[ ((abs(richfactor) > richfactorCut) & 
                                                  (padj < padjCut) & 
                                                  (size > sizeCut))][, kegg.ko] 

# Plotting
# col plot
ggplot(dfAll.long[mask, with = TRUE]) +
  geom_col(aes(x = pathway.name, y = richfactor, fill = padj)) +
  coord_flip() +
  geom_hline(yintercept=0) +
  facet_wrap(.~Bacteria)+
  scale_colour_gradientn(values = c(0, 0.05, 0.1, 1), 
                         breaks = c(0, 0.05, 0.5, 0.1, 1), 
                         colors = c('red', 'yellow', 'blue'), 
                         aesthetics = "fill", 
                         name = 'FDR')+
  xlab('')+
  ylab('Richfactor')
ggsave(filename = 'figures/richfactorColPlot.png')

