## Main workflow: PSconsortia
## Constributors: viktorht@yahoo.dk

library(here)
library(data.table)

if (use.data == "PS"){
  eggnog.files <- list("S" = c("data/raw/eGGNOGmapper/Bacillus.tsv"), 
                       "P" = c("data/raw/eGGNOGmapper/job_MM_skhrncty_annotations_pseudomonas.tsv"))
  outFolderName <- 'DESeq2based'
  outsuffix <- paste0('DESeq2_PS_all_shrunkenLFC')
  mainwd <- getwd()
  bacteriaList <- c(quote(S), quote(P)) # quote is used to the make use of the bac variable to subset data.table
}

degCountsAll <- data.table(NULL) # for collecting data
for (bac in bacteriaList) {
  # load deseq2 results
  dat.file <- list.files('data/clean/DESeq2/', pattern = paste0(use.data, 'vs', bac, '.*all.*shrunken'))
  deseq2Results <- fread(input = paste0('data/clean/DESeq2/', dat.file))
  setnames(deseq2Results, 'V1', 'geneID')
  
  # Make geneSets from eggnog
  eggnogTable <- readEggnogTable(file = eggnog.files[[bac]])
  geneSets <- makeKEGGGeneSets(eggnogTable = eggnogTable)
  
  # Calculate rich factor 
  degCounts <- calcRichfactor(degseq2Results = deseq2Results, geneSets = geneSets, padjCutOff = 0.05, lfcCutOff = 1) 
  degCounts[, Bacteria := as.character(bac)] # adds the bacteria 
  
  degCounts <- merge(degCounts, GeneSetHyperTest(degResults = deseq2Results, 
                                                 geneSet = geneSets, 
                                                 lfcCut = 1, 
                                                 pvalCut = 0.05))
  
  degCountsAll <- rbind(degCountsAll, degCounts) # collects data
}
