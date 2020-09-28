## Main workflow: PSconsortia
## Constributors: viktorht@yahoo.dk
# Functions for Gene set enrichment analysis

readEggnogTable <- function(file) {
  #### Loading the eggnog annotation file. 
  # Header in dataframe is incomplete there for the header is ommitted and added afterwards
  library(data.table)
  eggnogData <- read.delim(file = file, 
                           skip = 4, header = F, sep = "\t", stringsAsFactors = F)
  
  tableHeader <- c("Query.name","Seed.ortholog","e.value","Score","best.tax.lvl","Preferred.name","GO.terms","EC.number","KEGG.KO","KEGG.Pathway","KEGG.Module","KEGG.Reaction","KEGG.rclass","BRITE","KEGG.TC","CAZy","BiGG.Reaction","annot.lvl.max","eggNOG.OGs","BestOG","COG.cat","Description")
  # also found in file "data/raw/eGGNOGmapper/eggnogResultsHeader.txt"
  
  colnames(eggnogData) <- tableHeader # Adding header to the df
  
  #### Removing last 3 rows that contains metadata
  n <- dim(eggnogData)[1]
  
  print("Removed metadata: ")
  print(paste0(eggnogData[((n-2):n),][,1]))
  
  eggnogData <- eggnogData[(1:(n-3)),]
  
  return(as.data.table(eggnogData))
}


makeCogGeneSets <- function(eggnogTable) {
  
  cogSymbols <- c("A", "B", "C", "D", "E", "F", "G", "H","I", "J", "K","L","M", "N", "O", "P", "Q", "T", "U", "Y", "Z", "R", "S")
  geneSets <- list() # empty list
  
  for (cog in cogSymbols) {
    geneMask <- grepl(pattern = cog, eggnogTable[,"COG.cat"]) # Creates a mask for all occurances of a cog
    genesNames <- as.character(eggnogTable[geneMask ,"Query.name"]) # Get gene names
    geneSets[[cog]] <- genesNames # Writes to list
    
  }
  return(geneSets)
}

makeKEGGGeneSets <- function(eggnogTable) {
  '
  Takes an table output from eggnog-mapper and makes gene sets based on the predictions from eggnog.
  At the moment it is hardcoded to take KEGG.Pathway column, but could be rewritten to give more flexibility.
  '
  
  allTerms <- unlist(strsplit(eggnogTable[,KEGG.Pathway], split = ",")) # Each row in the kegg.pathway column is a string listing the pathway seperated by , this is made into one longe list
  uniqueTerms <- unique(allTerms) # As the same pathway appears for more genes we find the unique set of pathways present
  koPaths <- uniqueTerms[grepl("^ko", uniqueTerms)] # kegg.pathway contain both ko and map, which are redundant, map is removed here
  
  
  geneSets <- list() # empty list
  for (koPathway in koPaths) { # obs data.table syntax is used
    geneMask <- grepl(pattern = koPathway, eggnogTable[,KEGG.Pathway]) # Creates a mask for all occurances of a pathway
    genesNames <- as.character(eggnogTable[geneMask ,Query.name]) # Get gene names
    geneSets[[koPathway]] <- genesNames # Writes to list
    
  }
  return(geneSets)
}


makeRankedGeneList <- function(fpkm, method = "signal2noise") {
  # remove unused data, i.e. dicultures
  dicultures <- sapply(c(1:6), function(num, bac){paste0(bac, num)}, bac = c("BS", "AB", "PS", "AS"))
  fpkm.dat <- fpkm.dat[ ,!(names(fpkm.dat) %in% dicultures)]
  print("OBS dicultures data where removed!")
  
  # Setup masks to fetch only mono or triculture data
  triCultureMask <- grepl("^ABS", colnames(fpkm.dat)) # Makes a bool to select the columns with triculture data
  monoCulturePatterns <- "(^S\\d$)|(^A\\d$)|(^B\\d$)" # regexp pattern to catch any of the three monoculture 
  monoCultureMask <- grepl(pattern = monoCulturePatterns, x = colnames(fpkm.dat)) 
  
  # Calculate the mean and sd for mono and tricultures. uses mask to determine columns.
  fpkm.dat["mean.tri"] <- apply(fpkm.dat[, triCultureMask], MARGIN = 1, FUN = mean)
  fpkm.dat["sd.tri"] <- apply(fpkm.dat[, triCultureMask], MARGIN = 1, FUN = sd)
  fpkm.dat["mean.mono"] <- apply(fpkm.dat[, monoCultureMask], MARGIN = 1, FUN = mean)
  fpkm.dat["sd.mono"] <- apply(fpkm.dat[, monoCultureMask], MARGIN = 1, FUN = sd)
  
  if (method == "signal2noise"){
    # Calc signal to noise ratio 
    # Gives positive value if upregulated and negative downregulated
    fpkm.dat["sig2noi"] <- (fpkm.dat["mean.tri"] - fpkm.dat["mean.mono"]) / (fpkm.dat["sd.tri"] + fpkm.dat["sd.mono"])
    
    # Sorting gene list
    sig2noi <- fpkm.dat[["sig2noi"]] # Creating obj of class numeric with sig2noi values
    names(sig2noi) <- fpkm.dat[["geneID"]] # naming numeric with geneID 
    rankedGenes <- sort(sig2noi, decreasing = T) # Sorting for ranked list 
  } else {print(paste0(method, " is not a valid method."))}
  
  return(rankedGenes)
}

addCogDescription <- function(fgseaResults) {
  '
  Adds COG descriptions to fgsea output data table.
  '
  cogDescriptions <-  as.data.table(read.delim("data/raw/eGGNOGmapper/COG_descriptions.txt", header = T, skip = 3)) # reads in the descriptions and convert to data table to be compatible with output of fgsea-package
  names(cogDescriptions) <- c("pathway", "description") # Change colnames to match output of fgsea for easy merging 
  fgseaRes <- merge(cogDescriptions, fgseaRes, by = "pathway") # add the description column
  
}

writeGmtFile <- function(geneSets, filepath){
  gmtRows <- c()
  for (i in 1:length(geneSets)) {
    gmtRows <- c(gmtRows, paste(c(names(geneSets)[i], geneSets[[i]]), collapse ="\t"))
  }
  confile <- file(filepath)
  writeLines(gmtRows, con = confile, sep = "\n")
  close(confile)
  
}

fetchKeggNames <- function(koNumbers, file = "") {
  '
  Takes a character vector of koNumbers of kegg pathways and return the a dataframe with the name for these pathways.
  
  Currently takes only pathway IDs, but could be expanded to take others IDs
  '
  library(KEGGREST) # load package
  library(data.table)
  
  uniqueKoNumbers <- unique(koNumbers) # As some pathways are enriched in more than one strain
  keggEntry <- list()
  i <- 0
  while (i < length(uniqueKoNumbers)){ # you can only submit 10 queries at once (limitation from kegg server)
    i <- i + 10
    keggEntry <- append(keggEntry, keggGet(uniqueKoNumbers[(i-9):i]))
  }
  
  keggNames <- sapply(keggEntry, function(x){c(x$ENTRY, x$NAME)}) # getting koNumber and name of the pathway
  keggNames <- t(as.data.table(keggNames)) # making a dataframe, actually a matrix
  
  colnames(keggNames) <- c("kegg.ko", "pathway.name") # naming columns 
  
  if (length(keggNames[,'pathway.name']) != length(uniqueKoNumbers)){ # give warning 
    message("Number of input pathways does not match number of pathway in output. This can be due to pathways that were not found in KEGG")
  }
  
  if (length(file) > 1){ # if file path is given the data is saved to a file
    write.table(keggNames, file = file, row.names = F)
    print(paste("Written to file ", file)) # notify user
  }
  return(keggNames)
}

fpkm2tpm <- function(fpkm){
  '
  Takes fragments per kilobase transcript per milloin mapped reads (FPKM) and converts it into transcripts per million (TPM).
  This is doen simply by 
  '
  fpkmSampleSum <- lapply(fpkm[,-1], sum) # removes geneID column and calculates the sum of each column
  tpm <- (fpkm[,-1] / fpkmSampleSum) * 1e6 # Calculates transcripts per million
  tpm <- cbind(fpkm[,1], tpm) # Adds geneIDs to tpm files
  colnames(tpm)[1] <- colnames(fpkm)[1]
  return(tpm)
}

readKovasAnnotation <- function(file){
  '
  Reads output from kovas annotation online tool in .txt format. 
  returns data.table format
  '
  library(data.table)
  kovasAnnotation <- as.data.table(read.csv(file, skip = 3, sep = "|",  header = T, stringsAsFactors = F)) # Read in the data and convert to data.table fread does not work directly on data
  lastDataRow <- kovasAnnotation[X.Query.Gene.ID == "--------------------", which = T] - 1 # A lot of extra data is at the last part of the file. Row number is identified
  kovasAnnotation <- kovasAnnotation[1:lastDataRow] # removing the not needed data
  kovasAnnotation[, c("query", "keggID"):= tstrsplit(X.Query.Gene.ID, "\t")] # Splitting query id and kegg id 
  kovasAnnotation[, X.Query.Gene.ID := NULL] # Remove now splitted column
  setcolorder(kovasAnnotation, c("query", "keggID")) # reorder the columns
  return(kovasAnnotation)
}

convertGeneSet2keggId <- function(geneSet, kovasAnnotation, removeNoneAnnotations = TRUE){
  
  # Convert gene set from geneID to keggGeneID
  geneSets.kegg <- lapply(geneSet, function(set){kovasAnnotation[query %in% set][,keggID]})
  
  # Test if the structure of the two list are identical
  if (!class(geneSet) == class(geneSets.kegg)){ # test class
    return("Error: class' does not match.")
  } 
  if (0 != sum(!(unlist(lapply(geneSet, length)) == unlist(lapply(geneSets.kegg, length))))) { # Tests length a all gene sets
    return("Error: Gene set lengths does not match")
  } 
  
  if (removeNoneAnnotations){
    # Remove none values
    # Kovas does not find annotation for all genes. Some of these none annotated gene maybe present in the gene sets. These are removed.
    geneSets.kegg <- lapply(geneSets.kegg, function(set){Filter(function(set){!(set %in% "None")}, set)})
  }
  
  return(geneSets.kegg)
}

filterGeneSet <- function(geneSet, upperBound, lowerBound){
  '
  Takes a gene set of filter removes the sets with lengths outside the bounderies deffined 
  '
  length(Filter(function(x){length(x) < lowerBound}, geneSet)) # Number of gene sets with less than 3 genes
  length(Filter(function(x){length(x) > upperBound}, geneSet)) # Number of gene sets with more than 200 genes
  geneSet.filtered <- Filter(function(x){(length(x) > lowerBound & length(x) < upperBound)}, geneSet)
  return(geneSet.filtered)
}

gageOfLFC <- function(bac, degResults.file, eggnogTable, outsuffix, outFolderName, makePathviewFiles = FALSE, pathviewOnlyDEGCutoff = FALSE){
  geneSet <- makeKEGGGeneSets(eggnogTable = eggnogTable) # Make gene set from eggNOG mapper.
  
  # Loading DEG results from company results files
  degResults <- fread(degResults.file)
  setnames(degResults, "V1", "Gene_id")
  
  lfc <-  degResults[, .(Gene_id, log2FoldChange)]
  
  
  lfc.df <- data.frame("log2FoldChange" = lfc[,log2FoldChange], row.names = lfc[, Gene_id]) # Formats for gage
  
  # Running gage
  gageResults <- gage(lfc.df, gsets = geneSet, set.size = c(3,500)) # set.size creates limits for gene sets
  
  # Test that pairing makes no difference
  #gageResults.unpaired <- gage(lfc.df, gsets = geneSet.filtered, compare = "unpaired")
  #gageResults.paired <- gage(lfc.df, gsets = geneSet.filtered, compare = "paired")
  #gageResults.paired$less == gageResults.unpaired$less
  
  sigResults <- sigGeneSet(gageResults, cutoff = 0.1, heatmap = FALSE) # get significant results 
  pathIds.greater <- rownames(sigResults$greater)
  pathIds.less <- rownames(sigResults$less)
  
  # Identifying essential gene sets
  # q.val is used as this is also the cirteria used in later analysis
  setwd(paste0("C:/Users/Viktor/OneDrive - Danmarks Tekniske Universitet/SpecialeF20/data/clean/geneEnrichment/", outFolderName, "/"))
  if (length(pathIds.greater) > 1){ # esset.grp fails if not
    essential.greater <- esset.grp(setp = gageResults$greater, exprs = lfc.df, gsets = geneSet, outname = paste0(bac, '_essensGreat', outsuffix), use.q = T, cutoff = 0.1)
  }
  if (length(pathIds.less) > 1) { # esset.grp fails if not
    essential.less <- esset.grp(setp = gageResults$less, exprs = lfc.df, gsets = geneSet, outname = paste0(bac, "_essensLess", outsuffix), use.q = T, cutoff = 0.1)
    
  }
  
  
  # Preparing expression data with ko id for pathview and joyplot.dt
  #eggnogTable[(nchar(KEGG.KO) == 0 & nchar(KEGG.Pathway) > 0)] # Checks the we don't missing pathways by using ko 
  lfc.ko <- merge.data.table(lfc, eggnogTable[, .(Query.name, KEGG.KO)], by.x = "Gene_id", by.y = "Query.name") # adds ko terms to lfc expression data
  dt.split <- lfc.ko[, tstrsplit(KEGG.KO, ",")] # Some gene have more ko ids these are seperated by ,
  lfc.ko.split <- cbind(lfc.ko, dt.split) # columns containing the extra ko ids are appended
  lfc.ko.long <- melt.data.table(lfc.ko.split, id.vars = c("Gene_id", "log2FoldChange", "KEGG.KO"), na.rm = T) # a new row is made for all ko ids such that some genes appear multiple time in the list
  lfc.ko.long <- cbind(lfc.ko.long, lfc.ko.long[, tstrsplit(value, ":")[2]]) # the format of the ko ids are ko:XXXXX, but only XXXXX is requried for pathview
  lfc.ko.long[, c("KEGG.KO", "variable", "value"):=NULL] # Left over columns are removed
  setnames(lfc.ko.long, "V1", "kegg.ko") # remnaing column
  
  # Set non DEGs to fold change of 0
  if (is.numeric(pathviewOnlyDEGCutoff)){ # Only used when a value is given for pathviewOnlyDEGsCutoff
    degResults[, "significant" := padj < pathviewOnlyDEGCutoff] # Find significant DEGs
    lfc.ko.long.sig <- merge.data.table(lfc.ko.long, degResults[, .(Gene_id, significant)]) # Merges boolean to long formated data
    lfc.ko.long.sig[significant == FALSE, log2FoldChange := 0] # Set non significant DEGs lfc = 0
    
    exp.ko <- lfc.ko.long.sig[,log2FoldChange] # log2 fold change expression data
    names(exp.ko) <- lfc.ko.long.sig[, kegg.ko] # naming the vector to format correctly 
    
    outsuffix <- paste0(outsuffix, "_DEGsOnly") # Adds to outsuffix to be able to differentiate files
  }else{
    exp.ko <- lfc.ko.long[,log2FoldChange] # log2 fold change expression data
    names(exp.ko) <- lfc.ko.long[, kegg.ko] # naming the vector to format correctly 
    
  }
  
  if (makePathviewFiles){
    
    library(pathview)
    setwd(paste0("C:/Users/Viktor/OneDrive - Danmarks Tekniske Universitet/SpecialeF20/results/GeneEnrichment/", outFolderName, "/pathviewGAGE/"))
    
    dontPlot <- c("ko01110", "ko01120", "ko01230")
    pathIds.greater <- pathIds.greater[!(pathIds.greater %in% dontPlot)]
    pathIds.less <- pathIds.less[!(pathIds.less %in% dontPlot)]
    print("Omit plotting list:")
    print(paste(dontPlot, collapse = ","))
    
    
    
    for (id in pathIds.greater){
      pathview(
        gene.data = exp.ko, pathway.id = id, # Input data and pathway id
        species = "ko", out.suffix=paste0(bac, "_",outsuffix, "_greater"), kegg.native = T, # species ko is used as K-ids are used instead of gene names
        limit = list(gene = 5), bins = list(gene = 20), # setting scale of legend
        low = list(gene = "red"), mid = list(gene = "grey"), high = list(gene = "green")) # Seting colour of legends
      Sys.sleep(1)
    }
    
    for (id in pathIds.less){
      pathview(
        gene.data = exp.ko, pathway.id = id, # Input data and pathway id
        species = "ko", out.suffix=paste0(bac, "_", outsuffix, "_less"), kegg.native = T, # species ko is used as K-ids are used instead of gene names
        limit = list(gene = 5), bins = list(gene = 20), # setting scale of legend
        low = list(gene = "red"), mid = list(gene = "grey"), high = list(gene = "green")) # Seting colour of legends
      Sys.sleep(1)
    }
  }
  
  # Save data for joyplot
  geneSetLong <- geneSet2longformat(geneSet, bacteria = bac)
  joyplot.dt <- merge(lfc.ko.long, geneSetLong)
  
  return(list(gageResults = gageResults, joyplot.dt = joyplot.dt)) 
}


geneSet2longformat <- function(geneSet, bacteria){
  '
  Creates a long format data table with columns gene_id, pathway and bacteria.
  Used to make gage results in into correct format for joyplot.
  '
  geneSetNames <- names(geneSet)
  geneSetLong <- as.data.table(c(NULL, NULL))
  for (name in geneSetNames){
    dt.tmp <- as.data.table(geneSet[[name]])
    dt.tmp <- cbind(dt.tmp, name, as.character(bacteria))
    geneSetLong <- rbind(geneSetLong, dt.tmp)
  }
  names(geneSetLong) <- c("Gene_id", "pathway", "bacteria")
  return(geneSetLong)
}

createAllResultsDT <- function(bacterialist, allResults, pathwayNames, outpath){
  '
  Assembles gage results into one data table. Also adds pathwaynames to the table.
  Used in gage analysis
  '
  library(data.table)
  allResultsDT <- as.data.table(NULL)
  for (bac in bacteriaList){
    for (dir in c(quote(greater), quote(less))){
      allResultsDT <- rbind(allResultsDT, as.data.table(allResults[[bac]][[dir]], keep.rownames = TRUE)[, "bacteria" := as.character(bac)][,"direction" := as.character(dir)])
    }
  }
  allResultsDT <- na.omit(allResultsDT) 
  setnames(allResultsDT, "rn", "kegg.pathway")
  allResultsDT <- merge.data.table(x = allResultsDT, y = pathwayNames, by.x = "kegg.pathway", by.y = "kegg.pathway")
  setcolorder(allResultsDT, c("bacteria", "kegg.pathway", "pathway.name", "direction"))
  fwrite(allResultsDT, file = paste0(outpath, "gageAllResuls_", outsuffix, ".txt")) # Write file for use in other scripts
  print(paste0("allResultsDT saved to ", outpath, "gageAllResuls_", outsuffix, ".txt"))
  return(allResultsDT)
}

plotPathviewMapKO <- function(pathIds, deseq2Results, eggnogTable, pvOutsuffix, bac, outpath, pathviewOnlyDEGCutoff = FALSE){
  '
  Wrapper function for the pathview-function in pathview package.
  Plot a pathway map and colouring the enzymes according to their log2foldchange.
  If pathviewOnlyDEGCutoff is set to a number this is interpreted as the significance level to apply to the DEGs.
  The log2foldchange of non-significant genes is set to 0
  '
  lfc <-  deseq2Results[, .(Gene_id, log2FoldChange)]
  
  # Preparing expression data with ko id for pathview and joyplot.dt
  #eggnogTable[(nchar(KEGG.KO) == 0 & nchar(KEGG.Pathway) > 0)] # Checks the we don't missing pathways by using ko 
  lfc.ko <- merge.data.table(lfc, eggnogTable[, .(Query.name, KEGG.KO)], by.x = "Gene_id", by.y = "Query.name") # adds ko terms to lfc expression data
  dt.split <- lfc.ko[, tstrsplit(KEGG.KO, ",")] # Some gene have more ko ids these are seperated by ,
  lfc.ko.split <- cbind(lfc.ko, dt.split) # columns containing the extra ko ids are appended
  lfc.ko.long <- melt.data.table(lfc.ko.split, id.vars = c("Gene_id", "log2FoldChange", "KEGG.KO"), na.rm = T) # a new row is made for all ko ids such that some genes appear multiple time in the list
  lfc.ko.long <- cbind(lfc.ko.long, lfc.ko.long[, tstrsplit(value, ":")[2]]) # the format of the ko ids are ko:XXXXX, but only XXXXX is requried for pathview
  lfc.ko.long[, c("KEGG.KO", "variable", "value"):=NULL] # Left over columns are removed
  setnames(lfc.ko.long, "V1", "kegg.ko") # remnaing column
  
  # Set non DEGs to fold change of 0
  if (is.numeric(pathviewOnlyDEGCutoff)){ # Only used when a value is given for pathviewOnlyDEGsCutoff
    deseq2Results[, "significant" := padj < pathviewOnlyDEGCutoff] # Find significant DEGs
    lfc.ko.long.sig <- merge.data.table(lfc.ko.long, deseq2Results[, .(Gene_id, significant)]) # Merges boolean to long formated data
    lfc.ko.long.sig[significant == FALSE, log2FoldChange := 0] # Set non significant DEGs lfc = 0
    
    exp.ko <- lfc.ko.long.sig[,log2FoldChange] # log2 fold change expression data
    names(exp.ko) <- lfc.ko.long.sig[, kegg.ko] # naming the vector to format correctly 
    
    outsuffix <- paste0(outsuffix, "_DEGsOnly") # Adds to outsuffix to be able to differentiate files
  }else{
    exp.ko <- lfc.ko.long[,log2FoldChange] # log2 fold change expression data
    names(exp.ko) <- lfc.ko.long[, kegg.ko] # naming the vector to format correctly 
    
  }
  setwd(outpath)
  for (id in pathIds){
    pathview(
      gene.data = exp.ko, pathway.id = id, # Input data and pathway id
      species = "ko", out.suffix=paste0(bac, "_",pvOutsuffix), kegg.native = T, # species ko is used as K-ids are used instead of gene names
      limit = list(gene = 5), bins = list(gene = 20), # setting scale of legend
      low = list(gene = "red"), mid = list(gene = "grey"), high = list(gene = "yellow")) # Seting colour of legends
    Sys.sleep(1)
  }
}

calcRichfactor <- function(degseq2Results, geneSets, padjCutOff, lfcCutOff){
  geneSet.lengths <- lapply(geneSets, length) # Count number of DEGs in a gene set
  
  sigDEGs <- deseq2Results[padj < padjCutOff, geneID]
  sigDEGsUp <- deseq2Results[(padj < padjCutOff & log2FoldChange > lfcCutOff), geneID]
  sigDEGsDown <- deseq2Results[(padj < padjCutOff & log2FoldChange < -1*lfcCutOff), geneID]
  
  
  # Count up, down and total DEGs in each gene set
  countUp <- lapply(geneSets, function (geneSet) {sum(sigDEGsUp %in% geneSet)})
  countDown <- lapply(geneSets, function (geneSet) {sum(sigDEGsDown %in% geneSet)})
  countAll <- lapply(geneSets, function (geneSet) {sum(sigDEGs %in% geneSet)})
  
  
  degCounts <- data.table('kegg.ko' = names(countUp),'degUp' = unlist(countUp), 'degDown' = unlist(countDown), 'deg' = unlist(countAll), 'size' = unlist(geneSet.lengths))
  
  # calculate rich factor 
  degCounts$richfactor.up <- with(degCounts, degUp / size)
  degCounts$richfactor.down <- with(degCounts, -1 * (degDown / size))
  
  return(degCounts)
}

GeneSetHyperTest <- function(degResults, geneSet, lfcCut, pvalCut, padj.method = 'BH'){
  # Takes DESeq2results dt x and a list of gene sets and cut off values and return the probability 
  # x number of enriched genes from a give gene set, given the size of the gene set and the total number of 
  # enriched genes.
  # Inspired by https://stats.stackexchange.com/questions/16247/calculating-the-probability-of-gene-list-overlap-between-an-rna-seq-and-a-chip-c
  
  #### Principle in hypergeometric test
  # For gene set X calculate:
  # Prob of drawing [enriched and part of gene set X] or more white balls when drawing [# genes in X] times in a population 
  # with [# enriched genes] white balls and [# non-enriched genes] black balls.
  # Checking for down regulation
  df <- data.frame('kegg.ko' = names(geneSet), 'pval.down' = NA)
  
  enriched <- degResults[(log2FoldChange < -lfcCut & padj < pvalCut)]
  notEnriched <- degResults[!(log2FoldChange < -lfcCut & padj < pvalCut)]
  
  numEnriched <- dim(enriched)[1] # number of enriched genes (white balls)
  numNotEnriched <- dim(notEnriched)[1] # number of not enriched genes (Black balls)
  
  for (i in c(1:length(geneSet))){
    lenGeneSet <- length(geneSet[[i]]) # Number of times to draw
    numEnrichedInGeneSet <- sum(geneSet[[i]] %in% enriched$geneID) # Number of success/white balls 
    
    df$pval.down[i] <- phyper(numEnrichedInGeneSet-1, numEnriched, numNotEnriched, lenGeneSet, lower.tail = FALSE)
  }
  
  # Checking for up regulation
  df$pval.up <- NA
  enriched <- degResults[(log2FoldChange > lfcCut & padj < pvalCut)]
  notEnriched <- degResults[!(log2FoldChange > lfcCut & padj < pvalCut)]
  
  numEnriched <- dim(enriched)[1] # number of enriched genes (white balls)
  numNotEnriched <- dim(notEnriched)[1] # number of not enriched genes (Black balls)
  
  for (i in c(1:length(geneSet))){
    lenGeneSet <- length(geneSet[[i]]) # Number of times to draw
    numEnrichedInGeneSet <- sum(geneSet[[i]] %in% enriched$geneID) # Number of success/white balls 
    
    df$pval.up[i] <- phyper(numEnrichedInGeneSet-1, numEnriched, numNotEnriched, lenGeneSet, lower.tail = FALSE)
  }
  
  # FDR correction
  df$padj.down <- p.adjust(df$pval.down, padj.method)
  df$padj.up <- p.adjust(df$pval.up, padj.method)
  
  return(df)
}

indexDuplicates <- function(vector){
  forward <- duplicated(vector, fromLast = F)
  backward <- duplicated(vector, fromLast = T)
  return(forward | backward)
}

lfcOfGenesInPathway <- function(bac, pathway.id, joyplot.dt2, result.files){
  deseq2Results <- fread(paste0('data/clean/DESeq2/', result.files[[bac]]))
  setnames(deseq2Results, "V1", "Gene_id")
  
  genes <- joyplot.dt2[(bacteria == bac & pathway == pathway.id & direction == 'greater'), 
                       .(bacteria, pathway, Gene_id, kegg.ko, pathway.name)]
  
  sel <- deseq2Results$Gene_id %in% genes$Gene_id
  res <- merge.data.table(genes, deseq2Results[sel, .(Gene_id, log2FoldChange, padj)], by = 'Gene_id')
  res <- res[order(kegg.ko,decreasing=TRUE),]
  return(res)
}
