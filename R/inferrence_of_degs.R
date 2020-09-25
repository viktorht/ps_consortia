##########
# Viktor Hesselberg-Thomsen 
##########

# DESeq2 analysis of the count data

rm(list=ls())
setwd("C:/Users/Viktor/OneDrive - Danmarks Tekniske Universitet/SpecialeF20")
library(data.table)
library(DESeq2)

do.DESeq2 <- function(cnts.file, use.conditions, direction, outsuffix = "", do.MAplots = FALSE){
  # Read the count data
  cnts <- as.matrix(read.delim(cnts.file, 
                               sep = '\t', row.names = "geneID"))
  sort(colnames(cnts))
  # Choose samples
  if (length(use.conditions) == 1){
    if (use.conditions == "all"){
      cnts <- cnts[, sort(colnames(cnts))]
    } else {stop("Only one condition is given. Give atleast 2 condition or 'all'")}
  } else if (class(use.conditions) == "character") {
    samples.choice <- paste0(rep(use.conditions, each = 6), seq(1,6))
    cnts <- cnts[, sort(samples.choice)] # Often throughs error if samples are not found 
  } else {
    stop("use.conditions is invalid. Need to be of class character.")
  }
  head(cnts)
  # pre-filtering not done
  # keep <- rowSums(cnts[,-1], na.rm = T) != 0 # samples with zero in all 
  # cnts <- cnts[keep, ]
  
  # Setting experiment design
  if (length(use.conditions) == 1){
    use.conditions <- unique(gsub(pattern = "[0-9]", replacement = '', colnames(cnts)))
  }
  coldata <- data.frame("condition" = sort(rep(use.conditions, each = 6)))
  rownames(coldata) <- colnames(cnts)
  
  # test match between coldata and cols in data matrix 
  col.conditions <- gsub(pattern = '[0-9]', replacement = '', colnames(cnts))
  if (!all(col.conditions == coldata$condition)) {
    print("Conditions in coldata does not match columns of data matrix.")
    print(coldata)
    stop("Execution stopped.")
  }
  
  # Create DESeq2 object 
  dds <- DESeqDataSetFromMatrix(countData = cnts,
                                colData = coldata,
                                design = ~ condition)
  
  dds$condition <- relevel(dds$condition, ref = unlist(strsplit(direction, "vs"))[2]) # changes control condition last strain given in direction 
  
  
  # Calc. DEGs
  dds <- DESeq(dds)
  res <- results(dds, name = paste0("condition_", gsub(pattern = "vs", "_vs_", direction))) 
  
  # lfcShrinkage
  resLFC <- lfcShrink(dds, coef = paste0("condition_", gsub(pattern = "vs", "_vs_", direction)), type="apeglm")
  
  # plotting 
  if (do.MAplots) {
    png(paste0('results/DESeq2/MAplot_', direction,'_', outsuffix, '.png'))
    par(mfrow=c(1,2))
    plotMA(res, ylim=c(-5,5), main = "Normal")
    plotMA(resLFC, ylim = c(-5,5), main = "Shrunken")
    dev.off()
  }
  
  # Assigning global variables so results can be investigated futher outside function
  dds <<- dds
  res <<- res
  resLFC <<- resLFC
  
  # Exporting 
  write.csv(as.data.frame(res), 
            file=paste0("data/clean/DESeq2/","DESeq2Results_",direction,  '-', bac, '_', outsuffix, ".csv"))
  write.csv(as.data.frame(resLFC), 
            file=paste0("data/clean/DESeq2/","DESeq2Results_",direction,  '-', bac, '_', outsuffix,"_shrunkenLFC.csv"))
  
}

<<<<<<< HEAD
bacteria <- c("A", "S")
for (bac in bacteria){
  cnts.file <- paste0('data/raw/trancriptome data/Transcriptome data ABS/',bac, '/4.GeneExprQuatification/4.1.GeneExprQuatification/readcount.xls')
  use.conditions <- 'all' #c("ABS", bac)
  direction <- paste0("ABSvs", 'AS')
  =======
    bacteria <- c("A", "B", "S")
  bacteria <- 'S'
  for (bac in bacteria){
    cnts.file <- paste0('data/raw/trancriptome data/Transcriptome data ABS/',bac, '/4.GeneExprQuatification/4.1.GeneExprQuatification/readcount.xls')
    use.conditions <- 'all'
    direction <- paste0("ABSvs", bac)
    >>>>>>> growthStateEstMethod1
    outsuffix = paste0('conditions_', paste0(use.conditions, collapse = '.'))
    
    do.DESeq2(cnts.file, use.conditions, direction, outsuffix, do.MAplots = T)
    summary(res)
    summary(resLFC)
  }
  
  <<<<<<< HEAD
  =======
    saveRDS(dds, file = 'data/clean/DESeq2/RdataFiles/S_dds.rds')
  >>>>>>> growthStateEstMethod1
  
  # For PS sample
  bacteria <- c("P", "S")
  for (bac in bacteria){
    cnts.file <- paste0('data/raw/trancriptome data/Transcriptome data ABS/',bac, '/4.GeneExprQuatification/4.1.GeneExprQuatification/readcount.xls')
    use.conditions <- c("PS", bac)
    direction <- paste0("PSvs", bac)
    outsuffix = paste0('conditions_', paste0(use.conditions, collapse = '.'))
    
    do.DESeq2(cnts.file, use.conditions, direction, outsuffix, do.MAplots = T)
    summary(res)
    summary(resLFC)
  }
  
  