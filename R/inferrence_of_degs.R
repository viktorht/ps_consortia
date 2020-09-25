##########
# Viktor Hesselberg-Thomsen 
##########

# DESeq2 analysis of the count data

rm(list=ls())
library(here)
library(data.table)
library(DESeq2)

# First I create a function that prepares the counts file and runs the DESeq2 analysis
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
    png(paste0('figures/deseq2/MAplot_', direction,'_', outsuffix, '.png'))
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
  write.csv(as.data.frame(resLFC), 
            file=paste0("data/tidy/deseq2/","DESeq2Results_",direction,  '-', bac, '_', outsuffix,"_shrunkenLFC.csv"))
  
}

 

# Run DESeq2 analysis 
bacteria <- c("P", "S")
for (bac in bacteria){
  cnts.file <- paste0('data/raw/',bac, '/4.GeneExprQuatification/4.1.GeneExprQuatification/readcount.xls')
  use.conditions <- 'all'
  direction <- paste0("PSvs", bac)
  outsuffix = paste0('conditions_', paste0(use.conditions, collapse = '_'))
  
  do.DESeq2(cnts.file, use.conditions, direction, outsuffix, do.MAplots = T)
  summary(res)
  summary(resLFC)
}

  