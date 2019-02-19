## Install/load libraries
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
## Shorten name for Tx library
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
##
#source("https://bioconductor.org/biocLite.R")
#biocLite("genomation", suppressUpdates = TRUE)
# used for readBed()
library(genomation)

## ** NEED TO USE ENSEMBL TRANSCRIPT ID, NOT GENE ID ** ## ???? no

## The input file (2490 lines)
inputFile <- "C:\\Users\\jsk33\\Documents\\atac\\mydata\\snu61\\snu61_tfc"
## Read in the file
expData <- read.table(inputFile, header=TRUE, sep = ',')

## Show the available data (columns) from TxDb
columns(txdb)
colum <- c("TXID", "TXNAME")
## Show available database select keytypes from TxDb
keytypes(txdb)

## Make a column key to retrieve all info for each gene
#colKey <- c(columns(txdb)) # for all columns
## Make a gene ID based key for retrieving data for each gene of interest
geneKey <- c(as.character(expData[,3]))

## Use the select method to get mapping between tx_id and gene_id
annotData <- select(txdb, keys = geneKey, columns = columns(txdb), keytype = "GENEID")
## Make a vector with the tx_names
txNames <- c(annotData[,2])

## Get promoters for all transcripts from TxDb
## Upstream and downstream can be custom set, default seems to be -2000/+200 from TSS
promoters <- promoters(txdb)
## Trim the GRanges object to keep standard entries only
promoters <- keepStandardChromosomes(promoters, pruning.mode="coarse")
## Subset the promoters GRanges object using the generated index
## Note that with multiple transcript variants, this number will be much higher than the gene_id list
promoterData <- promoters[promoters$tx_name %in% txNames]

##
newTxName <- promoterData@elementMetadata@listData[["tx_name"]]
##
newGeneID <- c()
##
for (a in 1:length(newTxName)){
  idx <- which(newTxName[a] == annotData[,2])
  newGeneID[a] <- annotData[idx,1]
}

## Add the GeneIDs to the promoters GRanges
promoterData@elementMetadata@listData[["gene_id"]] <- c(newGeneID)

## Cleanup workspace
rm(promoters, a, geneKey, idx, inputFile, newGeneID, newTxName, txdb, txNames)


#### Scan for peaks in the SNU-61 replicate-merged file

## Read in the bed file
bedFile <- "C:\\Users\\jsk33\\Documents\\atac\\mydata\\snu61\\wt01\\peaks\\SNU61-WT-01.all_summits.bed"
snu61Peaks <- readBed(bedFile, track.line = FALSE, remove.unusual = FALSE, zero.based = TRUE)
## Keep standard ranges
snu61Peaks <- keepStandardChromosomes(snu61Peaks, pruning.mode="coarse")

## Find overlapping ranges between peaks and genes
overlaps <- findOverlaps(snu61Peaks, promoterData)

## How many of the overlaps referend unique annotated promoters?
idxuq <- unique(overlaps@to)
## Get a GRanges of the genes that have a peak
genePeaks <- promoterData[idxuq]

## How many of these are unique GeneIDs? 
length(unique(genePeaks@elementMetadata@listData[["gene_id"]]))
##
activeGenePromoters <- unique(genePeaks@elementMetadata@listData[["gene_id"]])
idx2 <- which(activeGenePromoters %in% expData[,3])

## Change -inf values to -11
idx3 <- which(expData[,6] == "-Inf")
expData[idx3,6] <- -12

## Sort the data
sortedExpData <- expData[order(expData$SNU61_LARGE_INTESTINE_log2),]
## Annotate the genes with peaks
sortedExpData$peak <- "RED"
idx4 <- which(expData[,3] %in% activeGenePromoters)
sortedExpData$peak[idx4] <- "BLACK"

## What does this plot tell us?
#plot(sortedExpData[,6], col = sortedExpData$peak)
plot(sortedExpData[,6])

##
hist(sortedExpData[,6], breaks = 40)

#### Make bins of the genes based on log2 expression
## Plot a histogram/dotplot overlay of the percentage of genes in that bin with a peak/ without a peak

##
binData <- list()
binData$"-11"$data <- sortedExpData[which(sortedExpData[,6] <= -10),]
binData$"-10"$data <- sortedExpData[which(sortedExpData[,6] <= -9 & sortedExpData[,6] > -10),]
binData$"-9"$data <- sortedExpData[which(sortedExpData[,6] <= -8 & sortedExpData[,6] > -9),]
binData$"-8"$data <- sortedExpData[which(sortedExpData[,6] <= -7 & sortedExpData[,6] > -8),]
binData$"-7"$data <- sortedExpData[which(sortedExpData[,6] <= -6 & sortedExpData[,6] > -7),]
binData$"-6"$data <- sortedExpData[which(sortedExpData[,6] <= -5 & sortedExpData[,6] > -6),]
binData$"-5"$data <- sortedExpData[which(sortedExpData[,6] <= -4 & sortedExpData[,6] > -5),]
binData$"-4"$data <- sortedExpData[which(sortedExpData[,6] <= -3 & sortedExpData[,6] > -4),]
binData$"-3"$data <- sortedExpData[which(sortedExpData[,6] <= -2 & sortedExpData[,6] > -3),]
binData$"-2"$data <- sortedExpData[which(sortedExpData[,6] <= -1 & sortedExpData[,6] > -2),]
binData$"-1"$data <- sortedExpData[which(sortedExpData[,6] <= 0 & sortedExpData[,6] > -1),]
binData$"0"$data <- sortedExpData[which(sortedExpData[,6] <= 1 & sortedExpData[,6] > 0),]
binData$"1"$data <- sortedExpData[which(sortedExpData[,6] <= 2 & sortedExpData[,6] > 1),]
binData$"2"$data <- sortedExpData[which(sortedExpData[,6] <= 3 & sortedExpData[,6] > 2),]
binData$"3"$data <- sortedExpData[which(sortedExpData[,6] <= 4 & sortedExpData[,6] > 3),]
binData$"4"$data <- sortedExpData[which(sortedExpData[,6] <= 5 & sortedExpData[,6] > 4),]
binData$"5"$data <- sortedExpData[which(sortedExpData[,6] <= 6 & sortedExpData[,6] > 5),]
binData$"6"$data <- sortedExpData[which(sortedExpData[,6] <= 7 & sortedExpData[,6] > 6),]
binData$"7"$data <- sortedExpData[which(sortedExpData[,6] <= 8 & sortedExpData[,6] > 7),]
binData$"8"$data <- sortedExpData[which(sortedExpData[,6] <= 9 & sortedExpData[,6] > 8),]
binData$"9"$data <- sortedExpData[which(sortedExpData[,6] <= 10 & sortedExpData[,6] > 9),]

## Calculate the ratio
binData$"-11"$ratio <- (length(which(binData[["-11"]][["data"]][["peak"]] == "RED")) / length(binData[["-11"]][["data"]][["X"]]))*100
binData$"-10"$ratio <- (length(which(binData[["-10"]][["data"]][["peak"]] == "RED")) / length(binData[["-10"]][["data"]][["X"]]))*100
binData$"-9"$ratio <- (length(which(binData[["-9"]][["data"]][["peak"]] == "RED")) / length(binData[["-9"]][["data"]][["X"]]))*100
binData$"-8"$ratio <- (length(which(binData[["-8"]][["data"]][["peak"]] == "RED")) / length(binData[["-8"]][["data"]][["X"]]))*100
binData$"-7"$ratio <- (length(which(binData[["-7"]][["data"]][["peak"]] == "RED")) / length(binData[["-7"]][["data"]][["X"]]))*100
binData$"-6"$ratio <- (length(which(binData[["-6"]][["data"]][["peak"]] == "RED")) / length(binData[["-6"]][["data"]][["X"]]))*100
binData$"-5"$ratio <- (length(which(binData[["-5"]][["data"]][["peak"]] == "RED")) / length(binData[["-5"]][["data"]][["X"]]))*100
binData$"-4"$ratio <- (length(which(binData[["-4"]][["data"]][["peak"]] == "RED")) / length(binData[["-4"]][["data"]][["X"]]))*100
binData$"-3"$ratio <- (length(which(binData[["-3"]][["data"]][["peak"]] == "RED")) / length(binData[["-3"]][["data"]][["X"]]))*100
binData$"-2"$ratio <- (length(which(binData[["-2"]][["data"]][["peak"]] == "RED")) / length(binData[["-2"]][["data"]][["X"]]))*100
binData$"-1"$ratio <- (length(which(binData[["-1"]][["data"]][["peak"]] == "RED")) / length(binData[["-1"]][["data"]][["X"]]))*100
binData$"0"$ratio <- (length(which(binData[["0"]][["data"]][["peak"]] == "RED")) / length(binData[["0"]][["data"]][["X"]]))*100
binData$"1"$ratio <- (length(which(binData[["1"]][["data"]][["peak"]] == "RED")) / length(binData[["1"]][["data"]][["X"]]))*100
binData$"2"$ratio <- (length(which(binData[["2"]][["data"]][["peak"]] == "RED")) / length(binData[["2"]][["data"]][["X"]]))*100
binData$"3"$ratio <- (length(which(binData[["3"]][["data"]][["peak"]] == "RED")) / length(binData[["3"]][["data"]][["X"]]))*100
binData$"4"$ratio <- (length(which(binData[["4"]][["data"]][["peak"]] == "RED")) / length(binData[["4"]][["data"]][["X"]]))*100
binData$"5"$ratio <- (length(which(binData[["5"]][["data"]][["peak"]] == "RED")) / length(binData[["5"]][["data"]][["X"]]))*100
binData$"6"$ratio <- (length(which(binData[["6"]][["data"]][["peak"]] == "RED")) / length(binData[["6"]][["data"]][["X"]]))*100
binData$"7"$ratio <- (length(which(binData[["7"]][["data"]][["peak"]] == "RED")) / length(binData[["7"]][["data"]][["X"]]))*100
binData$"8"$ratio <- (length(which(binData[["8"]][["data"]][["peak"]] == "RED")) / length(binData[["8"]][["data"]][["X"]]))*100
binData$"9"$ratio <- (length(which(binData[["9"]][["data"]][["peak"]] == "RED")) / length(binData[["9"]][["data"]][["X"]]))*100

## Plot it (how?)
ratios <- c(binData$"-11"$ratio, binData$"-10"$ratio, binData$"-9"$ratio, binData$"-8"$ratio, binData$"-7"$ratio, binData$"-6"$ratio,
            binData$"-5"$ratio, binData$"-4"$ratio, binData$"-3"$ratio, binData$"-2"$ratio, binData$"-1"$ratio, binData$"0"$ratio,
            binData$"1"$ratio, binData$"2"$ratio, binData$"3"$ratio, binData$"4"$ratio, binData$"5"$ratio, binData$"6"$ratio,
            binData$"7"$ratio, binData$"8"$ratio, binData$"9"$ratio)


plot(ratios)








#### Next step, subset based on only genes with peaks, then plot magnitude of peak against expression
#############
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("mygene", version = "3.8")
library(mygene)

##
gene <- getGene("1017", fields="all")
