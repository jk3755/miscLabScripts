
## Install the packages you will need
source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges", suppressUpdates = TRUE)
biocLite("EnsDb.Hsapiens.v86", suppressUpdates = TRUE) # This is the database from which you will retrieve the human gene annotations. An alternative would be to use TxDb.Hsapiens.UCSC.hg38.knownGene
biocLite("genomation", suppressUpdates = TRUE) # This package allows you to read in a .bed directly to a GRanges object

## Load the packages
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(genomation)

## Step 1 - load the .bed file and convert to GRanges
bedFile <- "C:\\Users\\Jordan\\Documents\\atac\\SNU61-WT-01.all_peaks.narrowPeak"
snu61Peaks <- readBed(bedFile, track.line = FALSE, remove.unusual = FALSE, zero.based = TRUE)
## Remove the entries that are not on standard chromosomes (chr1-22, X, Y), often helps prevent problems downstream
snu61Peaks <- keepStandardChromosomes(snu61Peaks, pruning.mode="coarse")

## Step 2 - Retrieve human gene annotations in GRanges format
edb <- EnsDb.Hsapiens.v86
## Annoyingly, the gene database annotates the chromosomes differently than the standard "chr1, etc...", using only a single integer
## to represent each chromosome in the GRanges object. In order to find the overlaps, you will have to edit the chromosome annotations first
seqlevelsStyle(edb) <- "UCSC"
## Now, retrieve the gene annotations from the database
hg38Genes <- genes(edb)
## Again helpful to keep only the standard entries
hg38Genes <- keepStandardChromosomes(hg38Genes, pruning.mode="coarse")

## Step 3 - Extend the ranges of the peaks in the GRanges objects (+/- 50,000-100,000 bp, depending on how far you want to search from each peak)
## Note that this will almost certainly result in ranges that extend beyond the ends of the chromosome
## This isn't a huge deal, but it will return warnings to you at some point
## There is a method to remove the out of bounds ranges but I can't find my code for it at the moment
extPeaks <- promoters(snu61Peaks, upstream = 50000, downstream = 50000, use.names=TRUE)

## Step 4 - Find the overlaps between the extended peaks and the gene annotations
## This will return a 'SortedByQueryHits' object that gives you the indices of all overlaps
## If you do View(peakGeneOverlaps), you can see a 'from' entry and a 'to' entry'
## The values in the 'from' entry refer to the indices in the first GRanges object (extPeaks)
## And the values in the 'to' entry refer to the indices in the second (hg38Genes)
peakGeneOverlaps <- findOverlaps(extPeaks, hg38Genes)



