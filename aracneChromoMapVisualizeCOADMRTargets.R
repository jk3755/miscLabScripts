## Packages
#install.packages("BiocManager")
#BiocManager::install("aracne.networks", version = "3.8")
#BiocManager::install("mygene", version = "3.8")
#BiocManager::install("chromoMap", version = "3.8")
library(ATACseqQC)
library(GenomicRanges)
library(aracne.networks)
library(mygene)
library(chromoMap)

## Get the COAD interactome
coadInteractome <- aracne.networks::reguloncoad
## List of COAD MRs (with motif binding sites)
## Entrez IDs in parentheses
## ASCL2 (430), ESRRA (2101), TCF7 (6932), POU5F1B (5462), HNF4A (3172), OVOL1 (5017)
## GMEB2 (26205), CBFA2T2 (9139), HOXA3 (3200), MNX1 (3110), ZSWIM1 (90204), CDX2 (1045)
ascl2Targets <- names(coadInteractome[["430"]][["tfmode"]])
esrraTargets <- names(coadInteractome[["2101"]][["tfmode"]])
tcf7Targets <- names(coadInteractome[["6932"]][["tfmode"]])
pou5f1b2Targets <- names(coadInteractome[["5462"]][["tfmode"]])
hnf4aTargets <- names(coadInteractome[["3172"]][["tfmode"]])
ovol1Targets <- names(coadInteractome[["5017"]][["tfmode"]])
gmeb2Targets <- names(coadInteractome[["26205"]][["tfmode"]])
cbfa2t2Targets <- names(coadInteractome[["9139"]][["tfmode"]])
hoxa3Targets <- names(coadInteractome[["3200"]][["tfmode"]])
mnx1Targets <- names(coadInteractome[["3110"]][["tfmode"]])
zswim1Targets <- names(coadInteractome[["90204"]][["tfmode"]])
cdx2Targets <- names(coadInteractome[["1045"]][["tfmode"]])

## Build functions
showTargetsChromoMap <- function(targets){

  numTargets <- length(targets)
  ## Retrieve annotation info for the regulator targets
  TargetsList <- list()
  for (a in 1:length(targets)){
    TargetsList[a] <- mygene::getGene(geneid = targets[a], fields = "all")}
  
  # count the number of gene locations we have
  locCount <- 0
  for (a in 1:numTargets){
    if (is.null(TargetsList[[a]][["genomic_pos"]])){next} else {
      if (is.list(TargetsList[[a]][["genomic_pos"]][[1]])){
        locCount <- (locCount + length(TargetsList[[a]][["genomic_pos"]]))
      } else {
        locCount <- (locCount + 1)}}}
  
  # Retrieve the genomic coordinates of all network targets
  # retrieve all genomic coords here, can later trim non standard ones more easily with GRanges functions
  TargetLocations <- matrix(data = NA, nrow = locCount, ncol = 4)
  colnames(TargetLocations) <- c("gene", "chr", "start", "end")
  idx <- 1
  for (a in 1:numTargets){
    if (is.null(TargetsList[[a]][["genomic_pos"]])){next} else {
      if (is.list(TargetsList[[a]][["genomic_pos"]][[1]])){
        for (b in 1:length(TargetsList[[a]][["genomic_pos"]])){
          TargetLocations[idx,1] <- TargetsList[[a]][["symbol"]]
          TargetLocations[idx,2] <- paste0("chr", TargetsList[[a]][["genomic_pos"]][[b]][["chr"]])
          TargetLocations[idx,3] <- TargetsList[[a]][["genomic_pos"]][[b]][["start"]]
          TargetLocations[idx,4] <- TargetsList[[a]][["genomic_pos"]][[b]][["end"]]
          idx <- (idx+1)
        } # end for (b in 1:length(target_list[[a]][["genomic_pos"]]))
      } else {
        TargetLocations[idx,1] <- TargetsList[[a]][["symbol"]]
        TargetLocations[idx,2] <- paste0("chr", TargetsList[[a]][["genomic_pos"]][["chr"]])
        TargetLocations[idx,3] <- TargetsList[[a]][["genomic_pos"]][["start"]]
        TargetLocations[idx,4] <- TargetsList[[a]][["genomic_pos"]][["end"]]
        idx <- (idx+1)}}}
  
  ## Convert the genomic locations of the targets into GRanges
  TargetLocationsGR <- GRanges(
    seqnames = TargetLocations[,2],
    ranges = IRanges(start = as.numeric(TargetLocations[,3]), end = as.numeric(TargetLocations[,4])))
  #strand = c(rep("+", times = num_names)),
  #seqinfo = Seqinfo(genome="hg38"),
  #score = c(rep(1, times = num_names)))
  ## prune to standard xsomes
  TargetLocationsGR <- keepStandardChromosomes(TargetLocationsGR, pruning.mode="coarse")
  
  
  ## Visualize with chromoMap
  ## If it doesnt show up, try reloading the library
  canno <- data.frame(
    name = c(1:length(TargetLocationsGR@elementMetadata@nrows)),
    chrom = as.vector(TargetLocationsGR@seqnames, mode="any"),
    start = TargetLocationsGR@ranges@start,
    data = TargetLocationsGR@ranges@width)
  
  return(canno)
}

## Get annotation objects
ascl2Anno <- showTargetsChromoMap(ascl2Targets)
esrraAnno <- showTargetsChromoMap(esrraTargets)
tcf7Anno <- showTargetsChromoMap(tcf7Targets)
pou5f1bAnno <- showTargetsChromoMap(pou5f1b2Targets)
hnf4aAnno <- showTargetsChromoMap(hnf4aTargets)
ovol1Anno <- showTargetsChromoMap(ovol1Targets)
gmeb2Anno <- showTargetsChromoMap(gmeb2Targets)
cbfa2t2Anno <- showTargetsChromoMap(cbfa2t2Targets)
hoxa3Anno <- showTargetsChromoMap(hoxa3Targets)
mnx1Anno <- showTargetsChromoMap(mnx1Targets)
zswim1Anno <- showTargetsChromoMap(zswim1Targets)
cdx2Anno <- showTargetsChromoMap(cdx2Targets)

## Build chromomap objects
cmapASCL2 <- chromoMap(ascl2Anno, type = "annotation")
cmapESRRA <- chromoMap(esrraAnno, type = "annotation")
cmapTCF7 <- chromoMap(tcf7Anno, type = "annotation")
cmapPOU5F1B <- chromoMap(pou5f1bAnno, type = "annotation")
cmapHNF4A <- chromoMap(hnf4aAnno, type = "annotation")
cmapOVOL1 <- chromoMap(ovol1Anno, type = "annotation")
cmapGMEB2 <- chromoMap(gmeb2Anno, type = "annotation")
cmapCBFA2T2 <- chromoMap(cbfa2t2Anno, type = "annotation")
cmapHOXA3 <- chromoMap(hoxa3Anno, type = "annotation")
cmapMNX1 <- chromoMap(mnx1Anno, type = "annotation")
cmapZSWIM1 <- chromoMap(zswim1Anno, type = "annotation")
cmapCDX2 <- chromoMap(cdx2Anno, type = "annotation")

## Visualize the chromomaps
## If it doesn't show up, try loading library again
library(chromoMap)
cmapASCL2
cmapESRRA
cmapTCF7
cmapPOU5F1B
cmapHNF4A
cmapOVOL1
cmapGMEB2
cmapCBFA2T2
cmapHOXA3
cmapMNX1
cmapZSWIM1
cmapCDX2

