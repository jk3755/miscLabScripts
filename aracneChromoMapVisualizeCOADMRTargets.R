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


## Retrieve annotation info for the regulator targets
ascl2TargetsList <- list()
for (a in 1:length(ascl2Targets)){
  ascl2TargetsList[a] <- mygene::getGene(geneid = ascl2Targets[a], fields = "all")}

# count the number of gene locations we have
loc_count <- 0
for (a in 1:length(mnx1_targets)){
  if (is.null(target_list[[a]][["genomic_pos"]])){next} else {
    if (is.list(target_list[[a]][["genomic_pos"]][[1]])){
      loc_count <- (loc_count + length(target_list[[a]][["genomic_pos"]]))
    } else {
      loc_count <- (loc_count + 1)}}}

# Retrieve the genomic coordinates of all network targets
# retrieve all genomic coords here, can later trim non standard ones more easily with GRanges functions
target_locations <- matrix(data = NA, nrow = loc_count, ncol = 4)
colnames(target_locations) <- c("gene", "chr", "start", "end")
idx <- 1
#
for (a in 1:loc_count){
  
  if (is.null(target_list[[a]][["genomic_pos"]])){next} else {
    
    if (is.list(target_list[[a]][["genomic_pos"]][[1]])){
      
      for (b in 1:length(target_list[[a]][["genomic_pos"]])){
        target_locations[idx,1] <- target_list[[a]][["symbol"]]
        target_locations[idx,2] <- paste0("chr", target_list[[a]][["genomic_pos"]][[b]][["chr"]])
        target_locations[idx,3] <- target_list[[a]][["genomic_pos"]][[b]][["start"]]
        target_locations[idx,4] <- target_list[[a]][["genomic_pos"]][[b]][["end"]]
        idx <- (idx+1)
      } # end for (b in 1:length(target_list[[a]][["genomic_pos"]]))
      
    } else {
      
      target_locations[idx,1] <- target_list[[a]][["symbol"]]
      target_locations[idx,2] <- paste0("chr", target_list[[a]][["genomic_pos"]][["chr"]])
      target_locations[idx,3] <- target_list[[a]][["genomic_pos"]][["start"]]
      target_locations[idx,4] <- target_list[[a]][["genomic_pos"]][["end"]]
      idx <- (idx+1)}}}

## Convert the genomic locations of the targets into GRanges
gr <- GRanges(
  seqnames = target_locations[,2],
  ranges = IRanges(start = as.numeric(target_locations[,3]), end = as.numeric(target_locations[,4])))
#strand = c(rep("+", times = num_names)),
#seqinfo = Seqinfo(genome="hg38"),
#score = c(rep(1, times = num_names)))

# prune to standard xsomes
gr <- keepStandardChromosomes(gr, pruning.mode="coarse")


## Intersect the targets GRanges with the binding sites list
intersection <- intersect(gr, wg_sites)
canno <- data.frame(
  name = c(1:148),
  chrom = as.vector(gr@seqnames, mode="any"),
  start = gr@ranges@start,
  data = gr@ranges@width)


cmap <- chromoMap(canno, type = "annotation")
cmap <- chromoMap(gr, type ="heatmap-single", HeatColRange = c("blue","white","red"))



gr1 <- GRanges(
  seqnames = c("chr1", "chr2"),
  ranges = IRanges(start = c(1000, 1000), end = c(2000, 2000))
)

gr2 <- GRanges(
  seqnames = c("chr1", "chr3"),
  ranges = IRanges(start = c(1900, 1000), end = c(3000, 2000))
)

gr3 <- intersect(gr1, gr2)
