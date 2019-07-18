# AR,SOX2,FOXM1,EZH2,MYCN,MYC
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(Biostrings))
load("C:\\Users\\jsk33\\Desktop\\lncap\\motifData.Rdata")
#x <- readRDS(outPath)

#### AR ####################################################################################
currentGene <- "AR"
com <- paste0("motifs <- motifData$", currentGene)
eval(parse(text = com))
score <- "90%"
numMotifs <- length(motifs)
genome <- Hsapiens
bindingSites <- list()
PWM <- motifs[[1]]
sites <- Biostrings::matchPWM(PWM, genome, min.score = score)
numSites <- length(sites)
tempSites <- list()
tempSites$Gene <- currentGene
tempSites$PWM <- PWM
tempSites$sites <- sites
tempSites$score <- score
tempSites$numSites <- numSites
bindingSites[[1]] <- tempSites
outPath <- "C:\\Users\\jsk33\\Desktop\\lncap\\AR.bindingSites.RDS"
saveRDS(bindingSites, file = outPath)

#### SOX2 ####################################################################################
currentGene <- "SOX2"
com <- paste0("motifs <- motifData$", currentGene)
eval(parse(text = com))
score <- "90%"
numMotifs <- length(motifs)
genome <- Hsapiens
bindingSites <- list()
PWM <- motifs[[1]]
sites <- Biostrings::matchPWM(PWM, genome, min.score = score)
numSites <- length(sites)
tempSites <- list()
tempSites$Gene <- currentGene
tempSites$PWM <- PWM
tempSites$sites <- sites
tempSites$score <- score
tempSites$numSites <- numSites
bindingSites[[1]] <- tempSites
outPath <- "C:\\Users\\jsk33\\Desktop\\lncap\\SOX2.bindingSites.RDS"
saveRDS(bindingSites, file = outPath)

#### FOXM1 ####################################################################################
currentGene <- "FOXM1"
com <- paste0("motifs <- motifData$", currentGene)
eval(parse(text = com))
score <- "95%"
numMotifs <- length(motifs)
genome <- Hsapiens
bindingSites <- list()
PWM <- motifs[[1]]
sites <- Biostrings::matchPWM(PWM, genome, min.score = score)
numSites <- length(sites)
tempSites <- list()
tempSites$Gene <- currentGene
tempSites$PWM <- PWM
tempSites$sites <- sites
tempSites$score <- score
tempSites$numSites <- numSites
bindingSites[[1]] <- tempSites
outPath <- "C:\\Users\\jsk33\\Desktop\\lncap\\FOXM1.bindingSites.RDS"
saveRDS(bindingSites, file = outPath)

#### EZH2 ####################################################################################
currentGene <- "EZH2"
com <- paste0("motifs <- motifData$", currentGene)
eval(parse(text = com))
score <- "99%"
numMotifs <- length(motifs)
genome <- Hsapiens
bindingSites <- list()
PWM <- motifs[[1]]
sites <- Biostrings::matchPWM(PWM, genome, min.score = score)
numSites <- length(sites)
tempSites <- list()
tempSites$Gene <- currentGene
tempSites$PWM <- PWM
tempSites$sites <- sites
tempSites$score <- score
tempSites$numSites <- numSites
bindingSites[[1]] <- tempSites
outPath <- "C:\\Users\\jsk33\\Desktop\\lncap\\EZH2.bindingSites.RDS"
saveRDS(bindingSites, file = outPath)

#### MYCN ####################################################################################
currentGene <- "MYCN"
com <- paste0("motifs <- motifData$", currentGene)
eval(parse(text = com))
score <- "95%"
numMotifs <- length(motifs)
genome <- Hsapiens
bindingSites <- list()
PWM <- motifs[[1]]
sites <- Biostrings::matchPWM(PWM, genome, min.score = score)
numSites <- length(sites)
tempSites <- list()
tempSites$Gene <- currentGene
tempSites$PWM <- PWM
tempSites$sites <- sites
tempSites$score <- score
tempSites$numSites <- numSites
bindingSites[[1]] <- tempSites
outPath <- "C:\\Users\\jsk33\\Desktop\\lncap\\MYCN.bindingSites.RDS"
saveRDS(bindingSites, file = outPath)

#### MYC ####################################################################################
currentGene <- "MYC"
com <- paste0("motifs <- motifData$", currentGene)
eval(parse(text = com))
score <- "95%"
numMotifs <- length(motifs)
genome <- Hsapiens
bindingSites <- list()
PWM <- motifs[[1]]
sites <- Biostrings::matchPWM(PWM, genome, min.score = score)
numSites <- length(sites)
tempSites <- list()
tempSites$Gene <- currentGene
tempSites$PWM <- PWM
tempSites$sites <- sites
tempSites$score <- score
tempSites$numSites <- numSites
bindingSites[[1]] <- tempSites
outPath <- "C:\\Users\\jsk33\\Desktop\\lncap\\MYC.bindingSites.RDS"
saveRDS(bindingSites, file = outPath)