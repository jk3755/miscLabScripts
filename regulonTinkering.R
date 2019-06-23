#2018 Nov 6
#code to tinker with interactome
library(viper)
setwd("~/pCloud Drive/pCloud Sync/Collaborations/Jordan/misc")
# save(dset,regulon_motifMI, file = "dset_motifMIregulon.RData")
load(file = "dset_motifMIregulon.RData")#get the dset and regulon_motifMI
geneId <- "4089" #smad4
interactome_tinker <- read.table(file = "network_coad_tfCotf_3col_sorted.txt", sep = "\t")
regulon_req <- interactome_tinker[which(interactome_tinker[,1]==geneId),]#getting req regulon to tinker
interactome_tinker <- interactome_tinker[-(which(interactome_tinker[,1]==geneId)),]# removing the req regulon from interactome
colnames(regulon_motifMI) <- colnames(regulon_req)
regulon_req <- rbind(regulon_req,regulon_motifMI)
regulon_req <- regulon_req[order(regulon_req[,3], decreasing = TRUE),]
write.table(interactome_tinker, file = "network_coad_tfCotf_3col_tinkered.txt", sep = "\t", append = FALSE, row.names = FALSE, col.names = FALSE)
write.table(regulon_req, file = "network_coad_tfCotf_3col_tinkered.txt", sep = "\t", append = TRUE, row.names = FALSE, col.names = FALSE)
regul <- aracne2regulon("network_coad_tfCotf_3col_tinkered.txt", dset, format = "3col", verbose = TRUE)#regulon
print(regul)
coadViper <- viper(dset, regul, method = "scale") #the viper algorithm
dim(coadViper)