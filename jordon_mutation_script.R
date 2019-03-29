#Hi Jordon,
#I used maftools and lollipopPlot. Please find attached the sample script. Just pass a maf file to it.
#Many thanks!
#Som

if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("maftools")

library(maftools)
# load and read MAF file

luad <- read.maf(maf = "gsc_LUAD_pairs.aggregated.capture.tcga.uuid.automated.somatic.maf")

#MAF object
luad

#Shows sample summry.
getSampleSummary(luad)

#Shows gene summary.
getGeneSummary(luad)

#shows clinical data associated with samples
getClinicalData(luad)

#plot samples with mutations

#Lollipop plots for amino acid changes
#Labelling points.
pdf("lollipop_luad_EGFR_labels.pdf")
lollipopPlot(maf = luad, gene = 'EGFR', AACol = 'Protein_Change', labelPos = 'all')
dev.off()
