# RNA Edit Known and Novel Combined Figure
library(data.table)
library(maftools)
load("~/mnt/workstation/mnt/data1/adam/jamieson/holm/rdata/MAF_files.RData")

maf <- subset(mpn, Cell.type %in% c("Progenitor", "Stem") & Variant_Classification == "Missense_Mutation" & Condition != "Diseased")
maf <- subset(maf, Cell.type == "Progenitor")
mpn.meta <- subset(mpn.meta, Cell.type == "Progenitor")
maf[grepl("Known", maf$Novelty) & !grepl("Non", maf$Novelty)]$Variant_Classification <- "Reported Alu"
maf[!grepl("Known", maf$Novelty) & !grepl("Non", maf$Novelty)]$Variant_Classification <- "Unreported Alu"
maf[grepl("Known", maf$Novelty) & grepl("Non", maf$Novelty)]$Variant_Classification <- "Reported Nonalu"
maf[!grepl("Known", maf$Novelty) & grepl("Non", maf$Novelty)]$Variant_Classification <- "Unreported Nonalu"

# oncoplot
Condition.colors <- RColorBrewer::brewer.pal(n = 8, name = 'Set2')
names(Condition.colors) <- c("Aged_Normal", "Young_Normal", "ET", "PV", "MF", "CML", "sAML_Tx", "sAML_UnTx")
colors <- list(Condition = Condition.colors)
vc.colors <- RColorBrewer::brewer.pal(n = 4, name = 'Set1')
names(vc.colors) <- c("Reported Alu", "Unreported Alu", "Reported Nonalu", "Unreported Nonalu")
maf.in <- read.maf(maf, clinicalData = mpn.meta, vc_nonSyn = c("Reported Alu", "Unreported Alu", "Reported Nonalu", "Unreported Nonalu"))
dev.off(dev.list()["RStudioGD"])
oncoplot(maf = maf.in, clinicalFeatures = c("Condition"), sortByAnnotation = TRUE, 
         showTumorSampleBarcodes = TRUE, annotationColor = colors, top = 25, fontSize = 0.6,
         colors = vc.colors, annotationOrder = c("Aged_Normal", "Young_Normal", "ET", "PV", "MF", "CML", "sAML_UnTx", "sAML_Tx"),
      bgCol="white", removeNonMutated = FALSE, barcode_mar = 6, gene_mar = 6, fill = TRUE)
pdf("~/mnt/workstation/mnt/data1/adam/jamieson/holm/figures_022020/RNAEdit.known.novel.missense.prog.combined_02202020.pdf")
