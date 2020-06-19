# sAML are a separate condition from AML 
library(data.table)
library(maftools)
library(ggplot2)
library(ggpubr)
library(plyr)
setwd("~/mnt/workstation/mnt/data1/adam")
source("software/shiny/rnaediting/RNAEditingFunctions.R")
load("jamieson/holm/rdata/ADAR.tx.lcpm_02192020.rdata")
load("jamieson/holm/expression/DGE_MPN_sampledata.rda")
load("software/shiny/rnaediting/MAF_files.RData")

known.maf <- subset(mpn, grepl("Known", Novelty) & !(Condition_code == "Diseased") & !(Variant_Classification == "IGR"))
known.maf[grepl("Mutation|Start|Splice", Variant_Classification)]$Variant_Classification <- "Missense_Mutation"

# Figures
prog.all <- subset(known.maf, Cell.type=="Progenitor")
ggBox(prog.all)
ggsave("jamieson/holm/figures_022020/rnaedit.prog.vaf.bycondition.region_02202020.pdf")
ggViolin(prog.all)
ggsave("jamieson/holm/figures_022020/rnaedit.violin.bycondition_02202020.pdf")

# Reads per million
DGE_MPN_sampledata$Sample <- gsub("\\.", "-", gsub("X|_ACAGTG|_GTGAAA|_S[0-9]|PB.|BM.", "", DGE_MPN_sampledata$Sample))
DGE_MPN_sampledata$EditsPerMillionReads <- sapply(DGE_MPN_sampledata$Sample, function(i) nrow(subset(prog.all, Sample==i)))/(DGE_MPN_sampledata$lib.size/1000000)
DGE_MPN_sampledata <- merge(mpn.meta, DGE_MPN_sampledata, by="Sample")
ggplot(DGE_MPN_sampledata, aes(y=EditsPerMillionReads, x=Condition.x, fill=Condition.x)) + geom_boxplot() + theme_bw() + xlab("Condition") + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_fill_brewer(palette="Set2")
ggsave("jamieson/holm/figures_022020/rnaedit.editsPerMillionReads_02202020.pdf")

# Statistical tests for differences in VAF between phenotypes: Kolmogorov Smirnov Test and T Test
ks.phenotype <- sapply(levels(prog.all$Condition)[2:8], function(i) ks.test(subset(prog.all, Condition=="Aged_Normal")$VAF, subset(prog.all, Condition==i)$VAF, alternative = "greater")$p.value) 
ttest.phenotype <- sapply(levels(prog.all$Condition)[2:8], function(i) t.test(subset(prog.all, Condition=="Aged_Normal")$VAF, subset(prog.all, Condition==i)$VAF, alternative = "less")$p.value)
pv.phenotype.df <- data.frame(Kolmogorov_Smirnov=ks.phenotype, T_test=ttest.phenotype)
write.table(pv.phenotype.df, file="~/mnt/workstation/mnt/data1/adam/jamieson/holm/figures_022020/vaf_phenotype_statistical_signif_022020.tsv", sep="\t", quote=FALSE, row.names=TRUE)

ks.treatment <- sapply(levels(prog.all$Treatment_type)[2:8], function(i) ks.test(subset(prog.all, Treatment_type=="None")$VAF, subset(prog.all, Treatment_type==i)$VAF, alternative = "greater")$p.value)
ttest.treatment <- sapply(levels(prog.all$Treatment_type)[2:8], function(i) t.test(subset(prog.all, Treatment_type=="None")$VAF, subset(prog.all, Treatment_type==i)$VAF, alternative = "less")$p.value)
pv.treatment.df <- data.frame(Kolmogorov_Smirnov=ks.treatment, T_test=ttest.treatment)
write.table(pv.treatment.df, file="jamieson/holm/figures_022020/vaf_treatment_statistical_signif_022020.tsv", sep="\t", quote=FALSE, row.names=TRUE)

# Stratified by Variant Classification
ks.phenotype.by.vc <- sapply(levels(prog.all$Condition)[2:8], function(i) lapply(unique(prog.all$Variant_Classification), function(j) ks.test(subset(prog.all, Variant_Classification==j & Condition=="Aged_Normal")$VAF, subset(prog.all, Variant_Classification==j & Condition==i)$VAF, alternative = "greater")$p.value)) 
row.names(ks.phenotype.by.vc) <- unique(prog.all$Variant_Classification)
write.table(ks.phenotype.by.vc, file="jamieson/holm/figures_022020/ks_phenotype_vc_statistical_signif_022020.tsv", sep="\t", quote=FALSE, row.names=TRUE)

# ADAR expression correlation
mean.vaf <- do.call(rbind, lapply(unique(known.maf$Sample), function(i) data.frame(Sample=i, mean.vaf=mean(subset(known.maf, Sample==i)$VAF))))
adar.exp <- do.call(cbind, ADAR.tx.lcpm_02192020)
adar.exp <- as.data.frame(t(as.matrix(adar.exp)))
adar.exp$Sample <- gsub("\\.", "-", gsub("X|_ACAGTG|_GTGAAA|_S[0-9]|.PB|.BM", "", row.names(adar.exp)))
adar.exp <- merge(adar.exp, mpn.meta, by="Sample")
adar.exp <- merge(adar.exp, mean.vaf, by="Sample")
# Stem=1, square Progenitor=2, triangle
adar.exp$Cell.type <- as.numeric(factor(adar.exp$Cell.type, levels=c("Stem", "Progenitor")))
adar.exp$Condition <- factor(adar.exp$Condition, levels=c("Aged_Normal", "Young_Normal", "ET", "PV", "MF", "CML", "sAML", "denovoAML"))
# both stem and prog
adar.exp.p150gt0 <- subset(adar.exp, ENST00000368474 > 0)
adar.exp.p110gt0 <- subset(adar.exp, ENST00000368471 > 0)
ggscatter(adar.exp.p150gt0, x="ENST00000368474", y="mean.vaf", color="Condition", size=c(4, 4)[adar.exp.p150gt0$Cell.type], shape=c(15, 17)[adar.exp.p150gt0$Cell.type], palette = "Set2") + ylab("Mean RNA editing VAF") + xlab("ADAR p150 log CPM expression") + guides(colour = guide_legend(override.aes = list(size=4))) + ylim(c(0,0.6)) + xlim(c(0,10))
ggsave("jamieson/holm/figures_022020/rnaedit.p150.p150gt0.meanvaf.corr_05182020.pdf")
ggscatter(adar.exp.p150gt0, x="ENST00000368471", y="mean.vaf", color="Condition", size=c(4, 4)[adar.exp.p150gt0$Cell.type], shape=c(15, 17)[adar.exp.p150gt0$Cell.type], palette = "Set2") + ylab("Mean RNA editing VAF") + xlab("ADAR p150 log CPM expression") + guides(colour = guide_legend(override.aes = list(size=4))) + ylim(c(0,0.6)) + xlim(c(0,10))
ggsave("jamieson/holm/figures_022020/rnaedit.p110.p150gt0.meanvaf.corr_05182020.pdf")
ggscatter(adar.exp.p110gt0, x="ENST00000368471", y="mean.vaf", color="Condition", size=c(4, 4)[adar.exp.p110gt0$Cell.type], shape=c(15, 17)[adar.exp.p110gt0$Cell.type], palette = "Set2") + ylab("Mean RNA editing VAF") + xlab("ADAR p110 log CPM expression") + guides(colour = guide_legend(override.aes = list(size=4))) + ylim(c(0,0.6)) + xlim(c(0,10))
ggsave("jamieson/holm/figures_022020/rnaedit.p110.p110gt0.meanvaf.corr_05182020.pdf")
# Correlations
corstem <- cor.test(formula = ~ ENST00000368474 + mean.vaf, data = subset(adar.exp.p150gt0, Cell.type==1))
corprog <- cor.test(formula = ~ ENST00000368474 + mean.vaf, data = subset(adar.exp.p150gt0, Cell.type==2))
corstem <- cor.test(formula = ~ ENST00000368471 + mean.vaf, data = subset(adar.exp.p110gt0, Cell.type==1))
corprog <- cor.test(formula = ~ ENST00000368471 + mean.vaf, data = subset(adar.exp.p110gt0, Cell.type==2))

# CDK13 plots
cdk13.maf <- subset(prog.all, Hugo_Symbol=="CDK13")
cdk13.maf$Hugo_Symbol <- cdk13.maf$Protein_Change
cdk13 <- read.maf(cdk13.maf, clinicalData = mpn.meta)
risk.colors <- RColorBrewer::brewer.pal(n = 8,name = "Paired")
names(risk.colors) <- c("Aged_Normal", "Young_Normal", "ET", "PV", "MF", "CML", "sAML_UnTx", "sAML_Tx")
risk.order <- as.numeric(factor(names(risk.colors), levels=c("Aged_Normal", "Young_Normal", "ET", "PV", "MF", "CML", "sAML_UnTx", "sAML_Tx")))
risk.colors = list(Condition_code = risk.colors)
dev.off(dev.list()["RStudioGD"])
oncoplot(maf =cdk13,
         clinicalFeatures = "Condition_code", sortByAnnotation = TRUE, 
         annotationOrder = c("Aged_Normal", "Young_Normal", "ET", "PV", "MF", "CML", "sAML_Tx", "sAML_UnTx"),
         annotationColor = risk.colors, drawColBar = FALSE)

# STAT plots
stat <- fread("jamieson/holm/expression/ADAR.tx.lcpm.exp.rdataMPN_STAT3_STAT5B_editingSites_AML_AN.csv", stringsAsFactors = FALSE)
stat$Hugo_Symbol <- paste(stat$Hugo_Symbol, stat$Location)
stat.in <- read.maf(stat, vc_nonSyn = unique(stat$Variant_Classification), clinicalData = mpn.meta)
oncoplot(stat.in, clinicalFeatures = "Condition_code")
stat <- fread("jamieson/holm/vcf/rnaedit/STAT3_STAT5B_rnaedit_sites_hg38pos_05292019.maf", stringsAsFactors = FALSE, sep = "\t")
stat <- subset(stat, Condition_code %in% c("Aged_Normal", "AML"))
stat$Hugo_Symbol <- paste0("chr", stat$Chromosome, ":", stat$Start_position)
#stat$published <- ifelse(stat$hg38_position %in% c(42318026, 42318064, 42318069, 42318091), yes="published", no="novel")
#stat.known <- subset(stat, published=="published" & Condition_code != "MDS")
stat.in <- read.maf(stat, vc_nonSyn = unique(stat$Variant_Classification))
clinicalData <- unique(stat[,c("Tumor_Sample_Barcode", "Condition_code")])
tryCatch({dev.off(dev.list()["RStudioGD"])}, error=function(e){})
oncoplot(stat.in, annotationDat = clinicalData, clinicalFeatures = "Condition_code", 
         genes=c("chr17:40470044", "chr17:40470109"), annotationOrder = c("AML", "Aged_Normal"), drawColBar = FALSE,
         sortByAnnotation = TRUE, showTumorSampleBarcodes = TRUE, removeNonMutated = FALSE)


