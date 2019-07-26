library(data.table)
library(maftools)
library(ggplot2)
library(ggpubr)
library(plyr)

# Plot
give.n <- name <- function(variables) {
  return(c(y = median(x)*1.05, label = length(x)))
  # experiment with the multiplier to find the perfect position
}

# Combined Meta Data
comb.meta <- read.csv("~/mnt/workstation/mnt/data1/adam/jamieson/holm/DNA_RNA_combined_meta06072018.csv", stringsAsFactors = FALSE, check.names = FALSE, sep = ",")
comb.meta$Sample <- gsub("\\.", "-", gsub("X|_ACAGTG|_GTGAAA|_S[0-9]|PB.|BM.", "", comb.meta$Sample))
comb.meta <- unique(comb.meta[c("Sample", "Sample_type", "Condition", "Cell_type", "Patient_ID", "Mutation", "Treatment_type")])
comb.meta$Tumor_Sample_Barcode <- gsub("A|B", "", sapply(strsplit(comb.meta$Sample, "-"), `[`, 1))
comb.meta$Sample_type <- factor(as.character(comb.meta$Sample_type),
                                levels = c("Aged normal bone marrow", "Young normal bone marrow",
                                           "PV", "ET", "MF", "CML", "AML"))
comb.meta$Condition <- factor(as.character(comb.meta$Condition),
                              levels=c("Normal", "ET", "PV", "LR_MF_PostET", "Int_1_MF", "Int_2_MF",
                                       "Int_2_MF_PostET", "Int_2_MF_PostPV", "HR_MF", "HR_MF_PostPV",
                                       "sAML", "AML", "AP_CML", "BC_CML", "CP_CML", "CML", "Diseased"))
comb.meta$Treatment_type <- factor(as.character(comb.meta$Treatment_type),
                                   levels=c("None", "Untreated", "Jak2 Inhibitor", "SHH treated", "hydroxyurea", "TKI", "vidaza"))

# DNA Meta
meta <- read.csv("~/mnt/workstation/mnt/data1/adam/jamieson/holm/somatic_metadata_0329_2018.txt", sep = "\t", stringsAsFactors = FALSE)
names(meta)[1] <- "Tumor_Sample_Barcode"
meta[grepl("CML", meta$Diagnosis),]$Diagnosis <- "CML"
meta$Tumor_Sample_Barcode <- as.character(meta$Tumor_Sample_Barcode)
#meta[(meta$Treatment_Type != ""),]$Diagnosis <- paste(meta[(meta$Treatment_Type != ""),]$Diagnosis, meta[(meta$Treatment_Type != ""),]$Treatment, sep = "_")

# RNA Meta
rna.meta <- read.csv("~/mnt/workstation/mnt/data1/adam/jamieson/holm/fheditsHolm_Jamieson_RNAseq_with_controls_meta_20180326.csv", stringsAsFactors = FALSE, check.names = FALSE)
rna.meta$Sample <- gsub("\\.", "-", gsub("X|_ACAGTG|_GTGAAA|_S[0-9]|PB.|BM.", "", rna.meta$Sample))
rna.meta$Tumor_Sample_Barcode <- gsub("A|B", "", sapply(strsplit(rna.meta$Sample, "-"), `[`, 1))
rna.meta <- unique(rna.meta[c("Sample", "Condition_code", "Condition_code2", "Treatment_code2", "Cell_type", "Mutation", "Patient_ID")])
rna.meta[grepl("CML", rna.meta$Condition_code),]$Condition_code <- "CML"
rna.meta[grepl("AML", rna.meta$Condition_code),]$Condition_code <- "AML"
rna.meta[grepl("CML", rna.meta$Condition_code2),]$Condition_code2 <- "CML"
rna.meta[grepl("AML", rna.meta$Condition_code2),]$Condition_code2 <- "AML"

# Strand Data
strand <- unique(fread("~/mnt/workstation/mnt/data1/adam/software/annotation/gene_strand_hg19.txt", stringsAsFactors = FALSE, sep = "\t", header=FALSE))
names(strand) <- c("Hugo_Symbol", "Strand")

# ADAR p150 quantile data
load("~/mnt/workstation/mnt/data1/tomw/Holm_Jamieson_Analysis/ADARp150stemexp.rdata")
load("~/mnt/workstation/mnt/data1/tomw/Holm_Jamieson_Analysis/ADARp150progexp.rdata")
adar <- rbind(ADAR.p150.Stem.quartile.df, ADAR.p150.Prog.quartile.df)
adar <- subset(adar, !(grepl(".BM.", adar$Sample)))
adar$Sample <- gsub("\\.", "-", gsub("X|_ACAGTG|_GTGAAA|_S[0-9]|.PB", "", adar$Sample))
names(adar) <- c("Sample", "p150")

# RNA Edit data
names <- c("Hugo_Symbol", "Chromosome", "Start_position", "End_position", "Variant_Classification", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "dbSNP_RS", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", "Protein_Change", "t_alt_count", "t_ref_count", "1000gp3_AF", "ExAC_AF", "dbNSFP_CADD_phred")
# Known MPN
# The file "rnaedit.95.known.alu.merged.oncotator.filt.maf" is equivalent to the file "disease_alu_merged_known_rnaedit_sites.oncotator.filt.maf" produced in the RNA Editing Analysis methods
known.alu <- fread("~/mnt/jam_db/projects/holm/vcf/rnaedit/rnaedit.95.known.alu.merged.oncotator.filt.maf", stringsAsFactors=FALSE, sep="\t", skip=107)
names(known.alu) <- names
known.alu$Novelty <- "Known Alu"
# The file "rnaedit.95.known.nonalu.merged.oncotator.filt.maf" is equivalent to the file "disease_nonalu_merged_known_rnaedit_sites.oncotator.filt.maf" produced in the RNA Editing Analysis methods
known.nonalu <- fread("~/mnt/jam_db/projects/holm/vcf/rnaedit/rnaedit.95.known.nonalu.merged.oncotator.filt.maf", stringsAsFactors=FALSE, sep="\t", skip=107)
names(known.nonalu) <- names
known.nonalu$Novelty <- "Known NonAlu"
mpn <- rbind(known.alu, known.nonalu)
mpn$Sample <- gsub("-BM|-PB|_S[0-9]", "", mpn$Matched_Norm_Sample_Barcode)
mpn$Tumor_Sample_Barcode <- as.integer(substr(mpn$Matched_Norm_Sample_Barcode, 1, 3))

# Known normal
# The file "rnaedit.jiang.alu.merged.oncotator.filt.maf" is equivalent to the file "normal_alu_merged_known_rnaedit_sites.oncotator.filt.maf" produced in the RNA Editing Analysis methods
normal.known.alu <- fread("~/mnt/jam_db/projects/jiang/vcf/rnaedit/aging/rnaedit.jiang.alu.merged.oncotator.filt.maf", stringsAsFactors=FALSE, sep="\t", skip=104)
names(normal.known.alu) <- names
normal.known.alu$Novelty <- "Known Alu"
# The file "rnaedit.jiang.nonalu.merged.oncotator.filt.maf" is equivalent to the file "normal_nonalu_merged_known_rnaedit_sites.oncotator.filt.maf" produced in the RNA Editing Analysis methods
normal.known.nonalu <- fread("~/mnt/jam_db/projects/jiang/vcf/rnaedit/aging/rnaedit.jiang.nonalu.merged.oncotator.filt.maf", stringsAsFactors=FALSE, sep="\t", skip=106)
names(normal.known.nonalu) <- names
normal.known.nonalu$Novelty <- "Known NonAlu"
normal <- rbind(normal.known.alu, normal.known.nonalu)
normal$Sample <- gsub("_ACAGTG|_S[0-9]|_GTGAAA", "", normal$Matched_Norm_Sample_Barcode)
normal$Tumor_Sample_Barcode <- normal$Sample
normal <- subset(normal, Tumor_Sample_Barcode %in% subset(rna.meta, Condition_code2 %in% c("Aged_Normal", "Young_Normal"))$Sample)


# Merge all deduplicate
maf <- rbind(mpn, normal)
maf <- unique(maf)

# ADAR specific
rnaedit.known.maf <- merge(maf, strand, by = "Hugo_Symbol")
rnaedit.known.maf$mutsig <- paste0(rnaedit.known.maf$Reference_Allele, ">", rnaedit.known.maf$Tumor_Seq_Allele2)
rnaedit.known.maf <- subset(rnaedit.known.maf, (mutsig == "A>G" & Strand=="+") | (mutsig == "T>C" & Strand=="-"))
rnaedit.known.maf$VAF <- rnaedit.known.maf$t_alt_count/(rnaedit.known.maf$t_alt_count + rnaedit.known.maf$t_ref_count)
rnaedit.known.maf <- subset(rnaedit.known.maf, (t_alt_count >= 2 & (t_alt_count + t_ref_count >= 10)) & VAF > 0.1 & VAF < 0.95)
rnaedit.known.maf <- merge(rnaedit.known.maf, rna.meta, by = "Sample")
rnaedit.known.maf$Condition_code <- factor(rnaedit.known.maf$Condition_code, levels=c("Normal", "PV", "ET",
                                                              "myelofibrosis", "high_risk_myelofibrosis", "int2_myelofibrosis", "int2_myelofibrosis_postET",
                                                              "int2_myelofibrosis_postPV", "myelofibrosis_fromMDS", "myelofibrosis_postET", "myelofibrosis_postPV",
                                                              "systemic_mastocytosis_myelofibrosis", "CML", "AML"))
rnaedit.known.maf$Condition_code2 <- factor(rnaedit.known.maf$Condition_code2, levels=c("Aged_Normal", "Young_Normal", "ET", "PV", "MF", "CML", "AML"))
#rnaedit.known.maf <- merge(rnaedit.known.maf, comb.meta, by = "Sample")

write.table(rnaedit.known.maf, file="~/mnt/jam_db/projects/holm/vcf/rnaedit/RNAedit_known_Normal_MPN_merged_filtered_09172018.maf", quote = FALSE, sep = "\t", row.names = FALSE)


prog.all <- subset(rnaedit.known.maf, Cell_type=="Progenitor")
ggplot(rnaedit.known.maf, aes(x=Variant_Classification, y=VAF, fill=Condition_code2)) + geom_boxplot()
#The below line generates figure 4a
ggviolin(rnaedit.known.maf, x="Condition_code2", y="VAF", fill = "Condition_code2", palette = "npg") + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.5)
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures/rnaedit.violin.bycondition.png")
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures/rnaedit.violin.bycondition.pdf")
ggplot(prog.all, aes(VAF, colour=Condition_code2)) + stat_ecdf(size=0.75, show.legend = TRUE) + theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures/rnaedit.ecdf.bycondition.png")
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures/rnaedit.ecdf.bycondition.pdf")
ggplot(prog.all, aes(x=Variant_Classification, y=VAF, fill=Condition_code2)) + geom_boxplot() + theme_classic() + theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures/rnaedit.boxplot.region.bycondition.png")
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures/rnaedit.boxplot.region.bycondition.pdf")
ggplot(prog.all, aes(VAF, colour=Condition_code2)) + facet_wrap(~Variant_Classification) + stat_ecdf(size=0.75, show.legend = TRUE) + theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures/rnaedit.ecdf.bycondition_code2.facetbyvc.png")
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures/rnaedit.ecdf.bycondition_code2.facetbyvc.pdf")
ggplot(prog.all, aes(VAF, colour=Condition_code2)) + stat_ecdf(size=0.75, show.legend = TRUE) + theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures/rnaedit.ecdf.bycondition_code2.png")
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures/rnaedit.ecdf.bycondition_code2.pdf")
ggplot(prog.all, aes(VAF, colour=Treatment_type)) + stat_ecdf(size=0.75, show.legend = TRUE) + theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures/rnaedit.ecdf.bytreatment.png")
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures/rnaedit.ecdf.bytreatment.pdf")
# The below line generates extended data figure 4b
ggboxplot(subset(rnaedit.known.maf, Cell_type.x=="Progenitor"), x="Condition", y="VAF", fill="Condition", title = "Progenitor Cells") + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + xlab(label = "Condition") + scale_fill_discrete(name="Condition")
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures/rnaedit.vafbycondition.prog.png")
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures/rnaedit.vafbycondition.prog.pdf")
ggboxplot(subset(rnaedit.known.maf, Cell_type.x=="Stem"), x="Condition", y="VAF", fill="Condition", title = "Stem Cells") + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + xlab(label = "Condition") + scale_fill_discrete(name="Condition")
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures/rnaedit.vafbycondition.stem.png")
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures/rnaedit.vafbycondition.stem.pdf")
ggboxplot(subset(rnaedit.known.maf, Cell_type.x=="Progenitor"), x="Treatment_type", y="VAF", fill="Treatment_type", title = "Progenitor Cells") + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + xlab(label = "Treatment") + scale_fill_discrete(name="Treatment")
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures/rnaedit.vafbytreatment.prog.png")
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures/rnaedit.vafbytreatment.prog.pdf")
ggboxplot(subset(rnaedit.known.maf, Cell_type.x=="Stem"), x="Treatment_type", y="VAF", fill="Treatment_type", title = "Stem Cells") + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + xlab(label = "Treatment") + scale_fill_discrete(name="Treatment")
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures/rnaedit.vafbytreatment.stem.png")
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures/rnaedit.vafbytreatment.stem.pdf")
maf.p150 <- merge(rnaedit.known.maf, adar, by = "Sample")
ggplot(maf.p150, aes(x=Variant_Classification, y=VAF, fill=p150)) + geom_boxplot() + theme_classic()
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures/rnaedit.p150byregion.png")
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures/rnaedit.p150byregion.pdf")
ggplot(maf.p150, aes(VAF, colour=p150)) + stat_ecdf() + theme_bw() + theme_classic()
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures/rnaedit.p150.edcf.png")
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures/rnaedit.p150.edcf.pdf")
# Only nonsynonymous and 3'UTR
ns.3utr <- subset(maf.1, Variant_Classification %in% c("3'UTR", "Nonsynonymous"))
ns.3utr <- merge(ns.3utr, comb.meta, by = "Sample")
ggplot(ns.3utr, aes(x=Variant_Classification, y=VAF, fill=Condition_code2)) + geom_boxplot()
ggviolin(ns.3utr, x="Condition_code2", y="VAF", fill = "Condition_code2", palette = "npg") + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.5)
ggviolin(ns.3utr, x="Condition_code2", y="VAF", fill = "Condition_code2", palette = "npg", xlab = "Phenotype", title = "VAF by Phenotype") + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.5)
ggboxplot(subset(ns.3utr, Cell_type.x=="Progenitor"), x="Condition", y="VAF", fill="Condition", title = "Progenitor Cells") + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + xlab(label = "Condition") + scale_fill_discrete(name="Condition")
ggboxplot(subset(ns.3utr, Cell_type.x=="Stem"), x="Condition", y="VAF", fill="Condition", title = "Stem Cells") + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + xlab(label = "Condition") + scale_fill_discrete(name="Condition")
ggboxplot(subset(ns.3utr, Cell_type.x=="Progenitor"), x="Treatment_type", y="VAF", fill="Treatment_type", title = "Progenitor Cells") + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + xlab(label = "Treatment") + scale_fill_discrete(name="Treatment")
ggboxplot(subset(ns.3utr, Cell_type.x=="Stem"), x="Treatment_type", y="VAF", fill="Treatment_type", title = "Stem Cells") + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + xlab(label = "Treatment") + scale_fill_discrete(name="Treatment")

# Statistical tests for differences in VAF between phenotypes: Kolmogorov Smirnov Test and T Test
ks.phenotype <- sapply(levels(prog.all$Condition_code2)[2:7], function(i) ks.test(subset(prog.all, Condition_code2=="Aged_Normal")$VAF, subset(prog.all, Condition_code2==i)$VAF, alternative = "greater")$p.value)
ttest.phenotype <- sapply(levels(prog.all$Condition_code2)[2:7], function(i) t.test(subset(prog.all, Condition_code2=="Aged_Normal")$VAF, subset(prog.all, Condition_code2==i)$VAF, alternative = "less")$p.value)
pv.phenotype.df <- data.frame(Kolmogorov_Smirnov=ks.phenotype, T_test=ttest.phenotype)
write.table(pv.phenotype.df, file="~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures/vaf_phenotype_statistical_signif.tsv", sep="\t", quote=FALSE, row.names=TRUE)

ks.risk <- sapply(levels(prog.all$Condition)[2:15], function(i) ks.test(subset(prog.all, Condition=="Normal")$VAF, subset(prog.all, Condition==i)$VAF, alternative = "greater")$p.value)
ttest.risk <- sapply(levels(prog.all$Condition)[2:15], function(i) t.test(subset(prog.all, Condition=="Normal")$VAF, subset(prog.all, Condition==i)$VAF, alternative = "less")$p.value)
pv.risk.df <- data.frame(Kolmogorov_Smirnov=ks.risk, T_test=ttest.risk)
write.table(pv.risk.df, file="~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures/vaf_risk_statistical_signif.tsv", sep="\t", quote=FALSE, row.names=TRUE)

ks.treatment <- sapply(levels(prog.all$Treatment_type)[2:7], function(i) ks.test(subset(prog.all, Treatment_type=="None")$VAF, subset(prog.all, Treatment_type==i)$VAF, alternative = "greater")$p.value)
ttest.treatment <- sapply(levels(prog.all$Treatment_type)[2:7], function(i) t.test(subset(prog.all, Treatment_type=="None")$VAF, subset(prog.all, Treatment_type==i)$VAF, alternative = "less")$p.value)
pv.treatment.df <- data.frame(Kolmogorov_Smirnov=ks.treatment, T_test=ttest.treatment)
write.table(pv.treatment.df, file="~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures/vaf_treatment_statistical_signif.tsv", sep="\t", quote=FALSE, row.names=TRUE)

# The below section generates extended data figure 4c (a table)
# Stratified by Variant Classification
ks.phenotype.by.vc <- sapply(levels(prog.all$Condition_code2)[2:7], function(i) lapply(unique(prog.all$Variant_Classification), function(j) ks.test(subset(prog.all, Variant_Classification==j & Condition_code2=="Aged_Normal")$VAF, subset(prog.all, Variant_Classification==j & Condition_code2==i)$VAF, alternative = "greater")$p.value))
row.names(ks.phenotype.by.vc) <- unique(prog.all$Variant_Classification)
write.table(ks.phenotype.by.vc, file="~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures/ks_phenotype_vc_statistical_signif.tsv", sep="\t", quote=FALSE, row.names=TRUE)


# Oncoplot
rkm <- read.maf(rnaedit.known.maf, clinicalData = rna.meta)
risk.colors <- RColorBrewer::brewer.pal(n = 12,name = "Paired")
names(risk.colors) <- c("AP_CML", "CML", "CP_CML", "CLL", "ET", "PV", "IR_MF", "HR_MF", "MF", "Aged_Normal", "Young_normal_bone_marrow", "AML")
risk.colors.sub <- risk.colors[names(risk.colors) %in% c("Aged_Normal", "Young_normal_bone_marrow", "ET", "PV", "MF", "AML", "CML")]
risk.order <- as.numeric(factor(names(risk.colors.sub), levels=c("Aged_Normal", "Young_normal_bone_marrow", "ET", "PV", "MF", "AML", "CML")))
risk.colors = list(Condition_code2 = risk.colors.sub)
dev.off(dev.list()["RStudioGD"])
# The below line generates figure 4d
oncoplot(maf = rkm,
         clinicalFeatures = "Condition_code2", sortByAnnotation = TRUE,
         annotationColor = risk.colors, top = 10,
         drawColBar = FALSE)

# The below section generates figure 4b and extended data figure 4a
# ADAR expression correlation
load("~/mnt/workstation/mnt/data1/adam/jamieson/holm/ADAR.tx.lcpm.exp.rdata")
mean.vaf <- do.call(rbind, lapply(unique(rnaedit.known.maf$Sample), function(i) data.frame(Sample=i, mean.vaf=mean(subset(rnaedit.known.maf, Sample==i)$VAF))))
adar.exp <- do.call(cbind, ADAR.tx.lcpm.exp)
adar.exp <- as.data.frame(t(as.matrix(adar.exp)))
adar.exp$Sample <- gsub("\\.", "-", gsub("X|_ACAGTG|_GTGAAA|_S[0-9]|.PB|.BM", "", row.names(adar.exp)))
adar.exp <- merge(adar.exp, mean.vaf, by="Sample")
adar.exp <- merge(adar.exp, rna.meta, by="Sample")
adar.exp <- merge(adar.exp, adar, by="Sample")
# Stem=1, circle Progenitor=2, triangle
adar.exp$Cell_type <- as.numeric(factor(adar.exp$Cell_type, levels=c("Stem", "Progenitor")))
adar.exp$Condition_code2 <- factor(adar.exp$Condition_code2, levels=c("Aged_Normal", "Young_Normal", "ET", "PV", "MF", "CML", "AML"))
# both stem and prog
ggscatter(adar.exp, x="ENST00000368474", y="mean.vaf", color="Condition_code2", size=c(4, 4)[adar.exp$Cell_type], shape=c(15, 17)[adar.exp$Cell_type]) + ylab("Mean RNA editing VAF") + xlab("ADAR p150 log CPM expression") + stat_cor(label.x = .5) + guides(colour = guide_legend(override.aes = list(size=4)))
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures/rnaedit.p150.meanvaf.corr.pdf")
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures/rnaedit.p150.meanvaf.corr.png")
ggscatter(adar.exp, x="ENST00000368471", y="mean.vaf", color="Condition_code2", size=c(4, 4)[adar.exp$Cell_type], shape=c(15, 17)[adar.exp$Cell_type]) + ylab("Mean RNA editing VAF") + xlab("ADAR p110 log CPM expression") + stat_cor(label.x = .5) + guides(colour = guide_legend(override.aes = list(size=4)))
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures/rnaedit.p110.meanvaf.corr.pdf")
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures/rnaedit.p110.meanvaf.corr.png")
ggplot(adar.exp, aes(p150, mean.vaf)) + geom_boxplot() + theme_classic() + ylab("Mean RNA editing VAF") + xlab("ADAR p150 expression quartile")
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures/rnaedit.p150.meanvaf.boxplot.pdf")
ggsave("~/mnt/workstation/mnt/data1/adam/jamieson/holm/pub_figures/rnaedit.p150.meanvaf.boxplot.png")
