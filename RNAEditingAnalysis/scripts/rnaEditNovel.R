### Novel RNA Edits from all 95 samples ###
library(data.table)
library(maftools)
library(tidyverse)

# Strand Data
strand <- unique(fread("~/mnt/workstation/mnt/data1/adam/software/annotation/gene_strand_hg19.txt", stringsAsFactors = FALSE, sep = "\t", header=FALSE))
names(strand) <- c("Hugo_Symbol", "Strand")

# Combined metadata
comb.meta <- read.csv("~/mnt/workstation/mnt/data1/adam/jamieson/holm/DNA_RNA_combined_meta06072018.csv", stringsAsFactors = FALSE, check.names = FALSE)
comb.meta$Sample <- gsub("\\.", "-", gsub("X|_GTGAAA|_ACAGTG|_S[0-9]|PB.|BM.", "", comb.meta$Sample))
comb.meta$Tumor_Sample_Barcode <- comb.meta$Sample
comb.meta$File_ID <- NULL
comb.meta$ID <- NULL
comb.meta <- unique(comb.meta)
comb.meta[grepl("AML", comb.meta$Sample_type),]$Sample_type <- "AML"
comb.meta$Sample_type <- factor(as.character(comb.meta$Sample_type),
                                levels = c("Aged normal bone marrow", "Young normal bone marrow", "MDS",
                                           "PV", "ET", "MF", "CML", "AML"))

# RNA Meta
rna.meta <- read.csv("~/mnt/workstation/mnt/data1/adam/jamieson/holm/fheditsHolm_Jamieson_RNAseq_with_controls_meta_20180326.csv", stringsAsFactors = FALSE, check.names = FALSE)
rna.meta$Sample <- gsub("\\.", "-", gsub("X|_GTGAAA|_ACAGTG|_S[0-9]|PB.|BM.", "", rna.meta$Sample))
rna.meta <- unique(rna.meta[c("Sample", "Condition_code", "Condition_code2", "Treatment_code2", "Cell_type", "Mutation", "Patient_ID")])
rna.meta[grepl("CML", rna.meta$Condition_code),]$Condition_code <- "CML"
rna.meta[grepl("AML", rna.meta$Condition_code),]$Condition_code <- "AML"
rna.meta[grepl("CML", rna.meta$Condition_code2),]$Condition_code2 <- "CML"
rna.meta[grepl("AML", rna.meta$Condition_code2),]$Condition_code2 <- "AML"
rna.meta$Tumor_Sample_Barcode <- sapply(strsplit(rna.meta$Sample, "-"), `[`, 1)

names <- c("Hugo_Symbol", "Chromosome", "Start_position", "End_position", "Variant_Classification", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "dbSNP_RS", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", "Protein_Change", "t_alt_count", "t_ref_count", "1000gp3_AF", "ExAC_AF", "dbNSFP_CADD_phred")

# Novel MPN
# The file "rnaedit.novel.alu.merged.oncotator.filt.maf" is equivalent to the file "disease_alu_merged_novel_rnaedit_sites.oncotator.filt.maf" produced in the RNA Editing Analysis methods
mpn.novel.alu <- fread("~/mnt/jam_db/projects/holm/vcf/rnaedit/rnaedit.novel.alu.merged.oncotator.filt.maf", stringsAsFactors = FALSE, sep = "\t")
names(mpn.novel.alu) <- names
mpn.novel.alu$Novelty <- "Novel Alu"
# The file "rnaedit.novel.nonalu.merged.oncotator.filt.maf" is equivalent to the file "disease_nonalu_merged_novel_rnaedit_sites.oncotator.filt.maf" produced in the RNA Editing Analysis methods
mpn.novel.nonalu <- fread("~/mnt/jam_db/projects/holm/vcf/rnaedit/rnaedit.novel.nonalu.merged.oncotator.filt.maf", stringsAsFactors = FALSE, sep = "\t")
names(mpn.novel.nonalu) <- names
mpn.novel.nonalu$Novelty <- "Novel NonAlu"

# Novel Normal
# The file "rnaedit.novel.normal.alu.merged.oncotator.filt.maf" is equivalent to the file "normal_alu_merged_novel_rnaedit_sites.oncotator.filt.maf" produced in the RNA Editing Analysis methods
normal.novel.alu <- fread("~/mnt/jam_db/projects/holm/vcf/rnaedit/novel/rnaedit.novel.normal.alu.merged.oncotator.filt.maf", stringsAsFactors = FALSE, sep = "\t")
names(normal.novel.alu) <- names
normal.novel.alu$Novelty <- "Novel Alu"
# The file "rnaedit.novel.normal.nonalu.merged.oncotator.filt.maf" is equivalent to the file "normal_nonalu_merged_novel_rnaedit_sites.oncotator.filt.maf" produced in the RNA Editing Analysis methods
normal.novel.nonalu <- fread("~/mnt/jam_db/projects/holm/vcf/rnaedit/novel/rnaedit.novel.normal.nonalu.merged.oncotator.filt.maf", stringsAsFactors = FALSE, sep = "\t")
names(normal.novel.nonalu) <- names
normal.novel.nonalu$Novelty <- "Novel NonAlu"

# Combine
maf <- do.call(rbind, list(mpn.novel.alu, mpn.novel.nonalu, normal.novel.alu, normal.novel.nonalu))
maf$Tumor_Sample_Barcode <- gsub("-BM|-PB|_S[0-9]|_ACAGTG|_GTGAAA", "", maf$Matched_Norm_Sample_Barcode)
maf$Sample <- maf$Tumor_Sample_Barcode
maf$t_alt_count <- as.numeric(maf$t_alt_count)
maf$t_ref_count <- as.numeric(maf$t_ref_count)
maf$VAF <- maf$t_alt_count/(maf$t_alt_count + maf$t_ref_count)
maf2 <- subset(maf, t_alt_count >= 2 & (t_alt_count + t_ref_count >= 5) & VAF >= 0.10)
maf3 <- merge(maf2, strand, by = "Hugo_Symbol")
maf3$mutsig <- paste0(maf3$Reference_Allele, ">", maf3$Tumor_Seq_Allele2)
maf4 <- subset(maf3, (mutsig == "A>G" & Strand=="+") | (mutsig == "T>C" & Strand=="-"))
maf5 <- subset(maf4, !(Tumor_Sample_Barcode %in% c(91,96,97)))
maf6 <- unique(maf5)
maf7 <- subset(maf6, !(Tumor_Sample_Barcode %in% c("586sp1", "22_10", "22_13", "22_14", "22_15", "22_16", "4689sp1")))
maf7$VAF <- maf7$t_alt_count/(maf7$t_alt_count + maf7$t_ref_count)
stroma.maf <- subset(maf6, Tumor_Sample_Barcode %in% c("586sp1", "22_10", "22_13", "22_14", "22_15", "22_16", "4689sp1"))

# oncoplot
# The below section generates figure 4c
novel <- read.maf(maf7, clinicalData = comb.meta)
Sample_type.colors <- RColorBrewer::brewer.pal(n = 7, name = 'Paired')
names(Sample_type.colors) <- c("Aged_normal_bone_marrow", "Young_normal_bone_marrow", "MF", "PV", "ET", "CML", "AML")
Treatment_type.colors <- RColorBrewer::brewer.pal(n = 7, name = 'Accent')
names(Treatment_type.colors) <- c("None", "Untreated", "hydroxyurea", "Jak2_Inhibitor", "TKI", "vidaza", "SHH_treated")
cell_type.colors <- gray.colors(2)
names(cell_type.colors) <- c("Stem", "Progenitor")
colors <- list(Sample_type = Sample_type.colors, Treatment_type=Treatment_type.colors, Cell_type=cell_type.colors)
dev.off(dev.list()["RStudioGD"])
oncoplot(maf = novel, clinicalFeatures = c("Sample_type", "Treatment_type", "Cell_type"), sortByAnnotation = TRUE,
        showTumorSampleBarcodes = TRUE, annotationColor = colors, top = 25, fontSize = 6)

# stroma
novel.stroma <- read.maf(stroma.maf, clinicalData = comb.meta)
Sample_type.colors <- RColorBrewer::brewer.pal(n = 7, name = 'Paired')
names(Sample_type.colors) <- c("Aged_normal_bone_marrow", "Young_normal_bone_marrow", "MF", "PV", "ET", "CML", "AML")
Treatment_type.colors <- RColorBrewer::brewer.pal(n = 7, name = 'Accent')
names(Treatment_type.colors) <- c("None", "Untreated", "hydroxyurea", "Jak2_Inhibitor", "TKI", "vidaza", "SHH_treated")
colors <- list(Sample_type = Sample_type.colors, Treatment_type=Treatment_type.colors)
dev.off(dev.list()["RStudioGD"])
oncoplot(maf = novel.stroma, clinicalFeatures = c("Sample_type", "Treatment_type"), sortByAnnotation = TRUE,
         showTumorSampleBarcodes = TRUE, annotationColor = colors)
