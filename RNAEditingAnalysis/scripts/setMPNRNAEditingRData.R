# Set up .RData file to load into Shiny app quickly
setwd("/mnt/data1/adam/jamieson/holm")
library(data.table)
library(plyr)

# MPN
mpn.meta <- read.csv("metadata/RNAseq_with_controls_meta_20200304.txt", stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
mpn.meta$Treatment_type <- factor(as.character(mpn.meta$Treatment_type),
                                  levels=c("None", "Untreated", "Jak2 Inhibitor", "SHH treated", "hydroxyurea", "TKI", "vidaza"))
mpn.meta <- subset(mpn.meta, Cell.type %in% c("Stem", "Progenitor"))
condition.levels <- c("Aged_Normal", "Young_Normal", "ET", "PV", "MF", "CML", "sAML", "denovoAML")
mpn.meta$Condition <- factor(mpn.meta$Condition_code2, levels=condition.levels)

mpn_maf_path <- "vcf/rnaedit/MPN_Normal_known_novel_merged_rnaedit_sites_03042020.maf"
mpn <- fread(mpn_maf_path, sep = "\t")
mpn$Condition <- mpn$Condition_code2
mpn <- subset(mpn, Condition %in% condition.levels & Cell.type %in% c("Progenitor", "Stem"))
mpn$Condition <- factor(mpn$Condition, levels=condition.levels)
mpn$Cell_type <- factor(mpn$Cell.type, levels=c("Progenitor", "Stem"))
mpn$VarID <- paste0(mpn$Chromosome, ":", mpn$Start_position)

save(mpn, mpn.meta, file="rdata/MPN_MAF_files.RData")
