library(tidyverse)
library(KEGGREST)
library(clusterProfiler)
library(ReactomePA)
library(mygene)

findDiffEditedGene <- function(gene, maf, cond1, cond2, n.cond1, n.cond2){
  dis.gene <- subset(maf, Hugo_Symbol == gene)

  n.possible.sites <- uniqueN(dis.gene$Start_position)
  n.possible.sites.cond1 <- n.cond1*n.possible.sites
  n.possible.sites.cond2 <- n.cond2*n.possible.sites

  n.cond1.gene <- nrow(subset(dis.gene, Condition_code2 == cond1))
  n.cond2.gene <- nrow(subset(dis.gene, Condition_code2 == cond2))

  mat <- rbind(Edited=c(n.cond1.gene, n.cond2.gene), notEdited=c((n.possible.sites.cond1-n.cond1.gene), (n.possible.sites.cond2 - n.cond2.gene)))
  Xsq <- chisq.test(mat)
  data.frame(p.value=Xsq$p.value, chi_square_statistic=Xsq$statistic,
             Gene=gene,
             n_cond1=n.cond1.gene, n_cond2=n.cond2.gene,
             n_samples_cond1=n.cond1, n_samples_cond2=n.cond2,
             n_possible_sites_cond1=n.possible.sites.cond1, n_possible_sites_cond2=n.possible.sites.cond2)
}

pAdjustWithUniverse <- function(p.value.df, gene_universe){
  n.other.genes <- length(gene_universe) - nrow(p.value.df)
  universe.df <- data.frame(p.value=integer(n.other.genes) + 1, chi_square_statistic=integer(n.other.genes) + 1,
                            Gene=NA,
                            n_cond1=integer(n.other.genes), n_cond2=integer(n.other.genes),
                            n_samples_cond1=integer(n.other.genes), n_samples_cond2=integer(n.other.genes),
                            n_possible_sites_cond1=integer(n.other.genes), n_possible_sites_cond2=integer(n.other.genes))
  sig.df <- rbind(p.value.df, universe.df)
  #sig.df.adj <- sig.df[p.adjust(sig.df$p.value, method = "BH") <= 0.05,]
  sig.df$p.adjust <- p.adjust(sig.df$p.value, method = "BH")
  sig.df <- sig.df[c("p.value", "p.adjust", "chi_square_statistic", "Gene", "n_cond1", "n_cond2", "n_samples_cond1", "n_samples_cond2", "n_possible_sites_cond1", "n_possible_sites_cond2")]
  arrange(sig.df, p.value)
}

# read in the known rnaedit maf produced by rnaedit.known.comparison.R
rnaedit.known.maf <- fread("~/mnt/jam_db/projects/holm/vcf/rnaedit/RNAedit_known_Normal_MPN_merged_filtered_09172018.maf", stringsAsFactors = FALSE, sep="\t")

## MF vs AN
prog.mf <- subset(rnaedit.known.maf, Cell_type=="Progenitor" & (Condition_code2 == "MF" | Condition_code2 == "Aged_Normal") & Matched_Norm_Sample_Barcode != "380")
prog.mf <- unique(prog.mf)
expressed.genes <- read.csv("~/Downloads/Limma_DE_AgedBMvsMF_Prog_Holm_Jamieson_RNASeq_TCW.csv", stringsAsFactors = FALSE)
universe.for.padj <- unique(prog.mf$Hugo_Symbol)

n.normal <- uniqueN(subset(prog.mf, Condition_code2 == "Aged_Normal")$Tumor_Sample_Barcode)
n.mpn <- uniqueN(subset(prog.mf, Condition_code2 != "Aged_Normal")$Tumor_Sample_Barcode)
diffedited.prog.mf.byve <- lapply(unique(prog.mf$Variant_Classification),
                                  function(j) do.call(rbind, lapply(unique(subset(prog.mf, Variant_Classification==j)$Hugo_Symbol),
                                                                    function(i) findDiffEditedGene(i, subset(prog.mf, Variant_Classification==j), "Aged_Normal", "MF", n.normal, n.mpn))))
names(diffedited.prog.mf.byve) <- unique(prog.mf$Variant_Classification)
diffedited.prog.mf.byve.adj <- lapply(diffedited.prog.mf.byve, function(i) pAdjustWithUniverse(i, universe.for.padj))
diffedited.prog.mf.byve.adj.df <- unique(do.call(rbind, diffedited.prog.mf.byve.adj))

# AML vs AN
prog.aml <- subset(rnaedit.known.maf, Cell_type=="Progenitor" & (Condition_code2 == "AML" | Condition_code2 == "Aged_Normal") & Matched_Norm_Sample_Barcode != "380")
prog.aml <- unique(prog.aml)
expressed.genes.aml <- read.csv("~/Downloads/Limma_DE_AgedBMvsAML_Prog_Holm_Jamieson_RNASeq_TCW.csv", stringsAsFactors = FALSE)
universe.for.padj.amlvsan <- unique(prog.aml$Hugo_Symbol)
n.normal <- uniqueN(subset(prog.aml, Condition_code2 == "Aged_Normal")$Tumor_Sample_Barcode)
n.aml <- uniqueN(subset(prog.aml, Condition_code2 = "AML")$Tumor_Sample_Barcode)
diffedited.prog.aml.byve <- lapply(unique(prog.aml$Variant_Classification),
                                  function(j) do.call(rbind, lapply(unique(subset(prog.aml, Variant_Classification==j)$Hugo_Symbol),
                                                                    function(i) findDiffEditedGene(i, subset(prog.aml, Variant_Classification==j), "Aged_Normal", "AML", n.normal, n.aml))))
diffedited.prog.aml.byve.adj <- lapply(diffedited.prog.aml.byve, function(i) pAdjustWithUniverse(i, universe.for.padj.amlvsan))
diffedited.prog.aml.byve.adj.df <- unique(do.call(rbind, diffedited.prog.aml.byve.adj))

# AML vs MF
prog.amlvsmf <- subset(rnaedit.known.maf, Cell_type=="Progenitor" & (Condition_code2 == "AML" | Condition_code2 == "MF") & Matched_Norm_Sample_Barcode != "380")
prog.amlvsmf <- unique(prog.amlvsmf)
expressed.genes.amlvsmf <- read.csv("~/Downloads/Limma_DE_AMLvsMF_Prog_Holm_Jamieson_RNASeq_TCW.csv", stringsAsFactors = FALSE)
universe.for.padj.amlvsmf <- unique(prog.amlvsmf$Hugo_Symbol)
n.mf <- uniqueN(subset(prog.amlvsmf, Condition_code2 == "MF")$Tumor_Sample_Barcode)
n.aml <- uniqueN(subset(prog.amlvsmf, Condition_code2 == "AML")$Tumor_Sample_Barcode)
diffedited.prog.amlvsmf.byve <- lapply(unique(prog.amlvsmf$Variant_Classification),
                                   function(j) do.call(rbind, lapply(unique(subset(prog.amlvsmf, Variant_Classification==j)$Hugo_Symbol),
                                                                     function(i) findDiffEditedGene(i, subset(prog.amlvsmf, Variant_Classification==j), "MF", "AML", n.mf, n.aml))))
diffedited.prog.amlvsmf.byve.adj <- lapply(diffedited.prog.amlvsmf.byve, function(i) pAdjustWithUniverse(i, universe.for.padj.amlvsmf))
diffedited.prog.amlvsmf.byve.adj.df <- unique(do.call(rbind, diffedited.prog.amlvsmf.byve.adj))


# Write to tables
# These output files were used as inputs to the creation of figure 4f, figure 4g, extended data figure 4g, and extended data figure 4i
write.table(diffedited.prog.mf.byve.adj.df, file="~/mnt/workstation/mnt/data1/adam/jamieson/holm/differentially_edited_genes_MFvsAN_progenitors09102018.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(diffedited.prog.aml.byve.adj.df, file="~/mnt/workstation/mnt/data1/adam/jamieson/holm/differentially_edited_genes_AMLvsAN_progenitors09102018.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(diffedited.prog.amlvsmf.byve.adj.df, file="~/mnt/workstation/mnt/data1/adam/jamieson/holm/differentially_edited_genes_AMLvsMF_progenitors09102018.txt", sep="\t", quote=FALSE, row.names=FALSE)
