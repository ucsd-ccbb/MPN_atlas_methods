# Filter SV annotated with AnnotSV
library(data.table)

raw.sv <- fread("~/mnt/workstation/mnt/data1/adam/jamieson/holm/sv/filtered/mpn/AnnotSV.MPNandControls.tsv", stringsAsFactors = FALSE, sep = "\t")
sv <- raw.sv
sv <- subset(sv, `AnnotSV ranking` == 5)
sv <- sv[sv$`1000g_max_AF` < 0.01,]
sv <- sv[sv$GD_POPMAX_AF < 0.01,]
sv <- sv[sv$`SV chrom` != "X"]

###################################################################################################################################
# Circos Plot Prep
# remove normal SV
normals <- c("743-PB", "780-PB", "792-PB", "795-PB")
normals <- as.matrix(sv)[,normals]
sv <- sv[rowSums(apply(normals, 2, function(i) i == "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN")) == 4,]
sv.mat <- as.matrix(sv)[,names(sv)[grepl("-PB", names(sv))]]
sv.mat <- sv.mat[,!colnames(sv.mat) %in% c("780-PB", "792-PB", "795-PB")]
samples <- colnames(sv.mat)
sample.mats <- lapply(samples, function(i) cbind(sv[,1:13], Gene=sv$`Gene name`, GT=sv[[i]], Sample=i))
sample.mats <- lapply(sample.mats, function(i) subset(i, GT != "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN"))
sv.for.circos <- unique(do.call(rbind, sample.mats)[,c("SV chrom", "SV start", "SV end", "SV type", "Sample")])
names(sv.for.circos) <- c("chromosome", "start", "end", "SVtype", "Sample")
sv.for.circos$Class <- "SV"

dna.meta <- read.csv("~/mnt/workstation/mnt/data1/adam/jamieson/holm/metadata/somatic_metadata_clean_12062018.txt", sep = "\t", stringsAsFactors = FALSE)
et.samples <- paste0(subset(dna.meta, Diagnosis == "ET")$Tumor_Sample_Barcode, "-PB")
pv.samples <- paste0(subset(dna.meta, Diagnosis == "PV")$Tumor_Sample_Barcode, "-PB")
mf.samples <- paste0(subset(dna.meta, Diagnosis == "MF")$Tumor_Sample_Barcode, "-PB")

load("~/mnt/workstation/mnt/data1/adam/jamieson/holm/cnvkit/pooled_normal/noLCR/CNVkit_circos_data.RData")
et.sv <- subset(sv.for.circos, Sample %in% et.samples)
et.sv <- rbind(et.cnv, et.sv)
pv.sv <- subset(sv.for.circos, Sample %in% pv.samples)
pv.sv <- rbind(pv.cnv, pv.sv)
mf.sv <- subset(sv.for.circos, Sample %in% mf.samples)
mf.sv <- rbind(mf.cnv, mf.sv)
save(et.sv, pv.sv, mf.sv, file="~/mnt/workstation/mnt/data1/adam/jamieson/holm/sv/filtered/mpn/SV_circos_data.RData")

###################################################################################################################################
# create sv.mat for top genes
sv.type <- sv$`SV type`
sv.mat <- apply(sv.mat, 2, function(i) ifelse(i == "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN", yes=NA, no=sv.type))
sv.mat[sv.mat=="DEL"] <- "SV_deletion"
sv.mat[sv.mat=="DUP"] <- "SV_duplication"
samples <- colnames(sv.mat)
sv.mat <- cbind(gene=sv$`Gene name`, sv.mat)

svList <- function(sample, sv.matrix){
  df <- sv.matrix[,c("gene", sample)]
  df <- data.frame(unique(subset(df, !is.na(df[,sample]))), stringsAsFactors=FALSE)
  names(df) <- c("gene", "SVTYPE")
  df
}

(sv.list <- lapply(samples, function(i) svList(i, sv.mat)))

sv.mat.exp <- matrix(ncol=length(samples), nrow=length(pc_genes))
colnames(sv.mat.exp) <- samples
row.names(sv.mat.exp) <- pc_genes
for(i in 1:length(samples)){sv.mat.exp[,i] <- sv.list[[i]][match(row.names(sv.mat.exp), sv.list[[i]]$gene), "SVTYPE"]}
save(sv.mat.exp, file="~/mnt/workstation/mnt/data1/adam/jamieson/holm/sv/filtered/mpn/MantaLumpy_merged_genotypes.pass.AnnotSV.sv_mat.RData")



