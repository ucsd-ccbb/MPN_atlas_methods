# Filter SV annotated with AnnotSV
library(data.table)

raw.sv <- fread("~/mnt/workstation/mnt/data1/adam/jamieson/holm/sv/merged/AnnotSV.zSpKPU1NZX.tsv", stringsAsFactors = FALSE, sep = "\t")
sv <- raw.sv
sv <- sv[sv$`1000g_max_AF` < 0.01,]
sv <- sv[sv$GD_POPMAX_AF < 0.01,]
sv <- subset(sv, dbVar_status != "")
sv <- subset(sv, !grepl("enign", dbVar_status))
sv <- sv[sv$pLI_ExAC > 0.9,] # Score computed by ExAC indicating the probability that a gene is intolerant to a loss of function variation (Nonsense, splice acceptor/donor variants due to SNV/indel). ExAC consider pLI>=0.9 as an extremely LoF intolerant gene
sv <- subset(sv, `AnnotSV ranking` >= 4)

# remove normal SV
normals <- c("780-PB", "792-PB", "795-PB")
normals <- as.matrix(sv)[,normals]
sv <- sv[rowSums(apply(normals, 2, function(i) i == "./.:.:.:.:.:.:.:.:.:.")) != 0,]
write.table(sv, file="~/mnt/workstation/mnt/data1/adam/jamieson/holm/sv/merged/SV2_merged_genotypes.pass.AnnotSV.tsv", row.names = FALSE, sep = "\t", quote = FALSE)

# create sv.mat
sv.mat <- as.matrix(sv)[,names(sv)[grepl("-PB", names(sv))]]
sv.type <- sv$`SV type`
sv.mat <- apply(sv.mat, 2, function(i) ifelse(i == "./.:.:.:.:.:.:.:.:.:.", yes=NA, no=sv.type))
sv.mat[sv.mat=="DEL"] <- "SV_deletion"
sv.mat[sv.mat=="DUP"] <- "SV_duplication"
sv.mat <- sv.mat[,!colnames(sv.mat) %in% c("780-PB", "792-PB", "795-PB")]
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

save(sv.mat.exp, file="~/mnt/workstation/mnt/data1/adam/jamieson/holm/sv/merged/SV2_merged_genotypes.pass.AnnotSV.sv_mat.RData")

# 4 43747338 81614087