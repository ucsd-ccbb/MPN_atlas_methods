# RNA Editing Functions
library(shiny)
library(data.table)
library(DT)
library(ggplot2)
library(ggpubr)
library(maftools)
library(Cairo)
options(shiny.maxRequestSize=1000*1024^2)
options(shiny.usecairo=T)

ggViolin <- function(data, column){
  ggviolin(data, x="Condition", y="VAF", fill = "Condition", palette = "Set2") + 
    stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.5)+ theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

ggBox <- function(data){
  ggplot(data, aes(x=Variant_Classification, y=VAF, fill=Condition)) + geom_boxplot() + theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_brewer(palette="Set2")
}

topGenesMPN <- function(data, metadata){
  Condition.colors <- RColorBrewer::brewer.pal(n = 8, name = 'Set2')
  names(Condition.colors) <- c("Aged_Normal", "Young_Normal", "ET", "PV", "MF", "CML", "sAML", "denovoAML")
  colors <- list(Condition = Condition.colors)
  maf.in <- read.maf(maf, clinicalData = metadata)
  oncoplot(maf = maf.in, clinicalFeatures = c("Condition"), sortByAnnotation = TRUE, 
           showTumorSampleBarcodes = TRUE, annotationColor = colors, top = 25, fontSize = 0.6,
           annotationOrder = levels(data),
           bgCol="white", removeNonMutated = FALSE, barcode_mar = 6, gene_mar = 6, fill = TRUE)
}

topGenes <- function(data, metadata){
  maf.in <- read.maf(data, clinicalData = metadata)
  oncoplot(maf.in, clinicalFeatures = "Condition", showTumorSampleBarcodes = TRUE, sortByAnnotation = TRUE)
}

editsPerGene <- function(data){
  EditsPerGene <- do.call(rbind, lapply(unique(data$Condition), function(i)
    cbind(Condition=i, data.frame(table(subset(data, Condition==i)$Hugo_Symbol)))))
  EditsPerGene$Condition <- factor(EditsPerGene$Condition, levels=levels(data$Condition))
  ggplot(EditsPerGene, aes(x=log(Freq), fill=Condition)) + geom_density(alpha=0.7) + theme_bw() + ylab("Log Edits Per Gene") + scale_fill_brewer(palette="Set2")
}

editedGenesPerSampleDensity <- function(data){
  EditedGenesPerSample <- do.call(rbind, lapply(unique(data$Condition), function(i)
    cbind(Condition=i, aggregate(Hugo_Symbol ~ Tumor_Sample_Barcode, data=subset(data, Condition==i), length))))
  EditedGenesPerSample$Condition <- factor(EditedGenesPerSample$Condition, levels=levels(data$Condition))
  ggplot(EditedGenesPerSample, aes(x=Hugo_Symbol, fill=Condition)) + geom_density(alpha=0.7) + theme_bw() + xlab("Density") + ylab("Edited Genes Per Sample") + scale_fill_brewer(palette="Set2")
}

editedGenesPerSampleBox <- function(data){
  EditedGenesPerSample <- do.call(rbind, lapply(unique(data$Condition), function(i)
    cbind(Condition=i, aggregate(Hugo_Symbol ~ Tumor_Sample_Barcode, data=subset(data, Condition==i), length))))
  EditedGenesPerSample$Condition <- factor(EditedGenesPerSample$Condition, levels=levels(data$Condition))
  ggplot(EditedGenesPerSample, aes(x=Condition, y=Hugo_Symbol, fill=Condition)) + geom_boxplot() + theme_bw() + ylab("Edited Genes Per Sample") + scale_fill_brewer(palette="Set2")
}


