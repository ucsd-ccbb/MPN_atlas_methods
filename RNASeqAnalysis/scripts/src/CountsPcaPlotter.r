# ---------------------------------------------------------------------------------
# Copyright (c) 2018 UC San Diego Center for Computational Biology & Bioinformatics
#
# Distributed under the terms of the MIT License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------------
# Initial author: Amanda Birmingham

library(edgeR)

expandDesignDf<-function(countsDf, designDf,
                         sampleNameColName = "sample_name"){
  aDgeList <- DGEList(counts=countsDf)
  mergedDesignDf = merge(x=designDf, y=aDgeList$samples,
                         by.y="row.names",
                         by.x=sampleNameColName)
  return(mergedDesignDf)
}

syncCountSampleOrderToDesignDf<-function(counts_df, designDf,
                                         sampleColName="sample_name"){

  sampleNamesInOrder = designDf[[sampleColName]]

  # check for samples in the design file that aren't in the counts file
  missingSamples = setdiff(sampleNamesInOrder, colnames(counts_df))
  if (length(missingSamples)>0){
    print(missingSamples)
    stop("Above samples are in design file but missing from counts file")
  }

  # ensure that the order of the samples in the counts table is
  # the same as the order of the samples in the design table
  reordered_counts_df = counts_df[sampleNamesInOrder]
  return(reordered_counts_df)
}

reformatDfForPca<-function(counts_df, designDf,
                           sampleColName="sample_name"){
  reordered_counts_df = syncCountSampleOrderToDesignDf(
    counts_df, designDf, sampleColName)

  # now transform the counts df so it is samples are in rows
  # (as in the design file) rather than in rows
  transformed_df = t(reordered_counts_df)
  return(transformed_df)
}

makeAndPrintRawCountsPca<-function(countsDf, designDf,
                                   pointShapeColName,
                                   designSampleNameColName="sample_name",
                                   libSizeColName = "lib.size",
                                   designColNameForLabels=NULL, labelOnlyOutliers=TRUE){

  if (!libSizeColName %in% colnames(designDf)){
    designDf = expandDesignDf(countsDf, designDf,
                              designSampleNameColName)
  }

  rawTitle = "PCA of Raw Counts"
  display_markdown(rawTitle)
  countsPca = doPcaFromSamplesAsColsDf(countsDf, designDf,
                                       designSampleNameColName)
  rawPlot = make2dPcaPlot(countsPca, designDf, pointShapeColName,
                          libSizeColName, designColNameForLabels,
                          labelOnlyOutliers)
  print(rawPlot + ggtitle(rawTitle))
}

makeAndPrintCpmsPca<-function(countsDf, designDf,
                              pointShapeColName,
                              designSampleNameColName="sample_name",
                              libSizeColName = "lib.size",
                              designColNameForLabels=NULL, labelOnlyOutliers=TRUE){

  if (!libSizeColName %in% colnames(designDf)){
    designDf = expandDesignDf(countsDf, designDf,
                              designSampleNameColName)
  }

  cpmsDf = getCpmsDf(countsDf)

  normTitle = "PCA of Normalized Counts"
  designAndPca = makeAndPrintPca(normTitle)

  display_markdown(normTitle)
  cpmsPca = doPcaFromSamplesAsColsDf(cpmsDf, designDf,
                                     designSampleNameColName)
  normPlot = make2dPcaPlot(cpmsPca, designDf, pointShapeColName,
                           libSizeColName, designColNameForLabels,
                           labelOnlyOutliers)
  print(normPlot + ggtitle(normTitle))

  return(list(rawcounts=countsDf, design=designDf, cpms=cpmsDf,
              cpmsPca=cpmsPca))
}

# TODO: come back and integrate this function with
# makeAndPrintRawCountsPca and makeAndPrintCpmsPca
makeAndPrintPca<-function(title, countsDf, designDf,
                          pointShapeColName,
                          designSampleNameColName="sample_name",
                          libSizeColName = "lib.size",
                          designColNameForLabels=NULL, labelOnlyOutliers=TRUE){

  if (!libSizeColName %in% colnames(designDf)){
    designDf = expandDesignDf(countsDf, designDf,
                              designSampleNameColName)
  }

  display_markdown(title)
  aPca = doPcaFromSamplesAsColsDf(countsDf, designDf, designSampleNameColName)
  aPlot = make2dPcaPlot(aPca, designDf, pointShapeColName,
                           libSizeColName, designColNameForLabels,
                           labelOnlyOutliers)
  print(aPlot + ggtitle(title))

  return(list(design=designDf, pca=aPca))
}

getCpmsDf<-function(counts_df){
  y <- DGEList(counts=counts_df)
  cpm_matrix = cpm(y)
  cpmDf = data.frame(cpm_matrix)
  return(cpmDf)
}
