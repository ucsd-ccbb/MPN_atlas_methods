# ---------------------------------------------------------------------------------
# Copyright (c) 2018 UC San Diego Center for Computational Biology & Bioinformatics
#
# Distributed under the terms of the MIT License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------------
# Initial author: Amanda Birmingham

library(cowplot)
library(ggplot2)
library(grid)
library(IRdisplay)

makeAndPrintPca <-function(data_df, design_df,
                           design_col_name_for_shapes=NULL, design_col_name_for_colors=NULL,
                           design_col_name_for_labels=NULL, labelOutliersOnly=TRUE,
                           shrink_viewport=FALSE) {

  pcaResults = doAndPrintScaledPcaOnSamplesAsRowsDf(data_df)
  makeAndPrintPcaPlot(pcaResults, design_df, design_col_name_for_shapes,
                      design_col_name_for_colors, design_col_name_for_labels,
                      labelOutliersOnly, TRUE, shrink_viewport)
}

makeAndPrintPcaPlot<-function(pcaResults, design_df=NULL,
                              design_col_name_for_shapes=NULL, design_col_name_for_colors=NULL,
                              design_col_name_for_labels=NULL, labelOutliersOnly=TRUE,
                              add_hotelling_ellipse=TRUE, shrink_viewport=FALSE) {

  pcaPlot = make2dPcaPlot(pcaResults, design_df,
                          design_col_name_for_shapes, design_col_name_for_colors,
                          design_col_name_for_labels, labelOutliersOnly,
                          add_hotelling_ellipse)
  printPlotInViewport(pcaPlot, shrink_viewport)
}

printPlotInViewport<-function(pcaPlot, shrink_viewport=FALSE){
  # NB that this method does NOT set the canvas back to the
  # default size after being called--when I try to do that,
  # the reset happens before the plot is rendered, thus
  # nullifying my attempts to resize the canvas to fit the
  # image (even if I try using Sys.sleep, etc).
  # Until I can spend more time exploring how to prevent that,
  # it is necessary to call resetPlotSize() after any run
  # of this method.

  viewport_val = NULL
  if (shrink_viewport==TRUE) {
    viewport_val = viewport(width=unit(0.8, "npc"))
  }

  startingWidth = getOption("repr.plot.width")
  startingHeight = getOption("repr.plot.height")

  # resize image to max width, appropriate height
  # to remove excessive whitespace in default square
  # image canvas if real image is not square
  aspectRatio = findPlotAspectRatio(pcaPlot)
  plotHeight = startingWidth/aspectRatio
  options(repr.plot.width=startingWidth, repr.plot.height=plotHeight)

  suppressWarnings(print(pcaPlot, vp=viewport_val))
}

make2dPcaPlot<-function(pca_result, design_df=NULL,
                        design_col_name_for_shapes=NULL, design_col_name_for_colors=NULL,
                        design_col_name_for_labels=NULL, label_only_outliers=TRUE,
                        add_hotelling_ellipse=TRUE){

  scores = data.frame(pca_result$x[,1:2])

  if (add_hotelling_ellipse){
    hotelling_ellipse = data.frame(getHotellingT2Ellipse(
      pca_result$x[,1], pca_result$x[,2]))
    colnames(hotelling_ellipse) = c("PC1", "PC2")
  }

  shape_values = rep("",nrow(pca_result$x))
  color_values = rep("",nrow(pca_result$x))
  label_values = rep("",nrow(pca_result$x))
  if (!is.null(design_col_name_for_shapes)){
    design_df[[design_col_name_for_shapes]] = factor(
      design_df[[design_col_name_for_shapes]])
    shape_values = design_df[[design_col_name_for_shapes]]
    scores = cbind(scores, shape_values)
  }
  if (!is.null(design_col_name_for_colors)){
    color_values = design_df[[design_col_name_for_colors]]
    scores = cbind(scores, color_values)
  }
  if (!is.null(design_col_name_for_labels)){
    if (label_only_outliers) {
      includeValues = getWhetherPointsAreOutliers(pca_result)
      label_values = ifelse(includeValues,
                            design_df[[design_col_name_for_labels]],'')
    } else {
      label_values = design_df[[design_col_name_for_labels]]
    }
    scores = cbind(scores, label_values)
  }

  pc1.2 = ggplot(scores, aes(x=PC1, y=PC2)) +
    geom_point(aes(shape=shape_values,
                   color=color_values), size = 4) + 
    scale_shape_manual(values=c(0:length(shape_values))) +
    coord_fixed(1/1) +
    labs(color=design_col_name_for_colors,
         shape=design_col_name_for_shapes)

  if (add_hotelling_ellipse){
    pc1.2 = pc1.2 + geom_path(data=hotelling_ellipse)
  }

  if (!is.null(design_col_name_for_labels)){
    pc1.2 = pc1.2 + geom_text(aes(label=label_values),
                              hjust=0, vjust=0)
  }

  if (is.numeric(color_values)) {
    pc1.2 = pc1.2 + scale_color_gradient(low="blue", high="red")
  }

  pc1.2 = pc1.2 + coord_fixed()
  return (pc1.2)
}

findPlotAspectRatio<-function(aGgplot){
  # get the x- and y-axis ranges actually used in the graph
  builtPlot = ggplot_build(aGgplot)

  # pre-ggplot2 version 2.2
  yRange <- builtPlot$panel$ranges[[1]]$y.range
  xRange <- builtPlot$panel$ranges[[1]]$x.range

  # ggplot2 version 2.2 and later
  if (is.null(yRange)){
    yRange = builtPlot$layout$panel_ranges[[1]]$y.range
    xRange <- builtPlot$layout$panel_ranges[[1]]$x.range
  }

  aspectRatio <- (max(xRange)-min(xRange))/(max(yRange)-min(yRange))
  return(aspectRatio)
}

doAndPrintScaledPcaOnSamplesAsRowsDf<-function(data_df){
  pcaResults = doScaledPcaOnSamplesAsRowsDf(data_df)
  display(summary(pcaResults)$importance)
  return(pcaResults)
}

doPcaFromSamplesAsColsDf<-function(samplesAsColsDf, designDf,
                                   sampleNameDesignColName = "sample_name"){

  transformedDf = reformatDfForPca(samplesAsColsDf, designDf,
                                   sampleNameDesignColName)
  pcaResults = doAndPrintScaledPcaOnSamplesAsRowsDf(transformedDf)
  return(pcaResults)
}

doScaledPcaOnSamplesAsRowsDf<-function(data_df){
  # remove any columns that are constant in order to allow scaling
  variable_df = data_df[,apply(data_df, 2, var, na.rm=TRUE) != 0]
  pca_result = prcomp(variable_df, scale = TRUE)
}

isPointOutsideEllipse<-function(x, y, ellipseCenterAndRadii){
  ellipseEqnValue = ((x - ellipseCenterAndRadii[1])^2)/(
    (ellipseCenterAndRadii[3])^2) +
    ((y - ellipseCenterAndRadii[2])^2)/((ellipseCenterAndRadii[4])^2)
  return(ellipseEqnValue > 1)
}

getWhetherPointsAreOutliers<-function(pcaResults){
  xVals = pcaResults$x[,1]
  yVals = pcaResults$x[,2]
  ellipseInfo = getHotellingT2EllipseCenterAndRadii(xVals, yVals)
  isOutsideEllipse = mapply(isPointOutsideEllipse, xVals, yVals,
                            MoreArgs=list(ellipseCenterAndRadii=ellipseInfo))
  return(isOutsideEllipse)
}

getHotellingT2Ellipse <-function (x, y, alfa = 0.95, len = 200) {
  ellipseInfo = getHotellingT2EllipseCenterAndRadii(x, y, alfa)
  mypi <- seq(0, 2 * pi, length = len)
  r1 = ellipseInfo[3]
  r2 = ellipseInfo[4]
  cbind(r1 * cos(mypi) + mean(x), r2 * sin(mypi) + mean(y))
}

getHotellingT2EllipseCenterAndRadii<-function(x, y, alfa = 0.95){
  # NOTE: this logic, except for the return statement,
  # is a trimmed COPY-PASTE of the simpleEllipse method in the
  # pcaMethods package.  However, although the authors included
  # this method in the documentation as public, they forgot to
  # *make* it public.  Also, that package
  # seems to bog down my notebook for unknown reasons.
  N <- length(x)
  r1 <- sqrt(var(x) * qf(alfa, 2, N - 2) * (2 * (N^2 - 1)/(N *
                                                             (N - 2))))
  r2 <- sqrt(var(y) * qf(alfa, 2, N - 2) * (2 * (N^2 - 1)/(N *
                                                             (N - 2))))
  return(c(mean(x), mean(y),r1,r2))
}


expandPlot<-function(aPlot, additiveExpandValue=25,
                     shrink_viewport=FALSE){
  aPlot = aPlot + scale_x_continuous(
    expand =(c(0.05,additiveExpandValue)))
  printPlotInViewport(aPlot, shrink_viewport=FALSE)
}


