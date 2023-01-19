# ---------------------------------------------------------------------------------
# Copyright (c) 2018 UC San Diego Center for Computational Biology & Bioinformatics
#
# Distributed under the terms of the MIT License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------------
# Initial author: Amanda Birmingham

makeRunName = function(gProjectName, gStepName){
	return(paste0(gProjectName, "_", gStepName, "_", gsub("[: -]", "", strptime(Sys.time(), "%Y-%m-%d %H:%M:%S"), perl=TRUE)))
}

writeWorkspaceImage = function(outputDir, runName){
    fileName = sprintf("%s.RData",runName)
    save.image(file=file.path(outputDir, fileName))
    print(paste0("Output file: ",fileName))
}

# from https://www.r-bloggers.com/safe-loading-of-rdata-files-2/
loadToEnvironment <- function(RData, env = new.env()){
  load(RData, env)
  return(env)
}
