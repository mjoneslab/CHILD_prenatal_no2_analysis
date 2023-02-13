rgSet = rgSet
compositeCellType = "Blood"
processMethod = "preprocessNoob"
probeSelect = c("auto", "any", "IDOL")
cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")
referencePlatform = c("IlluminaHumanMethylation450k", "IlluminaHumanMethylationEPIC",
                      "IlluminaHumanMethylation27k")
referenceset = NULL
IDOLOptimizedCpGs = NULL
returnAll = FALSE
meanPlot = FALSE
verbose = TRUE
lessThanOne = FALSE



library(ExperimentHub)
hub <- ExperimentHub()
query(hub, "FlowSorted.CordBloodCombined.450k")
FlowSorted.Blood.EPIC <- hub[["EH1136"]]
FlowSorted.Blood.EPIC
