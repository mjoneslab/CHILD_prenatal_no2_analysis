require(genefilter)
require(quadprog)

estimateCellCounts2v2 <- function (rgSet, 
                                   compositeCellType = "Blood", 
                                   processMethod = "preprocessNoob",
                                   probeSelect = "IDOL", 
                                   cellTypes = c("CD8T", "CD4T", "NK", 
                                                 "Bcell", "Mono", "Neu"),  
                                   referencePlatform = c("IlluminaHumanMethylation450k", 
                                                         "IlluminaHumanMethylationEPIC", 
                                                         "IlluminaHumanMethylation27k"), 
                                   referenceset = NULL, CustomCpGs = NULL, returnAll = FALSE, 
                                   meanPlot = FALSE, verbose = TRUE, lessThanOne = FALSE, cellcounts = NULL, 
                                   ...) 
{
  
  # if rgset is an RGChannelSetExtended
  if (is(rgSet, "RGChannelSetExtended")) {
    rgSet <- as(rgSet, "RGChannelSet")
  }
  
  # if rgset is NOT an RGChannelSet or MethylSet
  if ((!is(rgSet, "RGChannelSet")) && (!is(rgSet, "MethylSet"))) {
    stop(strwrap(sprintf("object is of class '%s', but needs to be of\n                                class 'RGChannelSet' 'RGChannelSetExtended' or\n                                'MethylSet' to use this function", 
                         class(rgSet)), width = 80, prefix = " ", initial = ""))
  }
  
  # if rgset is NOT an RGChannelSet and preprocessMethod is not preproprocessQuantile
  if (!is(rgSet, "RGChannelSet") && (processMethod[1] != "preprocessQuantile")) {
    stop(strwrap(sprintf("object is of class '%s', but needs to be of\n                                class 'RGChannelSet' or 'RGChannelSetExtended'\n                                to use other methods different to\n                                'preprocessQuantile'", 
                         class(rgSet)), width = 80, prefix = " ", initial = ""))
  }
  
  # if rgSet is an MethylSet and preprocessMethod is equal preprocessingQuantile
  if (is(rgSet, "MethylSet") && (processMethod[1] == "preprocessQuantile")) {
    message(strwrap("[estimateCellCounts2] The function will assume that\n                            no preprocessing has been performed. Using\n                            'preprocessQuantile' in prenormalized data is\n                            experimental and it should only be run under the\n                            user responsibility", 
                    width = 80, prefix = " ", initial = ""))
  }
  
  # get the user set referencePlatform
  referencePlatform <- match.arg(referencePlatform, 
                                 choices = c("IlluminaHumanMethylation450k", 
                                             "IlluminaHumanMethylationEPIC", 
                                             "IlluminaHumanMethylation27k"))
  
  # get the rgPlatform of user inputted data (e.g. 450k, EPIC, or 27k)
  rgPlatform <- sub("IlluminaHumanMethylation", "", annotation(rgSet)[which(names(annotation(rgSet)) == 
                                                                              "array")])
  # get the reference platform of the reference set
  platform <- sub("IlluminaHumanMethylation", "", referencePlatform)
  
  # check if cell type is related to cord blood
  if ((compositeCellType == "CordBlood" | compositeCellType == 
       "CordBloodNorway" | compositeCellType == "CordBloodCanada" | 
       compositeCellType == "CordTissueAndBlood")) {
    message(strwrap("[estimateCellCounts2] Consider using CordBloodCombined\n                        for cord blood estimation.\n", 
                    width = 80, prefix = " ", initial = ""))
  }
  
  # check if cell tpe is related to blood
  if ((compositeCellType == "Blood") && (referencePlatform == 
                                         "IlluminaHumanMethylation450k")) {
    message(strwrap("[estimateCellCounts2] Consider using referencePlatform\n                        IlluminaHumanMethylationEPIC for blood estimation,\n                        using the IDOLOptimizedCpGs450klegacy\n", 
                    width = 80, prefix = " ", initial = ""))
  }
  
  # check if cell type is related to blood extended and if EPIC is specified as it is only platform supported
  if ((compositeCellType == "BloodExtended") && (referencePlatform == 
                                                 "IlluminaHumanMethylation450k")) {
    message(strwrap("[estimateCellCounts2] Platform 450k is not available\n                        defaulting to platform EPIC for your estimation.\n", 
                    width = 80, prefix = " ", initial = ""))
    referencePlatform <- "IlluminaHumanMethylationEPIC"
    platform <- sub("IlluminaHumanMethylation", "", referencePlatform)
  }
  
  # check if cell type is cordbloodcombined and if reference platform is 450k
  if ((compositeCellType == "CordBloodCombined") && (referencePlatform == 
                                                     "IlluminaHumanMethylationEPIC")) {
    message(strwrap("[estimateCellCounts2] Platform EPIC is not available\n                        defaulting to platform 450k for your estimation.\n", 
                    width = 80, prefix = " ", initial = ""))
    referencePlatform <- "IlluminaHumanMethylation450k"
    platform <- sub("IlluminaHumanMethylation", "", referencePlatform)
  }
  
  # check if cell type is of type cord blood and if nRBCs are included
  if ((compositeCellType == "CordBlood" | compositeCellType == 
       "CordBloodCombined") && (!"nRBC" %in% cellTypes)) {
    message(strwrap("[estimateCellCounts2] Consider including 'nRBC' in\n                        argument 'cellTypes' for cord blood estimation.\n", 
                    width = 80, prefix = " ", initial = ""))
  }
  
  # check if cell type is of type blood and if it contains gran and neu
  if ((compositeCellType == "Blood" | compositeCellType == 
       "BloodExtended") && (referencePlatform == "IlluminaHumanMethylationEPIC") && 
      ("Gran" %in% cellTypes)) {
    message(strwrap("[estimateCellCounts2] Replace 'Gran' for 'Neu' in\n                        argument 'cellTypes' for EPIC blood estimation.\n", 
                    width = 80, prefix = " ", initial = ""))
  }
  
  # check if cell type is blood extended and if it contains new cell types
  if ((compositeCellType == "BloodExtended") && (referencePlatform == 
                                                 "IlluminaHumanMethylationEPIC") && ("CD4T" %in% cellTypes)) {
    message(strwrap("[estimateCellCounts2] Replace 'CD4T' for 'CD4nv',\n                        'CD4mem' and 'Treg' in argument 'cellTypes' for EPIC\n                        blood extended estimation.\n", 
                    width = 80, prefix = " ", initial = ""))
  }
  
  # check if cell type is blood extended and if it contains new cell types
  if ((compositeCellType == "BloodExtended") && (referencePlatform == 
                                                 "IlluminaHumanMethylationEPIC") && ("CD8T" %in% cellTypes)) {
    message(strwrap("[estimateCellCounts2] Replace 'CD8T' for 'CD8nv' and\n                        'CD8mem' in argument 'cellTypes' for EPIC blood\n                        extended estimation.\n", 
                    width = 80, prefix = " ", initial = ""))
  }
  
  # check if cell type is blood extended and if it contains new cell types
  if ((compositeCellType == "BloodExtended") && (referencePlatform == 
                                                 "IlluminaHumanMethylationEPIC") && ("Bcell" %in% cellTypes)) {
    message(strwrap("[estimateCellCounts2] Replace 'Bcell' for 'Bnv' and\n                        'Bmem' in argument 'cellTypes' for EPIC blood\n                        extended estimation.\n", 
                    width = 80, prefix = " ", initial = ""))
  }
  if ((compositeCellType != "Blood" && compositeCellType != 
       "BloodExtended") && ("Neu" %in% cellTypes)) {
    message(strwrap("[estimateCellCounts2] Check whether 'Gran' or 'Neu' is\n                        present in your reference and adjust argument\n                        'cellTypes' for your estimation.\n", 
                    width = 80, prefix = " ", initial = ""))
  }
  
  # check if cell type is blood extended and if it contains new cell types
  if ((compositeCellType == "Blood" | compositeCellType == 
       "BloodExtended") && ("Gran" %in% cellTypes)) {
    message(strwrap("[estimateCellCounts2] Check whether 'Gran' or 'Neu' is\n                        present in your reference and adjust argument\n                        'cellTypes' for your estimation.\n", 
                    width = 80, prefix = " ", initial = ""))
  }
  
  # check if cell type is blood extended and if it contains new cell types
  if ((compositeCellType != "BloodExtended") && ("Eos" %in% 
                                                 cellTypes)) {
    message(strwrap("[estimateCellCounts2] Check whether 'Eos' is\n                        present in your reference and adjust argument\n                        'cellTypes' for your estimation.\n", 
                    width = 80, prefix = " ", initial = ""))
  }
  
  # check if cell type is blood extended and if it contains new cell types
  if ((compositeCellType != "BloodExtended") && ("Bas" %in% 
                                                 cellTypes)) {
    message(strwrap("[estimateCellCounts2] Check whether 'Bas' is\n                        present in your reference and adjust argument\n                        'cellTypes' for your estimation.\n", 
                    width = 80, prefix = " ", initial = ""))
  }
  
  # check length of cell count (if not NULL) data vs number of samples (should be same)
  if (!is.null(cellcounts) && length(cellcounts) != dim(rgSet)[2]) {
    stop(strwrap(sprintf("length of cellcounts is '%s', but needs to be equal\n                        to the number of samples '%s'", 
                         length(cellcounts), dim(rgSet)[2], width = 80, prefix = " ", 
                         initial = "")))
  }
  
  # specify reference package
  referencePkg <- sprintf("FlowSorted.%s.%s", compositeCellType, 
                          platform)
  
  subverbose <- max(as.integer(verbose) - 1L, 0L)
  
  # if referenceset is not NULL
  if (!is.null(referenceset)) {
    referenceRGset <- get(referenceset)
    # if rgset is RGChannelSetExtended
    if (is(referenceRGset, "RGChannelSetExtended")) {
      referenceRGset <- as(referenceRGset, "RGChannelSet")
    }
    if (!is(rgSet, "RGChannelSet")) {
      referenceRGset <- preprocessRaw(referenceRGset)
    }
  }
  # if reference set is NULL
  else {
    if (!require(referencePkg, character.only = TRUE) && 
        referencePkg != "FlowSorted.BloodExtended.EPIC") {
      stop(strwrap(sprintf("Could not find reference data package for\n                                compositeCellType '%s' and referencePlatform\n                                '%s' (inferred package name is '%s')", 
                           compositeCellType, platform, referencePkg), width = 80, 
                   prefix = " ", initial = ""))
    }
    if (!require(referencePkg, character.only = TRUE) && 
        referencePkg == "FlowSorted.BloodExtended.EPIC") {
      stop(strwrap(sprintf("Could not find reference data package for\n                                compositeCellType '%s' and referencePlatform\n                                '%s' (inferred package name is '%s'),\n                                please contact\n                                Technology.Transfer@dartmouth.edu", 
                           compositeCellType, platform, referencePkg), width = 80, 
                   prefix = " ", initial = ""))
    }
    if ((referencePkg != "FlowSorted.Blood.EPIC") && (referencePkg != 
                                                      "FlowSorted.CordBloodCombined.450k")) {
      referenceRGset <- get(referencePkg)
    }
    else if ((referencePkg == "FlowSorted.Blood.EPIC") | 
             (referencePkg == "FlowSorted.CordBloodCombined.450k")) {
      referenceRGset <- libraryDataGet(referencePkg)
    }
    if (!is(rgSet, "RGChannelSet")) {
      referenceRGset <- preprocessRaw(referenceRGset)
    }
  }
  
  # if platform is not the same between user data and reference set
  if (rgPlatform != platform) {
    rgSet <- convertArray(rgSet, outType = referencePlatform, 
                          verbose = TRUE)
  }
  
  # if cell type not specified
  if (!"CellType" %in% names(colData(referenceRGset))) {
    stop(strwrap(sprintf("the reference sorted dataset (in this case '%s')\n                            needs to have a phenoData column called\n                            'CellType'"), 
                 names(referencePkg), width = 80, prefix = " ", initial = ""))
  }
  
  # if naming issue exists
  if (sum(colnames(rgSet) %in% colnames(referenceRGset)) > 
      0) {
    stop(strwrap("the sample/column names in the user set must not be in\n                    the reference data ", 
                 width = 80, prefix = " ", initial = ""))
  }
  
  # if cell types are not the same
  if (!all(cellTypes %in% referenceRGset$CellType)) {
    stop(strwrap(sprintf("all elements of argument 'cellTypes' needs to be\n                            part of the reference phenoData columns 'CellType'\n                            (containg the following elements: '%s')", 
                         paste(unique(referenceRGset$cellType), collapse = "', '")), 
                 width = 80, prefix = " ", initial = ""))
  }
  
  # if 1 or fewer cell types are provided
  if (length(unique(cellTypes)) < 2) {
    stop("At least 2 cell types must be provided.")
  }
  
  # if preprocessmethod is auto and cell types is of Blood of DLPFC
  if ((processMethod == "auto") && (compositeCellType %in% 
                                    c("Blood", "DLPFC"))) {
    processMethod <- "preprocessQuantile"
  }
  
  # if preprocessmethod is author and cell type is not blood or DLPFC and rgset is rgchannelset
  if ((processMethod == "auto") && (!compositeCellType %in% 
                                    c("Blood", "DLPFC")) && (is(rgSet, "RGChannelSet"))) {
    processMethod <- "preprocessNoob"
  }
  
  # get preprocessmethod
  processMethod <- get(processMethod)
  
  # if probe select is "auto" and cell types are of cord blood set any
  if ((probeSelect == "auto") && (compositeCellType %in% c("CordBloodCombined", 
                                                           "CordBlood", "CordBloodNorway", "CordTissueAndBlood"))) {
    probeSelect <- "any"
  }
  
  # if probe select is "auto" and cell types are of NOT cord blood set both
  if ((probeSelect == "auto") && (!compositeCellType %in% c("CordBloodCombined", 
                                                            "CordBlood", "CordBloodNorway", "CordTissueAndBlood"))) {
    probeSelect <- "both"
  }
  
  # if probe select is "IDOL" and cell types cordbloodcombined
  if ((probeSelect == "IDOL") && (compositeCellType %in% c("CordBloodCombined", 
                                                           "CordBlood", "CordBloodNorway", "CordTissueAndBlood"))) {
    CustomCpGs <- FlowSorted.Blood.EPIC::IDOLOptimizedCpGsCordBlood
  }
  
  # if probe select is "IDOL" and cell types blood
  if ((probeSelect == "IDOL") && (compositeCellType == "Blood")) {
    if ((rgPlatform == "450k")) {
      CustomCpGs <- FlowSorted.Blood.EPIC::IDOLOptimizedCpGs450klegacy
    }
    else {
      CustomCpGs <- FlowSorted.Blood.EPIC::IDOLOptimizedCpGs
    }
  }
  
  if (verbose) {
    message(strwrap("[estimateCellCounts2] Combining user data with\n                        reference (flow sorted) data.\n", 
                    width = 80, prefix = " ", initial = ""))
  }
  
  # create new data frame
  newpd <- DataFrame(sampleNames = c(colnames(rgSet), colnames(referenceRGset)), 
                     studyIndex = rep(c("user", "reference"), times = c(ncol(rgSet),ncol(referenceRGset))))
  
  # recode refence set cell types as character
  referenceRGset$CellType <- as.character(referenceRGset$CellType)
  
  # if cell type is null in rgset set NA
  if (is.null(rgSet$CellType)) {
    rgSet$CellType <- rep("NA", dim(rgSet)[2])
  }

  # if age is null in rgset set NA
  if (is.null(rgSet$Age)) {
    rgSet$Age <- rep("NA", dim(rgSet)[2])
  }
  
  # if sex is null in rgset set NA
  if (is.null(rgSet$Sex)) {
    rgSet$Sex <- rep("NA", dim(rgSet)[2])
  }
  
  # if sex is null in referenceset set NA
  if (is.null(referenceRGset$Sex)) {
    referenceRGset$Sex <- rep("NA", dim(referenceRGset)[2])
  }
  else {
    referenceRGset$Sex <- as.character(referenceRGset$Sex)
  }
  
  # if age is null in referenceset set NA
  if (is.null(referenceRGset$Age)) {
    referenceRGset$Age <- rep("NA", dim(referenceRGset)[2])
  }
  else {
    try(referenceRGset$Age <- as.numeric(referenceRGset$Age))
  }
  
  # get shared column names between reference set and rgset
  commoncolumn <- intersect(names(colData(rgSet)), names(colData(referenceRGset)))
  
  restry <- try({
    colData(rgSet)[commoncolumn] <- mapply(FUN = as, colData(rgSet)[commoncolumn], 
                                           vapply(colData(referenceRGset)[commoncolumn], class, 
                                                  FUN.VALUE = character(1)), SIMPLIFY = FALSE)
  }, silent = TRUE)
  
  if ("try-error" %in% class(restry)) {
    commoncolumn <- c("CellType", "Sex", "Age")
    colData(rgSet)[commoncolumn] <- mapply(FUN = as, colData(rgSet)[commoncolumn], 
                                           vapply(colData(referenceRGset)[commoncolumn], class, 
                                                  FUN.VALUE = character(1)), SIMPLIFY = FALSE)
  }
  else {
    colData(rgSet)[commoncolumn] <- mapply(FUN = as, colData(rgSet)[commoncolumn], 
                                           vapply(colData(referenceRGset)[commoncolumn], class, 
                                                  FUN.VALUE = character(1)), SIMPLIFY = FALSE)
  }
  
  rm(restry)
  colData(referenceRGset) <- colData(referenceRGset)[commoncolumn]
  colData(rgSet) <- colData(rgSet)[commoncolumn]
  referencePd <- colData(referenceRGset)
  combinedRGset <- combineArrays(rgSet, referenceRGset, outType = referencePlatform)
  colData(combinedRGset) <- newpd
  colnames(combinedRGset) <- newpd$sampleNames
  rm(referenceRGset)
  
  # output a message if set to verbose
  if (verbose) {
    message(strwrap("[estimateCellCounts2] Processing user and reference\n                        data together.\n", 
                    width = 80, prefix = " ", initial = ""))
  }
  
  # preprocess reference and rgset combined together
  if (compositeCellType == "CordBlood") {
    # if combined is not an rgset
    if (!is(combinedRGset, "RGChannelSet")) {
      combinedRGset@preprocessMethod["rg.norm"] <- "Raw (no normalization or bg correction)"
    }
    
    combinedMset <- preprocessNoob(combinedRGset, verbose = subverbose) 
    # preprocessNoob was formally "preprocessMethod"
    # was not calling appropriately
    rm(combinedRGset)
    gc()
    compTable <- get(paste0(referencePkg, ".compTable"))
    combinedMset <- combinedMset[which(rownames(combinedMset) %in% 
                                         rownames(compTable)), ]
  }
  else {
    if (!is(combinedRGset, "RGChannelSet")) {
      combinedRGset@preprocessMethod["rg.norm"] <- "Raw (no normalization or bg correction)"
    }
    combinedMset <- processMethod(combinedRGset)
    rm(combinedRGset)
    gc()
  }
  
  
  referenceMset <- combinedMset[, combinedMset$studyIndex == "reference"]
  colData(referenceMset) <- as(referencePd, "DataFrame")
  mSet <- combinedMset[, combinedMset$studyIndex == "user"]
  colData(mSet) <- as(colData(rgSet), "DataFrame")
  rm(combinedMset)
  
  # if probe select is not IDOL
  if (probeSelect != "IDOL") {
    if (verbose) {
      message(strwrap("[estimateCellCounts2] Picking probes for\n                            composition estimation.\n", 
                      width = 80, prefix = " ", initial = ""))
    }
    compData <<- pickCompProbes(referenceMset, cellTypes = cellTypes, 
                               compositeCellType = compositeCellType, probeSelect = probeSelect)
    coefs <- compData$coefEsts
    if (verbose) {
      message(strwrap("[estimateCellCounts2] Estimating  proportion\n                            composition (prop), if you provide cellcounts\n                            those will be provided as counts in the\n                            composition estimation.\n", 
                      width = 80, prefix = " ", initial = ""))
    }
    prop <- projectCellType_CP(getBeta(mSet)[rownames(coefs),], coefs, lessThanOne = lessThanOne)
    prop <- round(prop, 4)
    rownames(prop) <- colnames(rgSet)
    counts <- round(prop * cellcounts, 0)
    if (meanPlot) {
      smeans <- compData$sampleMeans
      smeans <- smeans[order(names(smeans))]
      sampleMeans <- c(colMeans(getBeta(mSet)[rownames(coefs), 
      ]), smeans)
      sampleColors <- c(rep(1, ncol(mSet)), 1 + as.numeric(factor(names(smeans))))
      plot(sampleMeans, pch = 21, bg = sampleColors)
      legend("bottomleft", c("blood", levels(factor(names(smeans)))), 
             col = seq_len(7), pch = 15)
    }
    if (returnAll) {
      list(prop = prop, counts = counts, compTable = compData$compTable, 
           normalizedData = mSet)
    }
    else {
      list(prop = prop, counts = counts)
    }
  }
  # if probe select is IDOL
  else {
    if (verbose) {
      message(strwrap("[estimateCellCounts2] Using IDOL L-DMR probes for\n                            composition estimation.\n", 
                      width = 80, prefix = " ", initial = ""))
    }
    p <- getBeta(referenceMset)
    pd <- as.data.frame(colData(referenceMset))
    rm(referenceMset)
    if (!is.null(cellTypes)) {
      if (!all(cellTypes %in% pd$CellType)) {
        stop(strwrap("elements of argument 'cellTypes' are not part of\n                            'referenceMset$CellType'", 
                     width = 80, prefix = " ", initial = ""))
      }
      keep <- which(pd$CellType %in% cellTypes)
      pd <- pd[keep, ]
      p <- p[, keep]
    }
    pd$CellType <- factor(pd$CellType, levels = cellTypes)
    ffComp <- genefilter::rowFtests(p, pd$CellType)
    tIndexes <- split(seq(along = pd$CellType), pd$CellType)
    prof <- vapply(tIndexes, function(i) rowMeans(p[, i]), 
                   FUN.VALUE = numeric(dim(p)[1]))
    r <- rowRanges(p)
    compTable <- cbind(ffComp, prof, r, abs(r[, 1] - r[, 
                                                       2]))
    names(compTable)[1] <- "Fstat"
    names(compTable)[c(-2, -1, 0) + ncol(compTable)] <- c("low", "high", "range")
    
    tstatList <- lapply(tIndexes, function(i) {
      x <- rep(0, ncol(p))
      x[i] <- 1
      return(rowttests(p, factor(x)))
    })
    trainingProbes <- CustomCpGs
    trainingProbes <- trainingProbes[trainingProbes %in% 
                                       rownames(p)]
    p <- p[trainingProbes, ]
    pMeans <- colMeans(p)
    names(pMeans) <- pd$CellType
    form <- as.formula(sprintf("y ~ %s - 1", paste(levels(pd$CellType), 
                                                   collapse = "+")))
    phenoDF <- as.data.frame(model.matrix(~pd$CellType - 
                                            1))
    colnames(phenoDF) <- sub("^pd\\$CellType", "", colnames(phenoDF))
    
    # if only 2 cell types specified
    if (ncol(phenoDF) == 2) {
      X <- as.matrix(phenoDF)
      coefEsts <- t(solve(t(X) %*% X) %*% t(X) %*% t(p))
      coefs <- coefEsts
    }
    # if >2 cell types specified
    else {
      ###########################
      # PROBLEM HERE ############
      ###########################
      
      tmp <- validationCellType(Y = p, 
                                pheno = phenoDF, 
                                modelFix = form,
                                modelBatch = NULL,
                                L.forFstat = NULL, 
                                verbose = FALSE)
      
      coefEsts <- tmp$coefEsts
      coefs <- coefEsts
    }
    compData <- list(coefEsts = coefEsts, compTable = compTable, 
                     sampleMeans = pMeans)
    if (verbose) {
      message(strwrap("[estimateCellCounts2] Estimating  proportion\n                            composition (prop), if you provide cellcounts\n                            those will be provided as counts in the\n                            composition estimation.\n", 
                      width = 80, prefix = " ", initial = ""))
    }
    prop <- projectCellType_CP(getBeta(mSet)[rownames(coefs), 
    ], coefs, lessThanOne = lessThanOne)
    prop <- round(prop, 4)
    rownames(prop) <- colnames(rgSet)
    counts <- round(prop * cellcounts, 0)
    if (meanPlot) {
      smeans <- compData$sampleMeans
      smeans <- smeans[order(names(smeans))]
      sampleMeans <- c(colMeans(getBeta(mSet)[rownames(coefs), 
      ]), smeans)
      sampleColors <- c(rep(1, ncol(mSet)), 1 + as.numeric(factor(names(smeans))))
      plot(sampleMeans, pch = 21, bg = sampleColors)
      legend("bottomleft", c("blood", levels(factor(names(smeans)))), 
             col = seq_len(7), pch = 15)
    }
    if (returnAll) {
      list(prop = prop, counts = counts, compTable = compTable, 
           normalizedData = mSet)
    }
    else {
      list(prop = prop, counts = counts)
    }
  }
}

projectCellType_CP <- function(Y, coefWBC, contrastWBC = NULL,
                               
                               nonnegative = TRUE, lessThanOne = FALSE) {
  
  if (is.null(contrastWBC)) {
    
    Xmat <- coefWBC
    
  } else {
    
    Xmat <- tcrossprod(coefWBC, contrastWBC)
    
  }
  
  nCol <- dim(Xmat)[2]
  
  if (nCol == 2) {
    
    Dmat <- crossprod(Xmat)
    
    mixCoef <- t(apply(Y, 2, function(x) {
      
      solve(Dmat, crossprod(Xmat, x))
      
    }))
    
    colnames(mixCoef) <- colnames(Xmat)
    
    return(mixCoef)
    
  } else {
    
    nSubj <- dim(Y)[2]
    
    mixCoef <- matrix(0, nSubj, nCol)
    
    rownames(mixCoef) <- colnames(Y)
    
    colnames(mixCoef) <- colnames(Xmat)
    
    if (nonnegative) {
      
      if (lessThanOne) {
        
        Amat <- cbind(rep(-1, nCol), diag(nCol))
        
        b0vec <- c(-1, rep(0, nCol))
        
      } else {
        
        Amat <- diag(nCol)
        
        b0vec <- rep(0, nCol)
        
      }
      
      for (i in seq_len(nSubj)) {
        
        obs <- which(!is.na(Y[, i]))
        
        Dmat <- crossprod(Xmat[obs, ])
        
        mixCoef[i, ] <- round(quadprog::solve.QP(
          
          Dmat, crossprod(
            
            Xmat[obs, ],
            
            Y[obs, i]
            
          ),
          
          Amat, b0vec
          
        )$sol, 4)
        
      }
      
    } else {
      
      for (i in seq_len(nSubj)) {
        
        obs <- which(!is.na(Y[, i]))
        
        Dmat <- crossprod(Xmat[obs, ])
        
        mixCoef[i, ] <- round(solve(Dmat, t(Xmat[obs, ]) %*% Y[obs, i]), 4)
        
      }
      
    }
    
    return(mixCoef)
    
  }
  
}

libraryDataGet <- function(title) {
  
  assign(title, ExperimentHub()[[query(
    
    ExperimentHub(),
    
    title
    
  )$ah_id]])
  
}

pickCompProbes <- function(mSet, cellTypes = NULL, numProbes = 50,
                           compositeCellType = compositeCellType,
                           probeSelect = probeSelect) {
  p <- getBeta(mSet)
  pd <- as.data.frame(colData(mSet))
  if (!is.null(cellTypes)) {
    if (!all(cellTypes %in% pd$CellType)) {
      stop(strwrap("elements of argument 'cellTypes' are not part of
                        'mSet$CellType'",
                   width = 80, prefix = " ",
                   initial = ""
      ))
    }
    keep <- which(pd$CellType %in% cellTypes)
    pd <- pd[keep, ]
    p <- p[, keep]
  }
  ## make cell type a factor
  pd$CellType <- factor(pd$CellType, levels = cellTypes)
  ffComp <- genefilter::rowFtests(p, pd$CellType)
  tIndexes <- split(seq(along = pd$CellType), pd$CellType)
  prof <- vapply(tIndexes, function(i) rowMeans(p[, i]),
                 FUN.VALUE = numeric(dim(p)[1])
  )
  r <- rowRanges(p)
  compTable <- cbind(ffComp, prof, r, abs(r[, 1] - r[, 2]))
  names(compTable)[1] <- "Fstat"
  names(compTable)[c(-2, -1, 0) + ncol(compTable)] <- c("low", "high", "range")
  tstatList <- lapply(tIndexes, function(i) {
    x <- rep(0, ncol(p))
    x[i] <- 1
    return(rowttests(p, factor(x)))
  })
  if (probeSelect == "any") {
    probeList <- lapply(tstatList, function(x) {
      y <- x[x[, "p.value"] < 1e-8, ]
      yAny <- y[order(abs(y[, "dm"]), decreasing = TRUE), ]
      c(rownames(yAny)[seq_len(numProbes * 2)])
    })
    probelist_all <<- probeList
  } else {
    probeList <- lapply(tstatList, function(x) {
      y <- x[x[, "p.value"] < 1e-8, ]
      yUp <- y[order(y[, "dm"], decreasing = TRUE), ]
      yDown <- y[order(y[, "dm"], decreasing = FALSE), ]
      c(
        rownames(yUp)[seq_len(numProbes)],
        rownames(yDown)[seq_len(numProbes)]
      )
    })
    probelist_all <<- probeList
  }
  trainingProbes <- unique(unlist(probeList))
  p <- p[trainingProbes, ]
  
  pMeans <- colMeans(p)
  names(pMeans) <- pd$CellType
  form <- as.formula(sprintf("y ~ %s - 1", paste(levels(pd$CellType),
                                                 collapse = "+"
  )))
  phenoDF <- as.data.frame(model.matrix(~ pd$CellType - 1))
  colnames(phenoDF) <- sub("^pd\\$CellType", "", colnames(phenoDF))
  if (ncol(phenoDF) == 2) { # two group solution
    X <- as.matrix(phenoDF)
    coefEsts <- t(solve(t(X) %*% X) %*% t(X) %*% t(p))
  } else { # > 2 group solution
    tmp <- validationCellType(Y = p, pheno = phenoDF, modelFix = form)
    coefEsts <- tmp$coefEsts
  }
  
  out <- list(
    coefEsts = coefEsts, compTable = compTable,
    sampleMeans = pMeans
  )
  return(out)
}

validationCellType <- function(Y, pheno, modelFix, modelBatch = NULL,
                               L.forFstat = NULL, verbose = FALSE) {
  print("validating")
  N <- dim(pheno)[1]
  pheno$y <- rep(0, N)
  xTest <- model.matrix(modelFix, pheno)
  sizeModel <- dim(xTest)[2]
  M <- dim(Y)[1]
  if (is.null(L.forFstat)) {
    print("NULL L.forFstat")
    L.forFstat <- diag(sizeModel)[-1, ] # All non-intercept coefficients
    colnames(L.forFstat) <- colnames(xTest)
    rownames(L.forFstat) <- colnames(xTest)[-1]
  }
  # Initialize various containers
  sigmaResid <- rep(NA, M)
  sigmaIcept <- rep(NA, M)
  nObserved <- rep(NA, M)
  nClusters <-  rep(NA, M)
  Fstat <- rep(NA, M)
  coefEsts <- matrix(NA, M, sizeModel)
  coefVcovs <- list()
  if (verbose) {
    cat("[validationCellType] ")
  }
  
  for (j in seq_len(M)) { # For each CpG
    print("seq_len(M)")
    
    ## Remove missing methylation values
    ii <- !is.na(Y[j, ])
    nObserved[j] <- sum(ii)
    pheno$y <- Y[j, ]
    
    if (j %% round(M / 10) == 0 && verbose) {
      cat(".")
    } # Report progress
    
    if (!is.null(modelBatch)) {
      print("null model batch")
      
      fit <- try(lme(modelFix, random = modelBatch, data = pheno[ii, ]))
      OLS <- inherits(fit, "try-error")
      # If LME can't be fit, just use OLS
    } else {
      print("set ols to true")
      
      OLS <- TRUE
    }
    
    if (OLS) {
      print("for true ols")
      
      fit <- lm(modelFix, data = pheno[ii, ])
      fitCoef <<- fit$coef
      sigmaResid[j] <- summary(fit)$sigma
      sigmaIcept[j] <- 0
      nClusters[j] <- 0
    } else {
      print("ols not tru")
      
      fitCoef <<- fit$coef$fixed
      sigmaResid[j] <- fit$sigma
      sigmaIcept[j] <- sqrt(getVarCov(fit)[1])
      nClusters[j] <- length(fit$coef$random[[1]])
    }
    coefEsts[j, ] <- fitCoef
    coefVcovs[[j]] <- vcov(fit)
    
    useCoef <- L.forFstat %*% fitCoef
    useV <- L.forFstat %*% coefVcovs[[j]] %*% t(L.forFstat)
    Fstat[j] <- (t(useCoef) %*% solve(useV, useCoef)) / sizeModel
    
  }
  if (verbose) {
    cat(" done\n")
  }
  ## Name the rows so that they can be easily matched to the target data set
  rownames(coefEsts) <- rownames(Y)
  colnames(coefEsts) <- names(fitCoef)
  degFree <- nObserved - nClusters - sizeModel + 1
  
  ## Get P values corresponding to F statistics
  Pval <- 1 - pf(Fstat, sizeModel, degFree)
  
  out <- list(
    coefEsts = coefEsts, coefVcovs = coefVcovs, modelFix = modelFix,
    modelBatch = modelBatch,
    sigmaIcept = sigmaIcept, sigmaResid = sigmaResid,
    L.forFstat = L.forFstat, Pval = Pval,
    orderFstat = order(-Fstat), Fstat = Fstat, nClusters = nClusters,
    nObserved = nObserved,
    degFree = degFree
  )
  
  out
}
