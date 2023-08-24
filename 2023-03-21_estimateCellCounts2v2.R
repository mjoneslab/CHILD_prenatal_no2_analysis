require(genefilter)
require(quadprog)

estimateCellCounts2v2 <- function(rgSet, compositeCellType = "Blood",
                                  processMethod = "preprocessNoob",
                                  probeSelect = "IDOL",
                                  cellTypes = c(
                                    "CD8T", "CD4T", "NK", "Bcell",
                                    "Mono", "Neu"
                                  ),
                                  referencePlatform = "IlluminaHumanMethylationEPIC",
                                  referenceset = NULL, CustomCpGs = NULL,
                                  returnAll = FALSE, meanPlot = FALSE,
                                  verbose = TRUE, lessThanOne = FALSE,
                                  cellcounts = NULL,
                                  ...) {
  
  # check what class rgSet is 
  if (is(rgSet, "RGChannelSetExtended")) {
    rgSet <- as(rgSet, "RGChannelSet")
  }
  
  if ((!is(rgSet, "RGChannelSet")) && (!is(rgSet, "MethylSet"))) {
    stop(strwrap(sprintf(
      "object is of class '%s', but needs to be of
                                class 'RGChannelSet' 'RGChannelSetExtended' or
                                'MethylSet' to use this function",
      class(rgSet)
    ),
    width = 80, prefix = " ",
    initial = ""
    ))
  }
  
  if (!is(rgSet, "RGChannelSet") &&
      (processMethod[1] != "preprocessQuantile")) {
    stop(strwrap(sprintf("object is of class '%s', but needs to be of
                                class 'RGChannelSet' or 'RGChannelSetExtended'
                                to use other methods different to
                                'preprocessQuantile'", class(rgSet)),
                 width = 80, prefix = " ", initial = ""
    ))
  }
  
  if (is(rgSet, "MethylSet") &&
      (processMethod[1] == "preprocessQuantile")) {
    message(strwrap("[estimateCellCounts2] The function will assume that
                            no preprocessing has been performed. Using
                            'preprocessQuantile' in prenormalized data is
                            experimental and it should only be run under the
                            user responsibility",
                    width = 80, prefix = " ",
                    initial = ""
    ))
  }
  
  # this line of code throws an error
  # fix by specifying options
  # referencePlatform <- match.arg(referencePlatform)
  referencePlatform <- match.arg(referencePlatform, 
                                 choices = c("IlluminaHumanMethylation450k", 
                                             "IlluminaHumanMethylationEPIC", 
                                             "IlluminaHumanMethylation27k"))
  rgPlatform <- sub(
    "IlluminaHumanMethylation", "",
    annotation(rgSet)[which(names(annotation(rgSet)) == "array")]
  )
  
  platform <- sub("IlluminaHumanMethylation", "", referencePlatform)
  
  # check which cell type is being deconvoluted
  # this suggestion does not make sense as raw cord blood betas may be 
  # processed with cordbloodcombined for deconvolution
  # i would ignore this message entirely if deconvoluting cord blood
  if ((compositeCellType == "CordBlood" |
       compositeCellType == "CordBloodNorway" |
       compositeCellType == "CordBloodCanada" |
       compositeCellType == "CordTissueAndBlood")) {
    message(strwrap("[estimateCellCounts2] Consider using CordBloodCombined
                        for cord blood estimation.\n",
                    width = 80, prefix = " ", initial = ""
    ))
  }
  
  if ((compositeCellType == "Blood") &&
      (referencePlatform == "IlluminaHumanMethylation450k")) {
    message(strwrap("[estimateCellCounts2] Consider using referencePlatform
                        IlluminaHumanMethylationEPIC for blood estimation,
                        using the IDOLOptimizedCpGs450klegacy\n",
                    width = 80, prefix = " ", initial = ""
    ))
  }
  
  if ((compositeCellType == "BloodExtended") &&
      (referencePlatform == "IlluminaHumanMethylation450k")) {
    message(strwrap("[estimateCellCounts2] Platform 450k is not available
                        defaulting to platform EPIC for your estimation.\n",
                    width = 80, prefix = " ", initial = ""
    ))
    referencePlatform <- "IlluminaHumanMethylationEPIC"
    platform <- sub("IlluminaHumanMethylation", "", referencePlatform)
  }
  
  if ((compositeCellType == "CordBloodCombined") &&
      (referencePlatform == "IlluminaHumanMethylationEPIC")) {
    message(strwrap("[estimateCellCounts2] Platform EPIC is not available
                        defaulting to platform 450k for your estimation.\n",
                    width = 80, prefix = " ", initial = ""
    ))
    referencePlatform <- "IlluminaHumanMethylation450k"
    platform <- sub("IlluminaHumanMethylation", "", referencePlatform)
  }
  
  if ((compositeCellType == "CordBlood" |
       compositeCellType == "CordBloodCombined") && (!"nRBC" %in% cellTypes)) {
    message(strwrap("[estimateCellCounts2] Consider including 'nRBC' in
                        argument 'cellTypes' for cord blood estimation.\n",
                    width = 80, prefix = " ", initial = ""
    ))
  }
  
  if ((compositeCellType == "Blood" |
       compositeCellType == "BloodExtended") &&
      (referencePlatform == "IlluminaHumanMethylationEPIC") &&
      ("Gran" %in% cellTypes)) {
    message(strwrap("[estimateCellCounts2] Replace 'Gran' for 'Neu' in
                        argument 'cellTypes' for EPIC blood estimation.\n",
                    width = 80, prefix = " ", initial = ""
    ))
  }
  
  if ((compositeCellType == "BloodExtended") &&
      (referencePlatform == "IlluminaHumanMethylationEPIC") &&
      ("CD4T" %in% cellTypes)) {
    message(strwrap("[estimateCellCounts2] Replace 'CD4T' for 'CD4nv',
                        'CD4mem' and 'Treg' in argument 'cellTypes' for EPIC
                        blood extended estimation.\n",
                    width = 80, prefix = " ", initial = ""
    ))
  }
  
  if ((compositeCellType == "BloodExtended") &&
      (referencePlatform == "IlluminaHumanMethylationEPIC") &&
      ("CD8T" %in% cellTypes)) {
    message(strwrap("[estimateCellCounts2] Replace 'CD8T' for 'CD8nv' and
                        'CD8mem' in argument 'cellTypes' for EPIC blood
                        extended estimation.\n",
                    width = 80, prefix = " ", initial = ""
    ))
  }
  
  if ((compositeCellType == "BloodExtended") &&
      (referencePlatform == "IlluminaHumanMethylationEPIC") &&
      ("Bcell" %in% cellTypes)) {
    message(strwrap("[estimateCellCounts2] Replace 'Bcell' for 'Bnv' and
                        'Bmem' in argument 'cellTypes' for EPIC blood
                        extended estimation.\n",
                    width = 80, prefix = " ", initial = ""
    ))
  }
  
  if ((compositeCellType != "Blood" &&
       compositeCellType != "BloodExtended") &&
      ("Neu" %in% cellTypes)) {
    message(strwrap("[estimateCellCounts2] Check whether 'Gran' or 'Neu' is
                        present in your reference and adjust argument
                        'cellTypes' for your estimation.\n",
                    width = 80, prefix = " ", initial = ""
    ))
  }
  
  if ((compositeCellType == "Blood" |
       compositeCellType == "BloodExtended") &&
      ("Gran" %in% cellTypes)) {
    message(strwrap("[estimateCellCounts2] Check whether 'Gran' or 'Neu' is
                        present in your reference and adjust argument
                        'cellTypes' for your estimation.\n",
                    width = 80, prefix = " ", initial = ""
    ))
  }
  
  if ((compositeCellType != "BloodExtended") &&
      ("Eos" %in% cellTypes)) {
    message(strwrap("[estimateCellCounts2] Check whether 'Eos' is
                        present in your reference and adjust argument
                        'cellTypes' for your estimation.\n",
                    width = 80, prefix = " ", initial = ""
    ))
  }
  
  if ((compositeCellType != "BloodExtended") &&
      ("Bas" %in% cellTypes)) {
    message(strwrap("[estimateCellCounts2] Check whether 'Bas' is
                        present in your reference and adjust argument
                        'cellTypes' for your estimation.\n",
                    width = 80, prefix = " ", initial = ""
    ))
  }
  
  if (!is.null(cellcounts) && length(cellcounts) != dim(rgSet)[2]) {
    stop(strwrap(sprintf("length of cellcounts is '%s', but needs to be equal
                        to the number of samples '%s'",
                         length(cellcounts), dim(rgSet)[2],
                         width = 80,
                         prefix = " ", initial = ""
    )))
  }
  
  # referencePkg should only be specified if referenceset is null
  # otherwise referncePkg should be the same as referenceset
  # otherwise raw cord blood is assumed to be cordbloodcombined
  # which doesnt make a lot of sense in how we specify parameters
  # ive modified the code to reflect this
  # referencePkg <- sprintf("FlowSorted.%s.%s", compositeCellType, platform)
  if(is.null(referenceset)){
    referencePkg <- sprintf("FlowSorted.%s.%s", compositeCellType, platform)
  } else {
    referencePkg <- referenceset
  }
  
  # setting subverbose 
  subverbose <- max(as.integer(verbose) - 1L, 0L)
  
  # retrieve reference set 
  if (!is.null(referenceset)) {
    referenceRGset <- get(referenceset)
    if (is(referenceRGset, "RGChannelSetExtended")) {
      referenceRGset <- as(referenceRGset, "RGChannelSet")
    }
    if (!is(rgSet, "RGChannelSet")) {
      referenceRGset <- preprocessRaw(referenceRGset)
    }
  } else {
    if (!require(referencePkg, character.only = TRUE) &&
        referencePkg != "FlowSorted.BloodExtended.EPIC") {
      stop(strwrap(sprintf(
        "Could not find reference data package for
                                compositeCellType '%s' and referencePlatform
                                '%s' (inferred package name is '%s')",
        compositeCellType, platform, referencePkg
      ),
      width = 80, prefix = " ", initial = ""
      ))
    }
    if (!require(referencePkg, character.only = TRUE) &&
        referencePkg == "FlowSorted.BloodExtended.EPIC") {
      stop(strwrap(sprintf(
        "Could not find reference data package for
                                compositeCellType '%s' and referencePlatform
                                '%s' (inferred package name is '%s'),
                                please contact
                                Technology.Transfer@dartmouth.edu",
        compositeCellType, platform, referencePkg
      ),
      width = 80, prefix = " ", initial = ""
      ))
    }
    if ((referencePkg != "FlowSorted.Blood.EPIC") &&
        (referencePkg != "FlowSorted.CordBloodCombined.450k")) {
      referenceRGset <- get(referencePkg)
    } else if ((referencePkg == "FlowSorted.Blood.EPIC") |
               (referencePkg == "FlowSorted.CordBloodCombined.450k")) {
      referenceRGset <- libraryDataGet(referencePkg)
    }
    if (!is(rgSet, "RGChannelSet")) {
      referenceRGset <- preprocessRaw(referenceRGset)
    }
  }
  
  # check that refernce set platform and data set platform are the same
  if (rgPlatform != platform) {
    rgSet <- convertArray(rgSet,
                          outType = referencePlatform,
                          verbose = TRUE
    )
  }
  
  # check that reference data set has a column called CellType
  if (!"CellType" %in% names(colData(referenceRGset))) {
    stop(strwrap(sprintf("the reference sorted dataset (in this case '%s')
                            needs to have a phenoData column called
                            'CellType'"), names(referencePkg),
                 width = 80, prefix = " ", initial = ""
    ))
  }
  
  # check that sample names are not the same between reference and user data
  if (sum(colnames(rgSet) %in% colnames(referenceRGset)) > 0) {
    stop(strwrap("the sample/column names in the user set must not be in
                    the reference data ",
                 width = 80, prefix = " ",
                 initial = ""
    ))
  }
  
  if (!all(cellTypes %in% referenceRGset$CellType)) {
    stop(strwrap(sprintf(
      "all elements of argument 'cellTypes' needs to be
                            part of the reference phenoData columns 'CellType'
                            (containg the following elements: '%s')",
      paste(unique(referenceRGset$cellType),
            collapse = "', '"
      )
    ),
    width = 80,
    prefix = " ", initial = ""
    ))
  }
  
  if (length(unique(cellTypes)) < 2) {
    stop("At least 2 cell types must be provided.")
  }
  
  if ((processMethod == "auto") &&
      (compositeCellType %in% c("Blood", "DLPFC"))) {
    processMethod <- "preprocessQuantile"
  }
  
  if ((processMethod == "auto") &&
      (!compositeCellType %in% c("Blood", "DLPFC")) &&
      (is(rgSet, "RGChannelSet"))) {
    processMethod <- "preprocessNoob"
  }
  
  processMethod <- get(processMethod)
  
  if ((probeSelect == "auto") &&
      (compositeCellType %in% c(
        "CordBloodCombined", "CordBlood",
        "CordBloodNorway", "CordTissueAndBlood"
      ))) {
    probeSelect <- "any"
  }
  
  if ((probeSelect == "auto") &&
      (!compositeCellType %in% c(
        "CordBloodCombined", "CordBlood",
        "CordBloodNorway", "CordTissueAndBlood"
      ))) {
    probeSelect <- "both"
  }
  
  if ((probeSelect == "IDOL") &&
      # this should also include cord blood???
      (compositeCellType %in% c(
        "CordBloodCombined", "CordBlood",
        "CordBloodNorway", "CordTissueAndBlood"
      ))) {
    CustomCpGs <- FlowSorted.Blood.EPIC::IDOLOptimizedCpGsCordBlood
  }
  
  if ((probeSelect == "IDOL") &&
      (compositeCellType == "Blood")) {
    if ((rgPlatform == "450k")) {
      CustomCpGs <- FlowSorted.Blood.EPIC::IDOLOptimizedCpGs450klegacy
    } else {
      CustomCpGs <- FlowSorted.Blood.EPIC::IDOLOptimizedCpGs
    }
  }
  
  if (verbose) {
    message(strwrap("[estimateCellCounts2] Combining user data with
                        reference (flow sorted) data.\n",
                    width = 80,
                    prefix = " ", initial = ""
    ))
  }
  
  newpd <- DataFrame(
    sampleNames = c(
      colnames(rgSet),
      colnames(referenceRGset)
    ),
    studyIndex = rep(c("user", "reference"),
                     times = c(
                       ncol(rgSet),
                       ncol(referenceRGset)
                     )
    )
  )
  
  referenceRGset$CellType <- as.character(referenceRGset$CellType)
  
  # check if column names exist
  if (is.null(rgSet$CellType)) {
    rgSet$CellType <- rep("NA", dim(rgSet)[2])
  }
  if (is.null(rgSet$Age)) {
    rgSet$Age <- rep("NA", dim(rgSet)[2])
  }
  if (is.null(rgSet$Sex)) {
    rgSet$Sex <- rep("NA", dim(rgSet)[2])
  }
  if (is.null(referenceRGset$Sex)) {
    referenceRGset$Sex <- rep("NA", dim(referenceRGset)[2])
  } else {
    referenceRGset$Sex <- as.character(referenceRGset$Sex)
  }
  if (is.null(referenceRGset$Age)) {
    referenceRGset$Age <- rep("NA", dim(referenceRGset)[2])
  } else {
    try(referenceRGset$Age <- as.numeric(referenceRGset$Age))
  }
  
  commoncolumn <- intersect(
    names(colData(rgSet)),
    names(colData(referenceRGset))
  )
  
  # this chunk produces the error message 
  # "In asMethod(object) : NAs introduced by coercion"
  # should as() function come from BiocGenerics or base R?
  restry <- try(
    {
      colData(rgSet)[commoncolumn] <- mapply(
        FUN = as,
        colData(rgSet)[commoncolumn],
        vapply(colData(referenceRGset)[commoncolumn],
               class,
               FUN.VALUE = character(1)
        ),
        SIMPLIFY = FALSE
      )
    },
    silent = F
  )
  
  # this chunk should catch and rectify the erorr but it is not maybe??
  if ("try-error" %in% class(restry)) {
    commoncolumn <- c("CellType", "Sex", "Age")
    colData(rgSet)[commoncolumn] <- mapply(
      FUN = as,
      colData(rgSet)[commoncolumn],
      vapply(colData(referenceRGset)[commoncolumn],
             class,
             FUN.VALUE = character(1)
      ),
      SIMPLIFY = FALSE
    )
  } else {
    colData(rgSet)[commoncolumn] <- mapply(
      FUN = as,
      colData(rgSet)[commoncolumn],
      vapply(colData(referenceRGset)[commoncolumn],
             class,
             FUN.VALUE = character(1)
      ),
      SIMPLIFY = FALSE
    )
  }
  
  rm(restry)
  colData(referenceRGset) <- colData(referenceRGset)[commoncolumn]
  colData(rgSet) <- colData(rgSet)[commoncolumn]
  referencePd <- colData(referenceRGset)
  combinedRGset <- combineArrays(rgSet, referenceRGset,
                                 outType = referencePlatform
  )
  
  colData(combinedRGset) <- newpd
  colnames(combinedRGset) <- newpd$sampleNames
  rm(referenceRGset)
  
  if (verbose) {
    message(strwrap("[estimateCellCounts2] Processing user and reference
                        data together.\n",
                    width = 80, prefix = " ",
                    initial = ""
    ))
  }
  
  
  if (compositeCellType == "CordBlood") {
    if (!is(combinedRGset, "RGChannelSet")) {
      combinedRGset@preprocessMethod["rg.norm"] <-
        "Raw (no normalization or bg correction)"
    }
    combinedMset <- processMethod(combinedRGset, verbose = subverbose)
    rm(combinedRGset)
    gc()
    compTable <- get(paste0(referencePkg, ".compTable"))
    combinedMset <- combinedMset[which(rownames(combinedMset) %in%
                                         rownames(compTable)), ]
  } else {
    if (!is(combinedRGset, "RGChannelSet")) {
      combinedRGset@preprocessMethod["rg.norm"] <-
        "Raw (no normalization or bg correction)"
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

  if (probeSelect != "IDOL") {
    if (verbose) {
      message(strwrap("[estimateCellCounts2] Picking probes for
                            composition estimation.\n",
                      width = 80,
                      prefix = " ", initial = ""
      ))
    }
    compData <- pickCompProbes(referenceMset,
                               cellTypes = cellTypes,
                               compositeCellType = compositeCellType,
                               probeSelect = probeSelect
    )
    coefs <- compData$coefEsts
    if (verbose) {
      message(strwrap("[estimateCellCounts2] Estimating  proportion
                            composition (prop), if you provide cellcounts
                            those will be provided as counts in the
                            composition estimation.\n",
                      width = 80,
                      prefix = " ", initial = ""
      ))
    }
    prop <- projectCellType_CP(getBeta(mSet)[rownames(coefs), ], coefs,
                               lessThanOne = lessThanOne
    )
    prop <- round(prop, 4)
    rownames(prop) <- colnames(rgSet)
    counts <- round(prop * cellcounts, 0)
    if (meanPlot) {
      smeans <- compData$sampleMeans
      smeans <- smeans[order(names(smeans))]
      sampleMeans <- c(
        colMeans(getBeta(mSet)[rownames(coefs), ]),
        smeans
      )
      sampleColors <- c(rep(1, ncol(mSet)), 1 +
                          as.numeric(factor(names(smeans))))
      plot(sampleMeans, pch = 21, bg = sampleColors)
      legend("bottomleft", c("blood", levels(factor(names(smeans)))),
             col = seq_len(7), pch = 15
      )
    }
    if (returnAll) {
      list(
        prop = prop, counts = counts, compTable = compData$compTable,
        normalizedData = mSet
      )
    } else {
      list(prop = prop, counts = counts)
    }
  } else {
    if (verbose) {
      message(strwrap("[estimateCellCounts2] Using IDOL L-DMR probes for
                            composition estimation.\n",
                      width = 80,
                      prefix = " ", initial = ""
      ))
    }
    p <- getBeta(referenceMset)
    pd <- as.data.frame(colData(referenceMset))
    rm(referenceMset)
    if (!is.null(cellTypes)) {
      if (!all(cellTypes %in% pd$CellType)) {
        stop(strwrap("elements of argument 'cellTypes' are not part of
                            'referenceMset$CellType'",
                     width = 80,
                     prefix = " ", initial = ""
        ))
      }
      keep <- which(pd$CellType %in% cellTypes)
      pd <- pd[keep, ]
      p <- p[, keep]
    }
    pd$CellType <- factor(pd$CellType, levels = cellTypes)
    ffComp <- rowFtests(p, pd$CellType)
    tIndexes <- split(seq(along = pd$CellType), pd$CellType)
    prof <- vapply(tIndexes, function(i) rowMeans(p[, i]),
                   FUN.VALUE = numeric(dim(p)[1])
    )
    r <- rowRanges(p)
    compTable <- cbind(ffComp, prof, r, abs(r[, 1] - r[, 2]))
    names(compTable)[1] <- "Fstat"
    names(compTable)[c(-2, -1, 0) + ncol(compTable)] <- c(
      "low",
      "high", "range"
    )
    tstatList <- lapply(tIndexes, function(i) {
      x <- rep(0, ncol(p))
      x[i] <- 1
      return(rowttests(p, factor(x)))
    })
    trainingProbes <- CustomCpGs
    trainingProbes <- trainingProbes[trainingProbes %in% rownames(p)]
    p <- p[trainingProbes, ]
    pMeans <- colMeans(p)
    names(pMeans) <- pd$CellType
    form <- as.formula(sprintf("y ~ %s - 1", paste(levels(pd$CellType),
                                                   collapse = "+"
    )))
    phenoDF <- as.data.frame(model.matrix(~ pd$CellType - 1))
    colnames(phenoDF) <- sub("^pd\\$CellType", "", colnames(phenoDF))
    if (ncol(phenoDF) == 2) {
      X <- as.matrix(phenoDF)
      coefEsts <- t(solve(t(X) %*% X) %*% t(X) %*% t(p))
      coefs <- coefEsts
    } else {
      tmp <- validationCellType(Y = p, pheno = phenoDF, modelFix = form)
      coefEsts <- tmp$coefEsts
      coefs <- coefEsts
    }
    compData <- list(
      coefEsts = coefEsts, compTable = compTable,
      sampleMeans = pMeans
    )
    if (verbose) {
      message(strwrap("[estimateCellCounts2] Estimating  proportion
                            composition (prop), if you provide cellcounts
                            those will be provided as counts in the
                            composition estimation.\n",
                      width = 80,
                      prefix = " ", initial = ""
      ))
    }
    prop <- projectCellType_CP(getBeta(mSet)[rownames(coefs), ], coefs,
                               lessThanOne = lessThanOne
    )
    prop <- round(prop, 4)
    rownames(prop) <- colnames(rgSet)
    counts <- round(prop * cellcounts, 0)
    if (meanPlot) {
      smeans <- compData$sampleMeans
      smeans <- smeans[order(names(smeans))]
      sampleMeans <- c(colMeans(getBeta(mSet)[rownames(coefs), ]), smeans)
      sampleColors <- c(rep(1, ncol(mSet)), 1 +
                          as.numeric(factor(names(smeans))))
      plot(sampleMeans, pch = 21, bg = sampleColors)
      legend("bottomleft", c("blood", levels(factor(names(smeans)))),
             col = seq_len(7), pch = 15
      )
    }
    if (returnAll) {
      list(
        prop = prop, counts = counts, compTable = compTable,
        normalizedData = mSet
      )
    } else {
      list(prop = prop, counts = counts)
    }
  }
}

# These are minfi internal functions here they are called to keep the function
# above running for the alternative selection ("auto" and "both")
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
  ffComp <- rowFtests(p, pd$CellType)
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
  N <- dim(pheno)[1]
  pheno$y <- rep(0, N)
  xTest <- model.matrix(modelFix, pheno)
  sizeModel <- dim(xTest)[2]
  M <- dim(Y)[1]
  if (is.null(L.forFstat)) {
    L.forFstat <- diag(sizeModel)[-1, ] # All non-intercept coefficients
    colnames(L.forFstat) <- colnames(xTest)
    rownames(L.forFstat) <- colnames(xTest)[-1]
  }
  # Initialize various containers
  sigmaResid <- sigmaIcept <- nObserved <- nClusters <- Fstat <- rep(NA, M)
  coefEsts <- matrix(NA, M, sizeModel)
  coefVcovs <- list()
  if (verbose) {
    cat("[validationCellType] ")
  }
  for (j in seq_len(M)) { # For each CpG
    ## Remove missing methylation values
    ii <- !is.na(Y[j, ])
    nObserved[j] <- sum(ii)
    pheno$y <- Y[j, ]
    
    if (j %% round(M / 10) == 0 && verbose) {
      cat(".")
    } # Report progress
    
    try({ # Try to fit a mixed model to adjust for plate
      if (!is.null(modelBatch)) {
        fit <- try(lme(modelFix, random = modelBatch, data = pheno[ii, ]))
        OLS <- inherits(fit, "try-error")
        # If LME can't be fit, just use OLS
      } else {
        OLS <- TRUE
      }
      
      if (OLS) {
        fit <- lm(modelFix, data = pheno[ii, ])
        fitCoef <- fit$coef
        sigmaResid[j] <- summary(fit)$sigma
        sigmaIcept[j] <- 0
        nClusters[j] <- 0
      } else {
        fitCoef <- fit$coef$fixed
        sigmaResid[j] <- fit$sigma
        sigmaIcept[j] <- sqrt(getVarCov(fit)[1])
        nClusters[j] <- length(fit$coef$random[[1]])
      }
      coefEsts[j, ] <- fitCoef
      coefVcovs[[j]] <- vcov(fit)
      
      useCoef <- L.forFstat %*% fitCoef
      useV <- L.forFstat %*% coefVcovs[[j]] %*% t(L.forFstat)
      Fstat[j] <- (t(useCoef) %*% solve(useV, useCoef)) / sizeModel
    })
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
