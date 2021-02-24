
rgSet=cb_rgset
compositeCellType="CordBloodCombined"
processMethod="preprocessNoob"
probeSelect="any"
cellTypes = c("CD8T", "CD4T", "NK", "Bcell","Mono")
referencePlatform="IlluminaHumanMethylation450k"
referenceset="FlowSorted.CordBloodCombined.450k"
IDOLOptimizedCpGs = NULL
returnAll = FALSE
meanPlot = FALSE 
verbose = TRUE
lessThanOne = FALSE

rgSet <- as(rgSet, "RGChannelSet")

rgPlatform <- sub("IlluminaHumanMethylation", "", 
                  annotation(rgSet)[which(names(annotation(rgSet)) == "array")])
platform <- sub("IlluminaHumanMethylation", "", 
                referencePlatform)

if (compositeCellType == "CordBloodCombined") 
  platform <= "450k"
referencePkg <- sprintf("FlowSorted.%s.%s", compositeCellType, 
                        platform)


  subverbose <- max(as.integer(verbose) - 1L, 0L)
  
  
  if ((referencePkg != "FlowSorted.Blood.EPIC") && 
      (referencePkg != "FlowSorted.CordBloodCombined.450k")) {
    referenceRGset <- get(referencePkg)
  }
  else if (referencePkg == "FlowSorted.Blood.EPIC") {
    hub <- ExperimentHub()
    referenceRGset <- hub[["EH1136"]]
  }
  else if (referencePkg == "FlowSorted.CordBloodCombined.450k") {
    hub <- ExperimentHub()
    referenceRGset <- hub[["EH2256"]]
  }
  
  processMethod <- get(processMethod)
  
  
  if (verbose) 
    message(strwrap("[estimateCellCounts2] Combining user data with \n                        reference (flow sorted) data.\n", 
                    width = 80, prefix = " ", initial = ""))
  newpd <- DataFrame(sampleNames = c(colnames(rgSet), colnames(referenceRGset)), 
                     studyIndex = rep(c("user", "reference"), 
                                      times = c(ncol(rgSet), ncol(referenceRGset))))
  
  
  
  
  
  
  referenceRGset <- get(referenceset)
  if (!is(rgSet, "RGChannelSet")) 
    referenceRGset <- preprocessRaw(referenceRGset)
  }
  
  
  if (!is.null(referenceset)) {
    referenceRGset <- get(referenceset)
    if (!is(rgSet, "RGChannelSet")) 
      referenceRGset <- preprocessRaw(referenceRGset)
  }
  else {
    if (!require(referencePkg, character.only = TRUE)) 
      stop(strwrap(sprintf("Could not find reference data package for \n                                compositeCellType '%s' and referencePlatform \n                                '%s' (inferred package name is '%s')", 
                           compositeCellType, platform, referencePkg), width = 80, 
                   prefix = " ", initial = ""))
    if ((referencePkg != "FlowSorted.Blood.EPIC") && 
        (referencePkg != "FlowSorted.CordBloodCombined.450k")) {
      referenceRGset <- get(referencePkg)
    }
    else if (referencePkg == "FlowSorted.Blood.EPIC") {
      hub <- ExperimentHub()
      referenceRGset <- hub[["EH1136"]]
    }
    else if (referencePkg == "FlowSorted.CordBloodCombined.450k") {
      hub <- ExperimentHub()
      referenceRGset <- hub[["EH2256"]]
    }
    if (!is(rgSet, "RGChannelSet")) 
      referenceRGset <- preprocessRaw(referenceRGset)
  }
  if (rgPlatform != platform) {
    rgSet <- convertArray(rgSet, outType = referencePlatform, 
                          verbose = TRUE)
  }
  if (!"CellType" %in% names(colData(referenceRGset))) 
    stop(strwrap(sprintf("the reference sorted dataset (in this case '%s') \n                            needs to have a phenoData column called \n                            'CellType'"), 
                 names(referencePkg), width = 80, prefix = " ", 
                 initial = ""))
  if (sum(colnames(rgSet) %in% colnames(referenceRGset)) > 
      0) 
    stop(strwrap("the sample/column names in the user set must not be in \n                    the reference data ", 
                 width = 80, prefix = " ", initial = ""))
  if (!all(cellTypes %in% referenceRGset$CellType)) 
    stop(strwrap(sprintf("all elements of argument 'cellTypes' needs to be \n                            part of the reference phenoData columns 'CellType' \n                            (containg the following elements: '%s')", 
                         paste(unique(referenceRGset$cellType), collapse = "', '")), 
                 width = 80, prefix = " ", initial = ""))
  if (length(unique(cellTypes)) < 2) 
    stop("At least 2 cell types must be provided.")
  if ((processMethod == "auto") && (compositeCellType %in% 
                                    c("Blood", "DLPFC"))) 
    processMethod <- "preprocessQuantile"
  if ((processMethod == "auto") && (!compositeCellType %in% 
                                    c("Blood", "DLPFC")) && (is(rgSet, "RGChannelSet"))) 
    processMethod <- "preprocessNoob"
  processMethod <- get(processMethod)
  if ((probeSelect == "auto") && (compositeCellType %in% 
                                  c("CordBloodCombined", "CordBlood", "CordBloodNorway", 
                                    "CordTissueAndBlood"))) {
    probeSelect <- "any"
  }
  if ((probeSelect == "auto") && (!compositeCellType %in% 
                                  c("CordBloodCombined", "CordBlood", "CordBloodNorway", 
                                    "CordTissueAndBlood"))) {
    probeSelect <- "both"
  }
  if (verbose) 
    message(strwrap("[estimateCellCounts2] Combining user data with \n                        reference (flow sorted) data.\n", 
                    width = 80, prefix = " ", initial = ""))
  newpd <- DataFrame(sampleNames = c(colnames(rgSet), colnames(referenceRGset)), 
                     studyIndex = rep(c("user", "reference"), 
                                      times = c(ncol(rgSet), ncol(referenceRGset))))
  
  #copy from here
  referenceRGset$CellType <- as.character(referenceRGset$CellType)
  
  
  if (is.null(rgSet$CellType)) 
    rgSet$CellType <- rep("NA", dim(rgSet)[2])
  if (is.null(rgSet$Age)) 
    rgSet$Age <- rep("NA", dim(rgSet)[2])
  if (is.null(rgSet$Sex)) 
    rgSet$Sex <- rep("NA", dim(rgSet)[2])
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
  } else {
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
  
  if (verbose) 
    message(strwrap("[estimateCellCounts2] Processing user and reference \n                        data together.\n", 
                    width = 80, prefix = " ", initial = ""))
  
  if (compositeCellType == "CordBlood") {
    if (!is(combinedRGset, "RGChannelSet")) 
      combinedRGset@preprocessMethod["rg.norm"] <- "Raw (no normalization or bg correction)"
    combinedMset <- processMethod(combinedRGset, verbose = subverbose)
    rm(combinedRGset)
    gc()
    compTable <- get(paste0(referencePkg, ".compTable"))
    combinedMset <- combinedMset[which(rownames(combinedMset) %in% 
                                         rownames(compTable)), ]
  } else {
    if (!is(combinedRGset, "RGChannelSet")) 
      combinedRGset@preprocessMethod["rg.norm"] <- "Raw (no normalization or bg correction)"
    combinedMset <- processMethod(combinedRGset)
    rm(combinedRGset)
    gc()
  }
  
  
  
  
  referenceMset <- combinedMset[, combinedMset$studyIndex == 
                                  "reference"]
  
  colData(referenceMset) <- as(referencePd, "DataFrame")
  
  mSet <- combinedMset[, combinedMset$studyIndex == "user"]
  
  colData(mSet) <- as(colData(rgSet), "DataFrame")
  
  rm(combinedMset)
  
  
  mSet = referenceMset
  cellTypes = cellTypes
  compositeCellType = compositeCellType
  probeSelect = probeSelect
  
  
  pickCompProbes <- function(mSet, cellTypes = NULL, numProbes = 50, 
                             compositeCellType = compositeCellType, 
                             probeSelect = probeSelect) {
    p <- getBeta(mSet)
    pd <- as.data.frame(colData(mSet))
    if(!is.null(cellTypes)) {
      if(!all(cellTypes %in% pd$CellType))
        stop(strwrap("elements of argument 'cellTypes' are not part of 
                        'mSet$CellType'", width = 80, prefix = " ", 
                     initial = ""))
      keep <- which(pd$CellType %in% cellTypes)
      pd <- pd[keep,]
      p <- p[,keep]
    }
    ## make cell type a factor 
    pd$CellType <- factor(pd$CellType, levels = cellTypes)
    ffComp <- rowFtests(p, pd$CellType)
    tIndexes <- split(seq(along=pd$CellType), pd$CellType)
    prof <- vapply(tIndexes, function(i) rowMeans(p[,i]), 
                   FUN.VALUE=numeric(dim(p)[1]))
    r <- rowRanges(p)
    compTable <- cbind(ffComp, prof, r, abs(r[,1] - r[,2]))
    names(compTable)[1] <- "Fstat"
    names(compTable)[c(-2,-1,0) + ncol(compTable)] <- c("low", "high", "range") 
    tstatList <- lapply(tIndexes, function(i) {
      x <- rep(0,ncol(p))
      x[i] <- 1
      return(rowttests(p, factor(x)))
    })
    if (probeSelect == "any"){
      probeList <- lapply(tstatList, function(x) {
        y <- x[x[,"p.value"] < 1e-8,]
        yAny <- y[order(abs(y[,"dm"]), decreasing=TRUE),]      
        c(rownames(yAny)[seq_len(numProbes*2)])
      })
    } else {
      probeList <- lapply(tstatList, function(x) {
        y <- x[x[,"p.value"] < 1e-8,]
        yUp <- y[order(y[,"dm"], decreasing=TRUE),]
        yDown <- y[order(y[,"dm"], decreasing=FALSE),]
        c(rownames(yUp)[seq_len(numProbes)], 
          rownames(yDown)[seq_len(numProbes)])
      })
    }
    
    #####################################
    #####HERE IS WERE WE LOSE PROBES#####
    #####################################
    trainingProbes <- unique(unlist(probeList))
    p <- p[trainingProbes,]
    
    pMeans <- colMeans(p)
    names(pMeans) <- pd$CellType
    form <- as.formula(sprintf("y ~ %s - 1", paste(levels(pd$CellType), 
                                                   collapse="+")))
    phenoDF <- as.data.frame(model.matrix(~pd$CellType-1))
    colnames(phenoDF) <- sub("^pd\\$CellType", "", colnames(phenoDF))
    if(ncol(phenoDF) == 2) { # two group solution
      X <- as.matrix(phenoDF)
      coefEsts <- t(solve(t(X) %*% X) %*% t(X) %*% t(p))
    } else { # > 2 group solution

      tmp <- validationCellType(Y = p, pheno = phenoDF, modelFix = form)
      coefEsts <- tmp$coefEsts
    }
    
    out <- list(coefEsts = coefEsts, compTable = compTable,
                sampleMeans = pMeans)
    return(out)
  }
  
  
  
  
  
  Y = p
  pheno = phenoDF
  modelFix = form
  modelBatch=NULL
  L.forFstat = NULL 
  verbose = FALSE
  
  validationCellType <- function(Y, pheno, modelFix, modelBatch=NULL,
                                 L.forFstat = NULL, verbose = FALSE){
    N <- dim(pheno)[1]
    pheno$y <- rep(0, N)
    xTest <- model.matrix(modelFix, pheno)
    sizeModel <- dim(xTest)[2]
    M <- dim(Y)[1]
    if(is.null(L.forFstat)) {
      L.forFstat <- diag(sizeModel)[-1,] # All non-intercept coefficients
      colnames(L.forFstat) <- colnames(xTest) 
      rownames(L.forFstat) <- colnames(xTest)[-1] 
    }
    # Initialize various containers
    sigmaResid <- sigmaIcept <- nObserved <- nClusters <- Fstat <- rep(NA, M)
    coefEsts <- matrix(NA, M, sizeModel)
    coefVcovs <- list()
    if(verbose)
      cat("[validationCellType] ")
    for(j in seq_len(M)) { # For each CpG
      ## Remove missing methylation values
      ii <- !is.na(Y[j,])
      nObserved[j] <- sum(ii)
      pheno$y <- Y[j,]
      
      if(j%%round(M/10)==0 && verbose)
        cat(".") # Report progress
      
      try({ # Try to fit a mixed model to adjust for plate
        if(!is.null(modelBatch)) {
          fit <- try(lme(modelFix, random=modelBatch, data=pheno[ii,]))
          OLS <- inherits(fit,"try-error") 
          # If LME can't be fit, just use OLS
        } else
          OLS <- TRUE
        
        if(OLS) {
          fit <- lm(modelFix, data=pheno[ii,])
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
        coefEsts[j,] <- fitCoef
        coefVcovs[[j]] <- vcov(fit)
        
        useCoef <- L.forFstat %*% fitCoef
        useV <- L.forFstat %*% coefVcovs[[j]] %*% t(L.forFstat)
        Fstat[j] <- (t(useCoef) %*% solve(useV, useCoef))/sizeModel
      })
    }
    if(verbose)
      cat(" done\n")
    ## Name the rows so that they can be easily matched to the target data set
    rownames(coefEsts) <- rownames(Y)
    colnames(coefEsts) <- names(fitCoef)
    degFree <- nObserved - nClusters - sizeModel + 1
    
    ## Get P values corresponding to F statistics
    Pval <- 1-pf(Fstat, sizeModel, degFree)
    
    out <- list(coefEsts=coefEsts, coefVcovs=coefVcovs, modelFix=modelFix, 
                modelBatch=modelBatch,
                sigmaIcept=sigmaIcept, sigmaResid=sigmaResid, 
                L.forFstat=L.forFstat, Pval=Pval,
                orderFstat=order(-Fstat), Fstat=Fstat, nClusters=nClusters, 
                nObserved=nObserved,
                degFree=degFree)
    
    out
  }
  
  projectCellType <- function(Y, coefCellType, contrastCellType=NULL, 
                              nonnegative=TRUE, lessThanOne=lessThanOne){ 
    if(is.null(contrastCellType))
      Xmat <- coefCellType
    else
      Xmat <- tcrossprod(coefCellType, contrastCellType) 
    nCol <- dim(Xmat)[2]
    if(nCol == 2) {
      Dmat <- crossprod(Xmat)
      mixCoef <- t(apply(Y, 2, function(x) {solve(Dmat, crossprod(Xmat, x))}))
      colnames(mixCoef) <- colnames(Xmat)
      return(mixCoef)
    } else {
      nSubj <- dim(Y)[2]
      mixCoef <- matrix(0, nSubj, nCol)
      rownames(mixCoef) <- colnames(Y)
      colnames(mixCoef) <- colnames(Xmat)
      if(nonnegative){
        if(lessThanOne) {
          Amat <- cbind(rep(-1, nCol), diag(nCol))
          b0vec <- c(-1, rep(0, nCol))
        } else {
          Amat <- diag(nCol)
          b0vec <- rep(0, nCol)
        }
        for(i in seq_len(nSubj)) {
          obs <- which(!is.na(Y[,i])) 
          Dmat <- crossprod(Xmat[obs,])
          mixCoef[i,] <- solve.QP(Dmat, crossprod(Xmat[obs,], Y[obs,i]), 
                                  Amat, b0vec)$sol
        }
      } else {
        for(i in seq_len(nSubj)) {
          obs <- which(!is.na(Y[,i])) 
          Dmat <- crossprod(Xmat[obs,])
          mixCoef[i,] <- solve(Dmat, t(Xmat[obs,]) %*% Y[obs,i])
        }
      }
      return(mixCoef)
    }
  }
  
  
  
  
  if (probeSelect != "IDOL") {
    if (verbose) 
      message(strwrap("[estimateCellCounts2] Picking probes for \n                            composition estimation.\n", 
                      width = 80, prefix = " ", initial = ""))
    compData <- pickCompProbes(referenceMset, cellTypes = cellTypes, 
                               compositeCellType = compositeCellType, probeSelect = probeSelect)
    coefs <- compData$coefEsts
    
    if (verbose) 
      message("[estimateCellCounts2] Estimating composition.\n")
    counts <- projectCellType(getBeta(mSet)[rownames(coefs), 
                                            ], coefs, lessThanOne = lessThanOne)
    rownames(counts) <- colnames(rgSet)
    
    if (meanPlot) {
      smeans <- compData$sampleMeans
      smeans <- smeans[order(names(smeans))]
      sampleMeans <- c(colMeans(minfi::getBeta(mSet)[rownames(coefs), 
                                                     ]), smeans)
      sampleColors <- c(rep(1, ncol(mSet)), 1 + as.numeric(factor(names(smeans))))
      plot(sampleMeans, pch = 21, bg = sampleColors)
      legend("bottomleft", c("blood", levels(factor(names(smeans)))), 
             col = seq_len(7), pch = 15)
    }
    if (returnAll) {
      list(counts = counts, compTable = compData$compTable, 
           normalizedData = mSet)
    }
    else {
      list(counts = counts)
    }
  } else {
    if (verbose) 
      message(strwrap("[estimateCellCounts2] Using IDOL L-DMR probes for \n                            composition estimation.\n", 
                      width = 80, prefix = " ", initial = ""))
    p <- getBeta(referenceMset)
    pd <- as.data.frame(colData(referenceMset))
    rm(referenceMset)
    if (!is.null(cellTypes)) {
      if (!all(cellTypes %in% pd$CellType)) 
        stop(strwrap("elements of argument 'cellTypes' is not part of \n                            'referenceMset$CellType'", 
                     width = 80, prefix = " ", initial = ""))
      keep <- which(pd$CellType %in% cellTypes)
      pd <- pd[keep, ]
      p <- p[, keep]
    }
    
    pd$CellType <- factor(pd$CellType, levels = cellTypes)
    ffComp <- rowFtests(p, pd$CellType)
    tIndexes <- split(seq(along = pd$CellType), pd$CellType)
    prof <- vapply(tIndexes, function(i) rowMeans(p[, i]), 
                   FUN.VALUE = numeric(dim(p)[1]))
    r <- rowRanges(p)
    compTable <- cbind(ffComp, prof, r, abs(r[, 1] - r[, 
                                                       2]))
    names(compTable)[1] <- "Fstat"
    names(compTable)[c(-2, -1, 0) + ncol(compTable)] <- c("low", 
                                                          "high", "range")
    tstatList <- lapply(tIndexes, function(i) {
      x <- rep(0, ncol(p))
      x[i] <- 1
      return(rowttests(p, factor(x)))
    })
    
    trainingProbes <- IDOLOptimizedCpGs
    trainingProbes <- trainingProbes[trainingProbes %in% 
                                       rownames(p)]
    p <- p[trainingProbes, ]
    pMeans <- colMeans(p)
    names(pMeans) <- pd$CellType
    form <- as.formula(sprintf("y ~ %s - 1", paste(levels(pd$CellType), 
                                                   collapse = "+")))
    phenoDF <- as.data.frame(model.matrix(~pd$CellType - 
                                            1))
    colnames(phenoDF) <- sub("^pd\\$CellType", "", 
                             colnames(phenoDF))
    if (ncol(phenoDF) == 2) {
      X <- as.matrix(phenoDF)
      coefEsts <- t(solve(t(X) %*% X) %*% t(X) %*% t(p))
      coefs <- coefEsts
    }
    else {
      tmp <- validationCellType(Y = p, pheno = phenoDF, 
                                modelFix = form)
      coefEsts <- tmp$coefEsts
      coefs <- coefEsts
    }
    compData <- list(coefEsts = coefEsts, compTable = compTable, 
                     sampleMeans = pMeans)
    if (verbose) 
      message("[estimateCellCounts2] Estimating composition.\n")
    counts <- projectCellType(getBeta(mSet)[rownames(coefs), 
                                            ], coefs, lessThanOne = lessThanOne)
    
    rownames(counts) <- colnames(rgSet)
    
    
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
      list(counts = counts, compTable = compTable, normalizedData = mSet)
    }
    else {
      list(counts = counts)
    }
  }
}