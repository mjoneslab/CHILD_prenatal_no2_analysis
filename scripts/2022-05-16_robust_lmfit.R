# 2022-05-16
# Robust lmFit
# SL

# This script modifies lmFit to output whether robust regressions 
# have converged or not. This requires modification of the internal mrlm
# function as well as lmFit. The new modified functions have been termed
# mrlm2 and lmFit2. Load these functions and then lmFit2 in place of lmFit.
# Remove CpGs that did not converge (if >= 100 consider increasing iterations
# up to 100 and/or modifying the model). Then proceed as normal in analysis! 

# Note: You must modify output terms based on the model. Current code supports
# output from three-way interaction model. Can add/remove terms as needed.

# modified mrlm function to output converged
mrlm2 <- function (M, design = NULL, ndups = 1, spacing = 1, weights = NULL, 
                   ...) 
{
  if (!requireNamespace("MASS", quietly = TRUE)) 
    stop("MASS package required but is not installed (or can't be loaded)")
  M <- as.matrix(M)
  narrays <- ncol(M)
  if (is.null(design)) 
    design <- matrix(1, narrays, 1)
  design <- as.matrix(design)
  coef.names <- colnames(design)
  nbeta <- ncol(design)
  if (!is.null(weights)) {
    weights <- asMatrixWeights(weights, dim(M))
    weights[weights <= 0] <- NA
    M[!is.finite(weights)] <- NA
  }
  if (ndups > 1) {
    M <- unwrapdups(M, ndups = ndups, spacing = spacing)
    design <- design %x% rep_len(1, ndups)
    if (!is.null(weights)) 
      weights <- unwrapdups(weights, ndups = ndups, spacing = spacing)
  }
  ngenes <- nrow(M)
  stdev.unscaled <- beta <- matrix(NA, ngenes, nbeta, dimnames = list(rownames(M), 
                                                                      coef.names))
  sigma <- rep_len(NA_real_, ngenes)
  df.residual <- rep_len(0, ngenes)
  
  # specify length of new vectors specific to analysis
  converged <- rep_len(0, ngenes)
  sig_ap <- rep_len(0, ngenes)
  # sig_ap_ses <- rep_len(0, ngenes)
  # sig_ap_pc1 <- rep_len(0, ngenes)
  # sig_ses_pc1 <- rep_len(0, ngenes)
  # sig_ap_ses_pc1 <- rep_len(0, ngenes)
  # 
  for (i in 1:ngenes) {
    y <- as.vector(M[i, ])
    obs <- is.finite(y)
    X <- design[obs, , drop = FALSE]
    y <- y[obs]
    if (is.null(weights)) 
      w <- rep_len(1, length(y))
    else w <- as.vector(weights[i, obs])
    if (length(y) > nbeta) {
      
      # specify M-estimator method within mrlm function
      # as not possible to pass two "method" parameters to lmFit
      # for clarity, specify Huber estimator and iterations and lmFit
      out <- MASS::rlm(x = X, y = y, weights = w, method = "M", ...)
      beta[i, ] <- coef(out)
      stdev.unscaled[i, ] <- sqrt(diag(chol2inv(out$qr$qr)))
      converged[i] <- out$converged
      df.residual[i] <- length(y) - out$rank
      
      # also need to specify p-values for each of the terms
      sig_ap[i] <- wald.test(Sigma = vcov(out), b = coef(out), Terms = 2)$result$chi2["P"]
      # sig_ap_ses[i] <- wald.test(Sigma = vcov(out), b = coef(out), Terms = 12)$result$chi2["P"]
      # sig_ap_pc1[i] <- wald.test(Sigma = vcov(out), b = coef(out), Terms = 13)$result$chi2["P"]
      # sig_ses_pc1[i] <- wald.test(Sigma = vcov(out), b = coef(out), Terms = 14)$result$chi2["P"]
      # sig_ap_ses_pc1[i] <- wald.test(Sigma = vcov(out), b = coef(out), Terms = 15)$result$chi2["P"]
      
      if (df.residual[i] > 0) 
        sigma[i] <- out$s
    }
  }
  
  QR <- qr(design)
  cov.coef <- chol2inv(QR$qr, size = QR$rank)
  est <- QR$pivot[1:QR$rank]
  dimnames(cov.coef) <- list(coef.names[est], coef.names[est])
  list(coefficients = beta, stdev.unscaled = stdev.unscaled, 
       sigma = sigma, df.residual = df.residual, cov.coefficients = cov.coef, 
       pivot = QR$pivot, rank = QR$rank, converged2 = converged,
       wald_sig_ap = sig_ap
       # wald_sig_ap_pc1 = sig_ap_pc1,
       # wald_sig_ap_ses = sig_ap_ses, wald_sig_ses_pc1 = sig_ses_pc1,
       # wald_sig_ap_ses_pc1 =sig_ap_ses_pc1
  )
}


# modified lmFit function to call mrlm2
lmFit2 <- function (object, design = NULL, ndups = NULL, spacing = NULL, 
                    block = NULL, correlation, weights = NULL, method = "ls", 
                    ...) 
{
  if (inherits(object, "data.frame")) {
    ColumnIsNumeric <- vapply(object, is.numeric, FUN.VALUE = TRUE)
    if (all(ColumnIsNumeric)) 
      y <- list(exprs = as.matrix(object))
    else {
      WhichNotNumeric <- which(!ColumnIsNumeric)
      if (identical(sum(WhichNotNumeric), 1L) && length(ColumnIsNumeric) > 
          1L) {
        y <- list()
        y$exprs <- as.matrix(object[, -1, drop = FALSE])
        y$probes <- object[, 1, drop = FALSE]
        message("Converting data.frame to matrix, treating first column as gene IDs.")
      }
      else {
        stop("Expression object should be numeric, instead it is a data.frame with ", 
             length(WhichNotNumeric), " non-numeric columns")
      }
    }
  }
  else {
    y <- getEAWP(object)
  }
  if (!nrow(y$exprs)) 
    stop("expression matrix has zero rows")
  if (is.null(design)) 
    design <- y$design
  if (is.null(design)) 
    design <- matrix(1, ncol(y$exprs), 1)
  else {
    design <- as.matrix(design)
    if (mode(design) != "numeric") 
      stop("design must be a numeric matrix")
    if (nrow(design) != ncol(y$exprs)) 
      stop("row dimension of design doesn't match column dimension of data object")
  }
  ne <- nonEstimable(design)
  if (!is.null(ne)) 
    cat("Coefficients not estimable:", paste(ne, collapse = " "), 
        "\n")
  if (is.null(ndups)) 
    ndups <- y$printer$ndups
  if (is.null(ndups)) 
    ndups <- 1
  if (is.null(spacing)) 
    spacing <- y$printer$spacing
  if (is.null(spacing)) 
    spacing <- 1
  if (is.null(weights)) 
    weights <- y$weights
  method <- match.arg(method, c("ls", "robust"))
  if (ndups > 1) {
    if (!is.null(y$probes)) 
      y$probes <- uniquegenelist(y$probes, ndups = ndups, 
                                 spacing = spacing)
    if (!is.null(y$Amean)) 
      y$Amean <- rowMeans(unwrapdups(as.matrix(y$Amean), 
                                     ndups = ndups, spacing = spacing), na.rm = TRUE)
  }
  if (method == "robust") 
    fit <- mrlm2(y$exprs, design = design, ndups = ndups, 
                 spacing = spacing, weights = weights, ...)
  else if (ndups < 2 && is.null(block)) 
    fit <- lm.series(y$exprs, design = design, ndups = ndups, 
                     spacing = spacing, weights = weights)
  else {
    if (missing(correlation)) 
      stop("the correlation must be set, see duplicateCorrelation")
    fit <- gls.series(y$exprs, design = design, ndups = ndups, 
                      spacing = spacing, block = block, correlation = correlation, 
                      weights = weights, ...)
  }
  if (NCOL(fit$coefficients) > 1) {
    n <- rowSums(is.na(fit$coefficients))
    n <- sum(n > 0 & n < NCOL(fit$coefficients))
    if (n > 0) 
      warning("Partial NA coefficients for ", n, " probe(s)", 
              call. = FALSE)
  }
  fit$genes <- y$probes
  fit$Amean <- y$Amean
  fit$method <- method
  fit$design <- design
  new("MArrayLM", fit)
}
