#' @importFrom sccore checkPackageInstalled
NULL

# Joint sample distance over cell types: the average of the per-type distance matrices. `weighting = "cells"` weights a
# pair's distance in a type by sqrt(min(n.cells)) of the two samples (types with many cells in both samples count more);
# `"equal"` gives every type in which both samples are present the same weight. Pairs that share no type are NA.
#' @keywords internal
prepareJointExpressionDistance <- function(p.dist.per.type, sample.groups=NULL, return.dists=TRUE, weighting=c("cells", "equal")) {
  weighting <- match.arg(weighting)
  pairWeight <- function(nc) if (weighting == "cells") sqrt(outer(nc, nc, FUN = 'pmin')) else 1 * outer(nc > 0, nc > 0, "&")
  checkPackageInstalled(c("abind"), cran=TRUE)
  # bring to a common set of cell types
  common.types <- lapply(p.dist.per.type, colnames) %>% unlist() %>% unique()

  p.dist.per.type %<>% lapply(function(x) {
    y <- matrix(0,nrow=length(common.types),ncol=length(common.types));  # can set the missing entries to zero, as they will carry zero weights
    rownames(y) <- colnames(y) <- common.types;
    y[rownames(x),colnames(x)] <- x;
    ycct <- setNames(rep(0,length(common.types)), common.types);
    ycct[colnames(x)] <- attr(x, 'n.cells')
    attr(y, 'n.cells') <- ycct
    y
  }) # reform the matrix to make sure all cell type have the same dimensions

  x <- abind::abind(lapply(p.dist.per.type, function(x) {
    return(x * pairWeight(attr(x, 'n.cells')))
  }), along = 3)

  # just the weights (for total sum of weights normalization)
  y <- abind::abind(lapply(p.dist.per.type, function(x) {
    pairWeight(attr(x, 'n.cells'))
  }), along = 3)

  # normalize by total weight sums
  xd <- apply(x, c(1, 2), sum) / apply(y, c(1, 2), sum)

  if (return.dists)
    return(xd)

  cross.factor <- outer(sample.groups[rownames(xd)], sample.groups[colnames(xd)], '==')
  diag(xd) <- NA # remove self pairs
  # restrict
  xd[!cross.factor] <- NA
  if (!any(!is.na(xd)))
    return(NULL)
  xmd2 <- na.omit(reshape2::melt(xd))
  xmd2 <- na.omit(xmd2)
  xmd2$type1 <- sample.groups[as.character(xmd2$Var1)]
  xmd2$type2 <- sample.groups[as.character(xmd2$Var2)]

  return(xmd2)
}


# x0 from factor-contrast weights applied to the predictor x (mirrors total/shift)
#' @keywords internal
contrastCenterX0 <- function(x, pair.level, w, tol = 1e-8) {
  pos.lv <- names(w)[w >  tol]
  neg.lv <- names(w)[w < -tol]
  if (!length(pos.lv)) stop("Factor centering contrast has no positive weights.")
  denom <- sum(w[pos.lv]); if (denom <= 0) stop("Sum of positive weights must be > 0.")

  mu.x.neg <- lapply(neg.lv, function(L) {
    idx <- (pair.level == L)
    if (!any(idx)) return(NA_real_)
    mean(x[idx], na.rm = TRUE)
  })
  names(mu.x.neg) <- neg.lv

  num <- sum(mapply(function(L, mu) { (-w[L]) * mu }, neg.lv, mu.x.neg))
  num / denom
}

#' @keywords internal
parseDistance <- function(dist, n.dims = NULL) {
  n.comps <- min(n.dims, Inf)
  if (is.null(dist)) {
    dist <- ifelse(n.comps < 20, 'l1', 'cor')
    return(dist)
  }

  dist %<>% tolower()
  if (dist == 'l2') {
    warning("Using dist='l2' is not recommended, as it may introduce unwanted dependency ",
            "on the number of cells per cluster. Please, consider using 'l1' instead.")
  } else if (dist == 'cor') {
    if (n.comps < 20) {
      warning("dist='cor' is not recommended for data with dimensionality < 20. ",
              "Please, consider using 'l1' instead.")
    }
  } else if (dist == 'l1') {
    if (n.comps > 30) {
      warning("dist='l1' is not recommended for data with dimensionality > 30. ",
              "Please, consider using 'cor' instead.")
    }
  } else {
    stop("Unknown dist: ", dist)
  }

  return(dist)
}
# sample-space names: "expression" (per cell type or all cell types combined), "pseudobulk" (all cells pooled per sample),
# "composition" (ILR of the cell-type counts); the earlier names are accepted
#' @keywords internal
normalizeSampleSpace <- function(space) {
  space <- space[1]
  alias <- c(expression.shifts = "expression", coda = "composition", pseudo.bulk = "pseudobulk")
  if (space %in% names(alias)) space <- alias[[space]]
  if (!space %in% c("expression", "pseudobulk", "composition"))
    stop("space must be one of 'expression', 'pseudobulk', 'composition' (got '", space, "')")
  space
}

# row labels of the combined sample spaces in the covariate screen
SAMPLE.SPACE.LABELS <- list(expression = "all cell types", expression.equal = "all cell types (equal weights)",
                            pseudobulk = "pooled pseudobulk", composition = "composition")
