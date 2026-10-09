# Option-grid helper for the kernel tests (engine_convergence.md §11): one metadata table per grid cell and a
# feature matrix with optional missing-data patterns.
gridMeta <- function(levels = 2, n.per.level = 6, batch = TRUE, age = FALSE, seed = 1) {
  set.seed(seed)
  n <- levels * n.per.level
  meta <- data.frame(group = factor(rep(paste0("G", seq_len(levels)), each = n.per.level)), row.names = sprintf("s%02d", seq_len(n)))
  if (batch) meta$batch <- factor(rep(c("b1", "b2"), length.out = n))
  if (age) meta$age <- round(rnorm(n, 50, 8))
  meta
}
gridY <- function(meta, m = 5, na = c("none", "pattern"), effect = 0, seed = 2) {
  na <- match.arg(na); set.seed(seed)
  n <- nrow(meta); Y <- matrix(rnorm(n * m), n, m, dimnames = list(rownames(meta), paste0("y", seq_len(m))))
  Y[meta$group == levels(meta$group)[2], ] <- Y[meta$group == levels(meta$group)[2], ] + effect
  if (na == "pattern") { Y[c(1, 3), 2] <- NA; Y[c(2, 5, n), 3] <- NA; if (m >= 4) Y[c(1, 3), 4] <- NA }   # columns 2 and 4 share a pattern
  Y
}
gridDesign <- function(meta, kind = c("two.level", "term"), formula = NULL) {
  kind <- match.arg(kind)
  levs <- levels(meta$group)
  if (is.null(formula)) formula <- if ("batch" %in% names(meta)) ~ group + batch else ~ group
  if (kind == "two.level") suppressMessages(buildDesignMatrices(meta, contrast = c("group", levs[2], levs[1]), formula = formula))
  else buildCacoaModel(meta, formula = formula, test = "group")
}
# brute-force reference: contrast statistic of the OLS fit of y on the permuted design X[q, ]
refContrastStat <- function(X, y, cvec, q = seq_len(nrow(X)), w = NULL) {
  Xq <- X[q, , drop = FALSE]
  if (is.null(w)) b <- solve(crossprod(Xq), crossprod(Xq, y)) else b <- solve(crossprod(Xq, w * Xq), crossprod(Xq, w * y))
  sum(cvec * b)
}
