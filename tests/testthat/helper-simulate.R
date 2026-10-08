# Simulation helpers shared by the test files.
#
# simulateIndividualModel(): data from the individual-level model
#   y_i = B' x_i + e_i,  with per-group residual scale (dispersion) and controlled angles between effect vectors.
# Returns list(meta, Y, D2) where D2 is the squared Euclidean sample-sample distance matrix.

unitVector <- function(p) { v <- rnorm(p); v / sqrt(sum(v^2)) }

# unit vector with a given cosine to `u`
mixedVector <- function(u, cosv, p) {
  w <- unitVector(p); w <- w - sum(w * u) * u; w <- w / sqrt(sum(w^2))
  cosv * u + sqrt(1 - cosv^2) * w
}

simulateIndividualModel <- function(n.per.group = c(A = 6, B = 6), p = 50,
                                    group.effect = 0.3, group.sd = NULL,
                                    batch.levels = NULL, batch.effect = 0, batch.cos = 0,
                                    age = FALSE, age.effect = 0, age.cos = 0,
                                    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  groups <- names(n.per.group)
  g <- factor(rep(groups, n.per.group), levels = groups)
  n <- length(g)
  meta <- data.frame(group = g, row.names = sprintf("s%02d", seq_len(n)))

  u.g <- unitVector(p)
  mu <- matrix(0, n, p)
  # group effects: reference level has none, others shift along u.g (scaled by index)
  for (k in seq_along(groups)[-1]) mu[g == groups[k], ] <- mu[g == groups[k], ] + rep((k - 1) * group.effect * u.g, each = sum(g == groups[k]))

  if (!is.null(batch.levels)) {
    meta$batch <- factor(sample(batch.levels, n, replace = TRUE), levels = batch.levels)
    u.b <- mixedVector(u.g, batch.cos, p)
    for (k in seq_along(batch.levels)[-1]) {
      idx <- meta$batch == batch.levels[k]
      mu[idx, ] <- mu[idx, ] + rep((k - 1) * batch.effect * u.b, each = sum(idx))
    }
  }
  if (age) {
    meta$age <- round(50 + rnorm(n, 0, 10))
    u.a <- mixedVector(u.g, age.cos, p)
    mu <- mu + outer((meta$age - mean(meta$age)) / 10, age.effect * u.a)
  }

  sdv <- rep(1, n)
  if (!is.null(group.sd)) sdv <- group.sd[as.character(g)]
  Y <- mu + matrix(rnorm(n * p), n, p) * sdv
  rownames(Y) <- rownames(meta)
  list(meta = meta, Y = Y, D2 = as.matrix(dist(Y))^2)
}

# raw pair-cell means of a squared distance matrix for two levels of a grouping factor
pairCellMeans <- function(D2, g, ref, alt) {
  m <- function(a, b) {
    D <- D2[g == a, g == b, drop = FALSE]
    if (a == b) mean(D[upper.tri(D)]) else mean(D)
  }
  c(RR = m(ref, ref), TT = m(alt, alt), RT = m(ref, alt))
}
