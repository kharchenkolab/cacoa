getColumnwiseCorrelations <- function(m1, m2, corrected=F) {
  m1 <- m1[, intersect(colnames(m1), colnames(m2))]
  m2 <- m2[, intersect(colnames(m1), colnames(m2))]

  if (corrected) {
    m1 <- m1 - getMinCorrectionVector(m1, m2)
    if (!all(abs(getMinCorrectionVector(m1, m1)) < 1e-10)) stop("Wrong correction")
  }
  return(lineup::corbetw2mat(m1, m2))
}

sn <- function(x) {setNames(x, x)}
