#' Estimate Gene Programmes with FABIA
#' @param n.programmes maximal number of gene programmes to find (parameter `p` for fabia).
#' @param n.sampled.cells number of sub-sampled cells for estimating the gene programmes. If 0, all cells are used.
#' it is interpreted as a vector of cell
#' @inheritDotParams fabia::fabia -p -X -cyc -alpha -random
estimateGeneProgrammesFabia <- function(z.scores, n.programmes, n.sampled.cells=15000, min.z=0.5, max.z=3,
                                        cyc=1500, alpha=0.2, random=-1, ...) {
  if (!requireNamespace("fabia", quietly=TRUE))
    stop("fabia package must be installed to run this function")

  z.scores@x[is.na(z.scores@x)] <- 0
  if (min.z > 1e-10) {
    z.scores.bin <- z.scores
    z.scores.bin@x %<>% {as.numeric(abs(.) >= min.z)}
    genes.filt <- which(colMeans(z.scores.bin) > 0.01)
    z.scores <- z.scores[, genes.filt]
  }

  if (n.sampled.cells > 0) {
    sample.ids <- unique(round(seq(1, nrow(z.scores), length.out=n.sampled.cells)))
    z.scores <- z.scores[sample.ids,]
  }

  z.scores@x %<>% pmax(-max.z) %<>% pmin(max.z)

  fabia.res <- z.scores %>% t() %>%
    fabia::fabia(p=n.programmes, alpha=alpha, cyc=cyc, random=random, ...)

  return(list(fabia=fabia.res, sample.ids=sample.ids))
}
