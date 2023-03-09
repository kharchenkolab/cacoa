## code to prepare `panel.preprocessed` dataset goes here

library(pagoda2)
library(magrittr)

sample.names <- c("CTRL1","CTRL2","DISEASE1","DISEASE2")

panel.preprocessed <- lapply(seq_len(4), \(x) {
  out <- rsparsematrix(5e2, 50, 0.1)
  out[out < 0] <- 1
  dimnames(out) <- list(sapply(1:5e2, \(x) paste0("gene",x)), sapply(1:50, \(x) paste0("cell",x)))

  return(out)
}) %>%
  {lapply(seq_len(4), \(s) .[[s]] %>% `colnames<-`(., paste0(sample.names[s],"!!",colnames(.))))} %>%
  lapply(basicP2proc, min.cells.per.gene = 0, min.transcripts.per.cell = 0, get.largevis = F, get.tsne = F, make.geneknn = F, n.cores = 1) %>%
  `names<-`(sample.names)

usethis::use_data(panel.preprocessed, overwrite = TRUE, compress = "xz")
