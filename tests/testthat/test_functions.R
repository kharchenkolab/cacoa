
library(cacoa)
library(testthat)
library(magrittr)
library(conos)

# Example panel was made this way:
#
# library(pagoda2)
# library(magrittr)
# panel.processed <- lapply(seq_len(4), \(x) {
#   out <- rsparsematrix(5e2, 50, 0.1)
#   out[out < 0] <- 1
#   dimnames(out) <- list(sapply(1:5e2, \(x) paste0("gene",x)), sapply(1:50, \(x) paste0("cell",x)))
#   
#   return(out)
# })
# 
# panel.processed %<>% 
#   names() %>% 
#   lapply(\(sname) testdata.cms[[sname]] %>% `colnames<-`(., paste0(sname,"!!",colnames(.))))
# 
# panel.preprocessed %<>% `names<-`(c("CTRL1","CTRL2"))
# 
# panel.preprocessed <- testdata.cms %>% 
#   lapply(basicP2proc, min.cells.per.gene = 0, min.transcripts.per.cell = 0, get.largevis = F, get.tsne = F, make.geneknn = F, n.cores = 1) %>% 
#   `names<-`(c("sample1","sample2","sample3","sample4"))

con <- Conos$new(panel.processed, n.cores = 1)
con$buildGraph()
con$findCommunities()
con$embedGraph()

sample.groups <- c("CTRL","CTRL","DISEASE","DISEASE") %>% `names<-`(names(con$samples))
cell.groups <- con$clusters$leiden$groups

cao <- Cacoa$new(data.object = con, sample.groups = sample.groups, cell.groups = cell.groups, ref.level = "CTRL", target.level = "DISEASE", n.cores = 1)

test_that("dummy test", {
	expect_equal(1, 1)
})
