
library(cacoa)
library(testthat)
library(magrittr)
library(conos)

con <- Conos$new(panel.preprocessed, n.cores = 1)
con$buildGraph(n.odgenes = 500)
con$findCommunities()
con$embedGraph()

sample.groups <- c("CTRL","CTRL","DISEASE","DISEASE") %>% 
  `names<-`(names(con$samples))
cell.groups <- con$clusters$leiden$groups

cao <- Cacoa$new(data.object = con, sample.groups = sample.groups, cell.groups = cell.groups, ref.level = "CTRL", target.level = "DISEASE", n.cores = 1)

test_that("dummy test", {
	expect_equal(1, 1)
})
