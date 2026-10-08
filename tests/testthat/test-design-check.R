# Track B.4: design check from metadata alone.

makeMeta4 <- function(n = 30, seed = 5) {
  set.seed(seed)
  cond <- factor(rep(c("control", "disease"), each = n / 2))
  data.frame(condition = cond,
             batch = factor(rep(c("b1", "b2", "b3"), length.out = n)),
             site = factor(ifelse(cond == "control", "A", "B")),                      # fully determines condition
             medication = factor(ifelse(cond == "disease" & runif(n) < 0.85, "yes", ifelse(runif(n) < 0.1, "yes", "no"))),  # strongly tied
             sex = factor(sample(c("F", "M"), n, replace = TRUE)),
             age = round(rnorm(n, 50, 10)),
             age2 = NA,
             sample_id = sprintf("S%02d", seq_len(n)),
             row.names = sprintf("s%02d", seq_len(n)))
}

test_that("association measures behave", {
  m <- makeMeta4(); m$age2 <- m$age + rnorm(nrow(m), 0, 1)
  expect_equal(covariateAssociation(m$condition, m$site), 1, tolerance = 1e-8)
  expect_gt(covariateAssociation(m$condition, m$medication), 0.5)
  expect_lt(covariateAssociation(m$condition, m$batch), 0.3)
  expect_gt(covariateAssociation(m$age, m$age2), 0.9)
  expect_true(covariateAssociation(m$condition, m$age) >= 0 && covariateAssociation(m$condition, m$age) <= 1)
  expect_true(is.na(covariateAssociation(factor(rep("a", 5)), 1:5)))
  # bias correction: independent variables give small V
  set.seed(1); x <- factor(sample(letters[1:3], 300, TRUE)); y <- factor(sample(LETTERS[1:4], 300, TRUE))
  expect_lt(cacoa:::cramersV(x, y), 0.15)
})

test_that("checkDesign finds aliasing, over-adjustment, df and permutation issues", {
  m <- makeMeta4(); m$age2 <- NULL
  chk <- checkDesign(m, test.variable = "condition")
  expect_s3_class(chk, "cacoaDesignCheck")
  expect_equal(dim(chk$associations), c(6, 6))
  expect_true(any(chk$issues$severity == "error" & grepl("fully determined by 'site'", chk$issues$message)))
  expect_true(any(grepl("'medication' is strongly associated", chk$issues$message)))
  expect_false(any(grepl("sample_id", chk$covariates)))
  expect_true(is.data.frame(chk$balance$batch))
  expect_output(print(chk), "Error")
  expect_output(print(chk), "association with 'condition'")
  # with a model: df budget, GVIF, permutation plan
  mod <- buildCacoaModel(m, formula = ~ condition + batch + sex + age, test = "condition")
  chk2 <- checkDesign(m, model = mod)
  expect_equal(chk2$df$parameters, 6)
  expect_equal(nrow(chk2$permutation), 1)
  expect_equal(chk2$permutation$scheme, "freedman-lane")
  expect_true(all(is.finite(chk2$gvif)))
  expect_true(any(grepl("Freedman-Lane", chk2$issues$message)))
  # few permutations are flagged
  small <- m[c(1:4, 16:19), ]; small$batch <- factor(rep(c("b1", "b2"), 4))
  mod3 <- buildCacoaModel(small, formula = ~ condition + batch, test = "condition")
  chk3 <- checkDesign(small, model = mod3)
  expect_true(any(grepl("distinct permutations", chk3$issues$message)))
  # plots render
  expect_s3_class(plotDesignCheck(chk, "associations"), "ggplot")
  expect_s3_class(plotDesignCheck(chk, "issues"), "ggplot")
  expect_true(inherits(plotDesignCheck(chk, "balance", meta = m, covariates = c("batch", "age")), "ggplot"))
})

test_that("Cacoa$checkDesign and plotDesign work with and without a stored model", {
  cao <- makeToyCacoa(contrast = NULL)
  expect_null(cao$model)
  chk <- cao$checkDesign(test = "group", verbose = FALSE)
  expect_equal(chk$test.variable, "group")
  expect_false(is.null(cao$cache$design.check))
  cao$setModel(~ group + batch, test = "group", verbose = FALSE)
  expect_output(chk2 <- cao$checkDesign(verbose = TRUE), "permutations for 'group: B vs A'")
  expect_equal(chk2$permutation$scheme, "block")
  expect_s3_class(cao$plotDesign("associations"), "ggplot")
  expect_true(inherits(cao$plotDesign("balance"), "ggplot"))
})
