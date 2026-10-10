# Default formula (D25), transformed covariates in endpoint rows (B8), endpoint settings, stored model fields.

makeMeta2 <- function(n = 24, seed = 11) {
  set.seed(seed)
  data.frame(group = factor(rep(c("A", "B"), each = n / 2)),
             batch = factor(rep(c("b1", "b2"), length.out = n)),
             age = round(rnorm(n, 50, 10)),
             sample_id = sprintf("S%02d", seq_len(n)),
             row.names = sprintf("s%02d", seq_len(n)))
}

test_that("default formula uses the contrast variables (and block.vars), never all metadata columns", {
  meta <- makeMeta2()
  expect_message(d <- cacoa:::buildDesignMatrices(meta, contrast = c("group", "B", "A")), "No formula supplied")
  expect_equal(colnames(d$F), c("groupA", "groupB"))
  d2 <- suppressMessages(cacoa:::buildDesignMatrices(meta, contrast = c("group", "B", "A"), blockVars = "batch"))
  expect_true(all(c("groupA", "groupB", "batchb2") %in% colnames(d2$F)))
  expect_false(any(grepl("sample_id", colnames(d2$F))))
  expect_error(cacoa:::buildDesignMatrices(meta, contrast = c(age = 1)), "coefficient-level contrast")
  expect_error(cacoa:::buildDesignMatrices(meta, contrast = c("sex", "M", "F")), "not found in sample metadata")
})

test_that("transformed covariates work in the design and in endpoint rows (B8)", {
  meta <- makeMeta2()
  for (f in list(~ group + log(age), ~ group + poly(age, 2), ~ group + I(age^2) + batch, ~ group + factor(batch))) {
    d <- cacoa:::buildDesignMatrices(meta, contrast = c("group", "B", "A"), formula = f)
    expect_true(all(is.finite(d$contrast.F)), info = deparse(f))
    expect_true(all(is.finite(d$contrast_endpoints_F$num)), info = deparse(f))
    expect_equal(unname(d$contrast_endpoints_F$num - d$contrast_endpoints_F$den), unname(d$contrast.F), info = deparse(f))
    # the group columns carry the contrast, nothing else does
    nz <- names(d$contrast.F)[abs(d$contrast.F) > 1e-12]
    expect_true(all(grepl("^group", nz)), info = deparse(f))
  }
  # poly() columns must come from the same basis as model.matrix on the full data
  d <- cacoa:::buildDesignMatrices(meta, contrast = c("group", "B", "A"), formula = ~ group + poly(age, 2))
  mm <- model.matrix(~ group + poly(age, 2), meta)
  expect_equal(unname(d$F[, grep("poly", colnames(d$F))]), unname(mm[, grep("poly", colnames(mm))]))
  # numerics not in the contrast sit at their mean in the endpoint rows
  d <- cacoa:::buildDesignMatrices(meta, contrast = c("group", "B", "A"), formula = ~ group + log(age))
  expect_equal(unname(d$contrast_endpoints_F$num["log(age)"]), log(mean(meta$age)))
})

test_that("endpoint settings reproduce the endpoint rows on the location design", {
  meta <- makeMeta2()
  meta$group <- factor(rep(c("A", "B", "C"), length.out = nrow(meta)))
  specs <- list(list(type = "simple", term = "group", num = "B", den = "A", at = list(batch = "b2")),
                list(type = "marginal", term = "group", num = "B", den = "A", over = "batch", weights = "proportional"),
                list(type = "simple", term = "group:batch", num = "B:b1", den = "A:b2"))
  for (sp in specs) {
    d <- cacoa:::buildDesignMatrices(meta, contrast = sp, formula = ~ group * batch + age)
    expect_false(is.null(d$contrast_endpoints_at))
    rows <- cacoa:::endpointRowsFromDesign(d$F, meta, d$contrast_endpoints_at, numericRef = d$numeric_ref_used)
    expect_equal(rows$num, d$contrast_endpoints_F$num)
    expect_equal(rows$den, d$contrast_endpoints_F$den)
    # weights sum to one on each side
    expect_equal(sum(sapply(d$contrast_endpoints_at$num, `[[`, "w")), 1)
  }
  d <- cacoa:::buildDesignMatrices(meta, contrast = c(age = 1), formula = ~ group + age)
  expect_null(d$contrast_endpoints_at)
  # evaluated on a different design (here: group only), the settings give that design's endpoint rows
  d <- cacoa:::buildDesignMatrices(meta, contrast = c("group", "B", "A"), formula = ~ group + batch + age)
  Fd <- cacoa:::buildFullDesign(~ group, meta)
  rows <- cacoa:::endpointRowsFromDesign(Fd, meta, d$contrast_endpoints_at)
  expect_equal(unname(rows$num - rows$den), c(0, 1, 0))
})

test_that("Cacoa object stores the model actually used and numeric.ref (B11, B20)", {
  cao <- makeToyCacoa()
  expect_s3_class(cao, "Cacoa")
  expect_true(inherits(cao$model$formula, "formula"))
  expect_equal(all.vars(cao$model$formula), "group")
  expect_equal(cao$numeric.ref, "auto")
  expect_equal(cao$model$tests[[1]]$contrast, c("group", "B", "A"))
  expect_null(cao$formula); expect_null(cao$contrast)                        # legacy fields are gone; the model holds them
  cao2 <- makeToyCacoa(formula = ~ group + batch, numeric.ref = list(age = 40))
  expect_setequal(all.vars(cao2$model$formula), c("group", "batch"))
  expect_equal(cao2$numeric.ref, list(age = 40))
})
