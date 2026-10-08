# Track B.1 / B.3: options, test grammar, reference heuristics, model object, metadata audit, legacy aliases.

makeMeta3 <- function(n = 24, seed = 3) {
  set.seed(seed)
  data.frame(condition = factor(rep(c("control", "disease"), each = n / 2)),
             batch = factor(rep(c("b1", "b2"), length.out = n)),
             sex = sample(c("F", "M"), n, replace = TRUE),
             age = round(rnorm(n, 50, 10)),
             stage = factor(rep(c("I", "II", "III"), length.out = n)),
             treated = rep(c(TRUE, FALSE), length.out = n),
             sample_id = sprintf("S%02d", seq_len(n)),
             library = sprintf("L%03d", seq_len(n)),
             pmi = c(rep(NA, 9), round(runif(n - 9, 2, 30))),
             protocol = "v3",
             row.names = sprintf("s%02d", seq_len(n)), stringsAsFactors = FALSE)
}

test_that("chooseReferenceLevel follows the D24 rules", {
  r <- chooseReferenceLevel(c("disease", "control", "disease", "control", "disease"))
  expect_equal(r$level, "control"); expect_equal(r$reason, "control-like")
  r <- chooseReferenceLevel(factor(c("IPF", "HC", "IPF", "HC")))
  expect_equal(r$level, "HC")
  r <- chooseReferenceLevel(factor(c("x", "y", "y", "y")))
  expect_equal(r$level, "y"); expect_equal(r$reason, "most-frequent")
  r <- chooseReferenceLevel(factor(c("x", "y", "y"), levels = c("y", "x")))   # explicit (non-alphabetical) order
  expect_equal(r$level, "y"); expect_equal(r$reason, "explicit-order")
  r <- chooseReferenceLevel(c("control", "ctrl", "x", "x"))                    # two control-like names: most frequent
  expect_equal(r$level, "x"); expect_equal(r$reason, "most-frequent")
  expect_equal(chooseReferenceLevel(c(TRUE, FALSE))$level, "FALSE")
})

test_that("every row of the test grammar resolves as specified", {
  meta <- makeMeta3()
  f <- ~ condition + batch + age
  one <- function(...) resolveTests(..., meta = meta, formula = f)
  t <- one("condition")[[1]]
  expect_equal(t$kind, "contrast"); expect_equal(t$contrast, c("condition", "disease", "control"))
  expect_equal(unname(t$levels), c("control", "disease")); expect_equal(t$reference.reason, "control-like")
  expect_match(t$label, "disease vs control")
  expect_match(t$interpretation[["shift"]], "disease samples differ from control samples")
  t <- one("stage")[[1]]
  expect_equal(t$kind, "term"); expect_equal(t$n.levels, 3)
  t <- one("age")[[1]]
  expect_equal(t$kind, "contrast"); expect_equal(t$step, 1)
  expect_equal(t$contrast$den, mean(meta$age)); expect_equal(t$contrast$num, mean(meta$age) + 1)
  t <- one("age", numeric.step = 10)[[1]]
  expect_equal(t$contrast$num - t$contrast$den, 10)
  t <- one("treated")[[1]]
  expect_equal(t$contrast, c("treated", "TRUE", "FALSE"))
  tt <- one(c("condition", "sex"))
  expect_length(tt, 2); expect_equal(vapply(tt, `[[`, character(1), "variable"), c("condition", "sex"))
  tt <- one("all")
  expect_equal(vapply(tt, `[[`, character(1), "variable"), c("condition", "batch", "age"))
  t <- one("condition: control vs disease")[[1]]
  expect_equal(t$contrast, c("condition", "control", "disease")); expect_equal(t$reference.reason, "user")
  t <- one(c("condition", "control", "disease"))[[1]]
  expect_equal(t$contrast, c("condition", "control", "disease"))
  t <- one(contrast = list(type = "simple", term = "condition", num = "disease", den = "control", at = list(batch = "b2")))[[1]]
  expect_equal(t$kind, "contrast"); expect_match(t$label, "at batch = b2")
  t <- one(contrast = c(age = 1))[[1]]
  expect_true(is.na(t$variable))
  expect_error(one("nope"), "not a column")
  expect_error(one(c("condition", "x", "control")), "not found")
  expect_error(one("condition", contrast = c("condition", "disease", "control")), "either")
})

test_that("describeMetadata assigns roles and notes", {
  meta <- makeMeta3()
  meta$skewed <- exp(rnorm(nrow(meta), 0, 2))
  d <- describeMetadata(meta)
  role <- setNames(d$role, d$column)
  expect_equal(role[["condition"]], "usable"); expect_equal(role[["age"]], "usable")
  expect_equal(role[["sample_id"]], "id-like"); expect_equal(role[["library"]], "id-like")
  expect_equal(role[["protocol"]], "constant")
  expect_equal(role[["pmi"]], "mostly-missing")
  expect_match(d$notes[d$column == "skewed"], "skew")
  expect_equal(d$n.missing[d$column == "pmi"], 9)
  meta$many <- factor(sprintf("x%d", rep(1:12, 2)))
  expect_equal(describeMetadata(meta)$role[describeMetadata(meta)$column == "many"], "high-cardinality")
  s <- cacoa:::formatMetadataSummary(d, meta)
  expect_match(s, "24 samples"); expect_match(s, "id-like")
})

test_that("buildCacoaModel builds designs per test, drops incomplete samples and prints", {
  meta <- makeMeta3()
  m <- buildCacoaModel(meta, formula = ~ condition + batch, test = "condition")
  expect_s3_class(m, "cacoaModel")
  expect_length(m$tests, 1)
  expect_equal(m$tests[[1]]$design$contrast_spec$type, "simple")
  expect_equal(colnames(m$F), colnames(m$tests[[1]]$design$F))        # primary design at top level
  expect_equal(deparse(m$dispersion.formula), "~condition")
  expect_equal(m$tests[[1]]$permutation$scheme, "block")
  expect_equal(m$tests[[1]]$permutation$n.strata, 2)
  expect_output(print(m), "disease vs control")
  expect_output(print(m), "adjusted for batch")
  expect_output(print(m), "reference 'control'")
  # default formula from the test variables
  m2 <- buildCacoaModel(meta, test = "condition", block.vars = "batch")
  expect_setequal(all.vars(m2$formula), c("condition", "batch")); expect_true(m2$formula.default)
  # numeric covariate -> Freedman-Lane in the plan
  m3 <- buildCacoaModel(meta, formula = ~ condition + age, test = "condition")
  expect_equal(m3$tests[[1]]$permutation$scheme, "freedman-lane")
  # missing covariate values -> listwise deletion with a warning issue
  m4 <- buildCacoaModel(meta, formula = ~ condition + pmi, test = "condition")
  expect_equal(length(m4$samples$dropped), 9)
  expect_true(any(grepl("dropped", m4$issues$message)))
  # multiple tests and a term test
  m5 <- buildCacoaModel(meta, formula = ~ condition + stage + batch, test = c("condition", "stage"))
  expect_length(m5$tests, 2)
  expect_equal(m5$tests[[2]]$kind, "term")
  expect_equal(m5$tests[[2]]$design$term.cols, c("stageII", "stageIII"))
  # ID-like column in the formula is an error
  expect_error(buildCacoaModel(meta, formula = ~ condition + sample_id, test = "condition"), "unique value per sample")
  # aliasing: condition fully determined by a nested covariate
  meta$site <- ifelse(meta$condition == "control", "A", "B")
  m6 <- buildCacoaModel(meta, formula = ~ condition + site, test = "condition")
  expect_true(any(m6$issues$severity == "error"))
  expect_match(m6$issues$message[m6$issues$severity == "error"], "fully determined by 'site'")
})

test_that("Cacoa options: defaults, setOptions validation, per-call override precedence", {
  cao <- makeToyCacoa()
  expect_equal(cao$getOptions("n.permutations"), 999)
  expect_equal(cao$getOptions("seed"), 1)
  cao$setOptions(n.permutations = 49, seed = NULL, n.cores = 2)
  expect_equal(cao$getOptions("n.permutations"), 49)
  expect_null(cao$getOptions("seed"))
  expect_equal(cao$n.cores, 2)
  expect_error(cao$setOptions(alpha = 2), "alpha")
  expect_error(cao$setOptions(nonsense = 1), "unknown option")
  expect_error(cao$setOptions(dist = "cosine"), "dist")
})

test_that("Cacoa constructor: optional model, setModel after exploration, legacy aliases and sample groups", {
  cao <- makeToyCacoa()
  expect_s3_class(cao$model, "cacoaModel")
  expect_equal(cao$ref.level, "A"); expect_equal(cao$target.level, "B")
  expect_equal(names(cao$sample.groups.palette), c("B", "A"))
  sg <- cao$getSampleGroups()
  expect_s3_class(sg, "factor"); expect_equal(levels(sg), c("A", "B")); expect_length(sg, 8)
  expect_equal(unname(cao$sample.groups["B1"]), factor("B", levels = c("A", "B")))
  # no model at construction: metadata audit is printed; methods needing a model say so
  expect_message(cao0 <- makeToyCacoa(contrast = NULL, verbose = TRUE, suppress.messages = FALSE), "Sample metadata")
  expect_null(cao0$model)
  expect_error(cao0$estimateExpressionShiftMagnitudes(n.permutations = 9, verbose = FALSE), "no model is set")
  cao0$setModel(test = "group", verbose = FALSE)
  expect_equal(cao0$model$tests[[1]]$contrast, c("group", "B", "A"))
  expect_equal(cao0$target.level, "B")
  # changing the model after exploration
  cao0$setModel(~ group + batch, test = "group: A vs B", verbose = FALSE)
  expect_equal(cao0$ref.level, "B"); expect_equal(cao0$target.level, "A")
  expect_true("batch" %in% all.vars(cao0$formula))
  # the old two-group constructor (vignette) still runs
  sg.old <- setNames(rep(c("ctrl", "dis"), each = 4), rownames(cao$sample.meta))
  expect_warning(cao.old <- makeToyCacoa(contrast = NULL, sample.groups = sg.old, ref.level = "ctrl", target.level = "dis", suppress.warnings = FALSE), "deprecated")
  expect_equal(cao.old$model$tests[[1]]$contrast, c("condition", "dis", "ctrl"))
  expect_equal(cao.old$ref.level, "ctrl")
  expect_true("condition" %in% names(cao.old$sample.meta))
})
