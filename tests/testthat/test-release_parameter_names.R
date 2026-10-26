.release_name_calls <- function() {
  dat <- data.frame(g = rep(1:2, each = 5), XB = seq_len(10))
  for (j in 1:4) {
    dat[[paste0("T", j)]] <- j - 1
    for (prefix in c("Y", "M", "X")) dat[[paste0(prefix, j)]] <- j + seq_len(10)
  }
  uni <- list(dat = dat, t_var = "T", y_var = "Y", curveFun = "LIN",
              intrinsic = FALSE, records = 1:4, paramOut = TRUE)
  lcsm <- uni
  lcsm$curveFun <- "NonP"
  list(
    getLGCM = uni,
    getLCSM = lcsm,
    getTVCmodel = c(uni, list(y_model = "LGCM", TVC = "M", decompose = 0)),
    getMGM = list(dat = dat, t_var = c("T", "T"), y_var = c("Y", "M"),
                  curveFun = "LIN", intrinsic = FALSE, records = rep(list(1:4), 2),
                  y_model = "LGCM", paramOut = TRUE),
    getMediation = list(dat = dat, t_var = c("T", "T"), y_var = "Y", m_var = "M",
                        x_var = "XB", x_type = "baseline", curveFun = "LIN",
                        records = rep(list(1:4), 2), paramOut = TRUE),
    getMGroup = c(uni, list(grp_var = "g", sub_Model = "LGCM")),
    getMIX = c(uni, list(prop_starts = 1, sub_Model = "LGCM"))
  )
}

test_that("all public fitting functions reject recyclable names before fitting", {
  # A regression must reject names rather than merely fail during estimation.
  local_mocked_bindings(mxRun = function(...) stop("optimizer was reached"),
                        mxTryHard = function(...) stop("optimizer was reached"),
                        .package = "nlpsem")
  for (fun in names(.release_name_calls())) {
    args <- .release_name_calls()[[fun]]
    args$starts <- list()
    for (nms in list("only_name", letters[1:2], letters[1:3])) {
      args$names <- nms
      expect_error(do.call(get(fun, asNamespace("nlpsem")), args),
                   "names must contain exactly", info = fun)
    }
    for (nms in list(character(), c("a", "a"), c("a", NA_character_), c("a", " "))) {
      args$names <- nms
      expect_error(do.call(get(fun, asNamespace("nlpsem")), args),
                   "non-missing, non-empty, unique", info = fun)
    }
    args$names <- NULL
    args$paramOut <- NA
    expect_error(do.call(get(fun, asNamespace("nlpsem")), args),
                 "paramOut must be a single logical", info = fun)
  }
})

test_that("standalone LCSM accepts only its documented complete or core names", {
  e <- .load_extdata("getLCSM_examples.RData")
  fit <- e$Read_LCSM_NonP
  full <- fit@Estimates$Name
  core <- head(full, length(full) - 6L * 8L)
  named <- nlpsem:::getLCSM.output(fit@mxOutput, "NonP", NULL, full)
  core_named <- nlpsem:::getLCSM.output(fit@mxOutput, "NonP", NULL, core)
  expect_identical(core_named, named)
  expect_error(nlpsem:::getLCSM.output(fit@mxOutput, "NonP", NULL, "one"),
               "names must contain exactly")
  expect_error(nlpsem:::getLCSM.output(fit@mxOutput, "NonP", NULL,
                                     c("slp_val_est1", core[-1])),
               "remain unique")
})

test_that("single and class-specific output cannot silently recycle parameter labels", {
  e <- .load_extdata("getLGCM_examples.RData")
  expect_error(nlpsem:::getLGCM.output(e$Math_LGCM_BLS_r@mxOutput, NULL, "one"),
               "names must contain exactly")
  mix <- .load_extdata("getMIX_examples.RData")$Math_BLS_LGCM3
  nms <- nlpsem:::.auto_names_LGCM("BLS", FALSE, NULL)
  extract <- function(nms) nlpsem:::getMIX.output(
    mix@mxOutput, 3L, "LGCM", NULL, "M", "BLS", NULL, 1:9, NULL, NULL, NULL, nms)
  expect_error(extract("one"), "names must contain exactly")
  expect_error(extract(rep("same", length(nms))), "unique")
  actual <- extract(nms)
  expect_equal(sum(startsWith(actual$Name, "c1")), length(nms))
  expect_identical(anyDuplicated(actual$Name), 0L)
})
