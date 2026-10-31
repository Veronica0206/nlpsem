# Regression tests driven by the fitted example objects shipped in inst/extdata. They exercise
# the plotting, output, factor-score, likelihood-ratio and summary paths without re-fitting
# any model, so they run in seconds.

xs <- 60

test_that("getFigure draws single-process LGCM, LCSM and TVC fits without console noise", {
  eL <- .load_extdata("getLGCM_examples.RData")
  eC <- .load_extdata("getLCSM_examples.RData")
  eT <- .load_extdata("getTVCmodel_examples.RData")

  expect_silent(f1 <- getFigure(
    eL$Math_LGCM_TIC_BLS_f@mxOutput, sub_Model = "LGCM", y_var = "M", curveFun = "BLS",
    y_model = "LGCM", t_var = "T", records = 1:9, xstarts = xs, xlab = "Month", outcome = "Math"
  ))
  expect_s4_class(f1, "figOutput")
  expect_length(f1@figures, 1)
  expect_length(f1@figures[[1]], 1)
  expect_true(inherits(f1@figures[[1]][[1]], "ggplot"))

  f2 <- getFigure(
    eC$Read_LCSM_NonP@mxOutput, sub_Model = "LCSM", y_var = "R", curveFun = "NonP", y_model = "LCSM",
    t_var = "T", records = 1:9, xstarts = xs / 12, xlab = "Year", outcome = "Reading"
  )
  expect_length(f2@figures[[1]], 2)
  f3 <- getFigure(
    eC$Read_LCSM_JB_r@mxOutput, sub_Model = "LCSM", y_var = "R", curveFun = "JB", y_model = "LCSM",
    t_var = "T", records = 1:9, xstarts = xs / 12, xlab = "Year", outcome = "Reading"
  )
  expect_length(f3@figures[[1]], 2)
  f4 <- getFigure(
    eT$Math_TVCslp_BLS_f@mxOutput, sub_Model = "TVC", y_var = "M", curveFun = "BLS", y_model = "LGCM",
    t_var = "T", records = 1:9, xstarts = xs, xlab = "Year", outcome = "Math"
  )
  expect_true(inherits(f4@figures[[1]][[1]], "ggplot"))
})

test_that("getFigure draws multivariate, mediation, multiple-group and mixture fits", {
  eG <- .load_extdata("getMGM_examples.RData")
  eD <- .load_extdata("getMediation_examples.RData")
  eR <- .load_extdata("getMGroup_examples.RData")
  eM <- .load_extdata("getMIX_examples.RData")

  f5 <- getFigure(
    eG$RM_PLGCM.f@mxOutput, sub_Model = "MGM", y_var = c("R", "M"), curveFun = "BLS", y_model = "LGCM",
    t_var = c("T", "T"), records = list(1:9, 1:9), xstarts = xs, xlab = "Month", outcome = c("Reading", "Math")
  )
  expect_length(f5@figures, 2)

  expect_silent(f6 <- getFigure(
    eD$Med2_LGCM_BLS@mxOutput, sub_Model = "MED", y_var = "M", m_var = "R", x_type = "baseline",
    x_var = "ex1", curveFun = "BLS", t_var = rep("T", 2), records = list(1:9, 1:9),
    xstarts = xs, xlab = "Month", outcome = c("Math", "Reading")
  ))
  expect_length(f6@figures, 2)
  expect_silent(f7 <- getFigure(
    eD$Med3_LGCM_BLS@mxOutput, sub_Model = "MED", y_var = "S", m_var = "M", x_type = "longitudinal",
    x_var = "R", curveFun = "BLS", t_var = rep("T", 3), records = list(2:9, 1:9, 1:9),
    xstarts = xs, xlab = "Month", outcome = c("Science", "Math", "Reading")
  ))
  expect_length(f7@figures, 3)
  for (fig in c(f6@figures, f7@figures)) expect_true(inherits(fig[[1]], "ggplot"))

  f8 <- getFigure(
    eR$MGroup_Math_BLS_LGCM_f@mxOutput, nClass = 2, grp_var = "SEX", sub_Model = "LGCM", y_var = "M",
    curveFun = "BLS", y_model = "LGCM", t_var = "T", records = 1:9, xstarts = xs, xlab = "Month", outcome = "Math"
  )
  expect_true(inherits(f8@figures[[1]][[1]], "ggplot"))
  f9 <- getFigure(
    eM$Math_BLS_LGCM3@mxOutput, nClass = 3, sub_Model = "LGCM", y_var = "M", curveFun = "BLS", y_model = "LGCM",
    t_var = "T", records = 1:9, xstarts = xs, xlab = "Month", outcome = "Math"
  )
  expect_true(inherits(f9@figures[[1]][[1]], "ggplot"))
  f10 <- getFigure(
    eM$RM_BLS_PLGCM3@mxOutput, nClass = 3, cluster_TIC = c("gx1", "gx2"), sub_Model = "MGM",
    y_var = c("R", "M"), curveFun = "BLS", y_model = "LGCM", t_var = c("T", "T"),
    records = list(1:9, 1:9), xstarts = xs, xlab = "Month", outcome = c("Reading", "Math")
  )
  expect_length(f10@figures, 2)
})

test_that("mediation trajectories are built from the role-named growth factors and the knot", {
  eD <- .load_extdata("getMediation_examples.RData")
  m <- eD$Med2_LGCM_BLS@mxOutput
  fig <- getFigure(
    m, sub_Model = "MED", y_var = "M", m_var = "R", x_type = "baseline", x_var = "ex1", curveFun = "BLS",
    t_var = rep("T", 2), records = list(1:9, 1:9), xstarts = xs, xlab = "Month", outcome = c("Math", "Reading")
  )
  for (role in c("Y", "M")) {
    gf <- c(OpenMx::mxEvalByName(paste0(role, "_mean"), m))   # slope 1, value at knot, slope 2
    knot <- c(OpenMx::mxEvalByName(paste0(role, "_mug"), m))
    est <- fig@figures[[match(role, c("Y", "M"))]][[1]]$data
    t0 <- est$time[1] - xs
    expect_equal(est$Y_hat[1], gf[2] + gf[1] * (t0 - knot))
    at_knot <- approx(est$time, est$Y_hat, xout = knot + xs)$y
    expect_lt(abs(at_knot - gf[2]), 0.05 * (abs(gf[1]) + abs(gf[3])) + 1e-6)
  }
})

test_that("output helpers reproduce the shipped estimate tables at full precision", {
  eL <- .load_extdata("getLGCM_examples.RData")
  eC <- .load_extdata("getLCSM_examples.RData")
  eT <- .load_extdata("getTVCmodel_examples.RData")
  eG <- .load_extdata("getMGM_examples.RData")
  eD <- .load_extdata("getMediation_examples.RData")

  obj <- eL$Math_LGCM_TIC_BLS_f
  expect_silent(out <- nlpsem:::getLGCM.output(obj@mxOutput, growth_TIC = c("ex1", "ex2"), names = obj@Estimates$Name))
  expect_equal((out), obj@Estimates)

  lcsm <- list(Read_LCSM_NonP = list("nonparametric", NULL),
               Read_LCSM_NonP_TIC = list("nonparametric", c("ex1", "ex2")),
               Read_LCSM_QUAD = list("quadratic", NULL),
               Read_LCSM_EXP_r = list("negative exponential", NULL),
               Read_LCSM_JB_r = list("Jenss-Bayley", NULL))
  for (nm in names(lcsm)) {
    obj <- eC[[nm]]
    expect_equal((nlpsem:::getLCSM.output(obj@mxOutput, curveFun = lcsm[[nm]][[1]], growth_TIC = lcsm[[nm]][[2]],
                                         names = obj@Estimates$Name)), obj@Estimates, info = nm)
  }

  tvc <- c(Math_TVCslp_BLS_f = 1, Math_TVCslp_BLS_r = 1, Math_TVCchg_BLS_r = 2, Math_TVCchgBL_BLS_r = 3)
  for (nm in names(tvc)) {
    obj <- eT[[nm]]
    expect_equal((nlpsem:::getTVC.output(obj@mxOutput, curveFun = "BLS", records = 1:9, y_model = "LGCM",
                                        decompose = tvc[[nm]], growth_TIC = "ex1", names = obj@Estimates$Name)),
                 obj@Estimates, info = nm)
  }

  for (nm in c("RM_PLGCM.r", "RM_PLGCM.f")) {
    obj <- eG[[nm]]
    expect_equal((nlpsem:::getMGM.output(obj@mxOutput, y_var = c("R", "M"), records = list(1:9, 1:9), curveFun = "BLS",
                                        y_model = "LGCM", names = obj@Estimates$Name)), obj@Estimates, info = nm)
  }

  obj <- eD$Med2_LGCM_BLS
  expect_equal((nlpsem:::getMED.output(obj@mxOutput, y_var = "M", m_var = "R", x_type = "baseline", x_var = "ex1",
                                      curveFun = "bilinear spline", names = obj@Estimates$Name)), obj@Estimates)
  obj <- eD$Med3_LGCM_BLS
  expect_equal((nlpsem:::getMED.output(obj@mxOutput, y_var = "S", m_var = "M", x_type = "longitudinal", x_var = "R",
                                      curveFun = "bilinear spline", names = obj@Estimates$Name)), obj@Estimates)
})

test_that("automatic parameter names match the shipped estimate tables", {
  eL <- .load_extdata("getLGCM_examples.RData")
  eC <- .load_extdata("getLCSM_examples.RData")
  eT <- .load_extdata("getTVCmodel_examples.RData")
  eG <- .load_extdata("getMGM_examples.RData")
  eD <- .load_extdata("getMediation_examples.RData")
  eR <- .load_extdata("getMGroup_examples.RData")
  eM <- .load_extdata("getMIX_examples.RData")

  expect_equal(nlpsem:::.auto_names_LGCM("bilinear spline", TRUE, c("ex1", "ex2")), eL$Math_LGCM_TIC_BLS_f@Estimates$Name)
  expect_equal(nlpsem:::.auto_names_LCSM("nonparametric", FALSE, NULL, 9), eC$Read_LCSM_NonP@Estimates$Name)
  expect_equal(nlpsem:::.auto_names_LCSM("nonparametric", FALSE, c("ex1", "ex2"), 9), eC$Read_LCSM_NonP_TIC@Estimates$Name)
  expect_equal(nlpsem:::.auto_names_LCSM("quadratic", FALSE, NULL, 9), eC$Read_LCSM_QUAD@Estimates$Name)
  expect_equal(nlpsem:::.auto_names_LCSM("negative exponential", FALSE, NULL, 9), eC$Read_LCSM_EXP_r@Estimates$Name)
  expect_equal(nlpsem:::.auto_names_LCSM("Jenss-Bayley", FALSE, NULL, 9), eC$Read_LCSM_JB_r@Estimates$Name)
  expect_equal(nlpsem:::.auto_names_TVC("BLS", TRUE, 1:9, "LGCM", 1, "ex1"), eT$Math_TVCslp_BLS_f@Estimates$Name)
  expect_equal(nlpsem:::.auto_names_TVC("BLS", FALSE, 1:9, "LGCM", 2, "ex1"), eT$Math_TVCchg_BLS_r@Estimates$Name)
  expect_equal(nlpsem:::.auto_names_TVC("BLS", FALSE, 1:9, "LGCM", 3, "ex1"), eT$Math_TVCchgBL_BLS_r@Estimates$Name)
  expect_equal(nlpsem:::.auto_names_MGM("BLS", TRUE, c("R", "M"), list(1:9, 1:9), "LGCM"), eG$RM_PLGCM.f@Estimates$Name)
  expect_equal(nlpsem:::.auto_names_MED("bilinear spline", "baseline", "ex1", "R", "M"), eD$Med2_LGCM_BLS@Estimates$Name)
  expect_equal(nlpsem:::.auto_names_MED("bilinear spline", "longitudinal", "R", "M", "S"), eD$Med3_LGCM_BLS@Estimates$Name)

  # Class-specific fits: the auto names feed the output helpers, which add class prefixes
  # (and the cluster-TIC coefficients for mixture models).
  grp <- eR$MGroup_Read_EXP_LGCM_TIC_r
  nms <- nlpsem:::.auto_names_subModel("LGCM", "EXP", FALSE, c("ex1", "ex2"), 1:9, NULL, NULL, "R", NULL, NULL, NULL)
  expect_equal((nlpsem:::getMGroup.output(grp@mxOutput, nClass = 2, sub_Model = "LGCM", y_var = "R", curveFun = "EXP",
                                         x_type = NULL, records = 1:9, growth_TIC = c("ex1", "ex2"), y_model = NULL,
                                         decompose = NULL, names = nms)), grp@Estimates)
  # RM_BLS_PLGCM3 (three classes with cluster TICs): the whole stored table, including the cluster-TIC rows,
  # must be reproduced, and every cluster-TIC label must carry the value the model stores under that label
  # (0.4 listed the labels class by class but filled the values column by column).
  mix <- eM$RM_BLS_PLGCM3
  nms <- nlpsem:::.auto_names_subModel("MGM", "BLS", FALSE, NULL, list(1:9, 1:9), "LGCM", NULL, c("R", "M"), NULL, NULL, NULL)
  expect_equal(3 * length(nms) + 3 * 2, nrow(mix@Estimates))
  expect_equal(mix@Estimates$Name[seq_along(nms)], paste0("c1", nms))
  out <- nlpsem:::getMIX.output(mix@mxOutput, nClass = 3, sub_Model = "MGM",
                                cluster_TIC = c("gx1", "gx2"), y_var = c("R", "M"), curveFun = "BLS",
                                x_type = NULL, records = list(1:9, 1:9), growth_TIC = NULL,
                                y_model = "LGCM", decompose = NULL, names = nms)
  expect_equal((out), mix@Estimates)
  cb <- mix@mxOutput$classbeta
  beta_rows <- out[-seq_len(3 * length(nms)), ]
  expect_equal(beta_rows$Name, c(t(cb$labels[-1, , drop = FALSE])))
  for (lab in beta_rows$Name) {
    expect_equal(beta_rows$Estimate[beta_rows$Name == lab], cb$values[cb$labels == lab], info = lab)
  }
})

test_that("getIndFS reduces the data to the model variables and keeps regression scores unchanged", {
  eL <- .load_extdata("getLGCM_examples.RData")
  m <- eL$Math_LGCM_BLS_r@mxOutput
  used <- intersect(nlpsem:::.model_variables(m), colnames(m$data$observed))
  expect_gt(ncol(m$data$observed), length(used))
  prepared <- nlpsem:::.prepare_FS_model(m)
  expect_false(grepl("[^A-Za-z0-9_]", prepared@name))
  expect_setequal(colnames(prepared$data$observed), used)

  fs <- suppressMessages(getIndFS(m, FS_type = "Regression"))
  expect_s4_class(fs, "FSOutput")
  expect_equal(dim(fs@scores_est), c(500, 3))
  expect_true(all(is.finite(fs@scores_est)))
  direct <- suppressMessages(OpenMx::mxFactorScores(OpenMx::mxRename(m, "direct"), type = "Regression"))
  expect_equal(unname(fs@scores_est), unname(direct[, , 1]))
})

test_that("getIndFS supports ML and weighted ML scores on models whose data carry unused columns", {
  skip_on_cran()   # each scoring run re-optimises one model per row (roughly 20 s each)
  eL <- .load_extdata("getLGCM_examples.RData")
  m <- eL$Math_LGCM_BLS_r@mxOutput
  for (type in c("ML", "WeightedML")) {
    fs <- suppressWarnings(suppressMessages(getIndFS(m, FS_type = type)))
    expect_s4_class(fs, "FSOutput")
    expect_equal(dim(fs@scores_est), c(500, 3), info = type)
    expect_true(all(is.finite(fs@scores_est)), info = type)
  }
})

test_that("getLRT compares nested fits", {
  eL <- .load_extdata("getLGCM_examples.RData")
  full <- eL$Math_LGCM_BLS_f@mxOutput
  reduced <- eL$Math_LGCM_BLS_r@mxOutput
  lrt <- getLRT(full, reduced)
  expect_equal(rownames(lrt), c("Full Model", "Reduced Model"))
  expect_equal(lrt[["Diff in DoF"]][2], 4)
  expect_type(lrt[["p.values"]], "double")
  expect_lt(lrt[["p.values"]][2], 1e-4)
  expect_equal(lrt$BIC, c(BIC(full), BIC(reduced)))
  expect_error(getLRT(full, reduced, boot = TRUE), "bootstrap replications")
})

test_that("conf.level is passed to likelihood-based intervals", {
  expect_equal(nlpsem:::.ci_probs(0.90), c(0.05, 0.95))
  expect_error(nlpsem:::.ci_probs(1.2), "between 0 and 1")

  captured <- NULL
  local_mocked_bindings(
    mxTryHard = function(model, ...) {
      captured <<- model
      list(output = list(confidenceIntervals = matrix(c(0.8, 1, 1.2), nrow = 1,
                                                     dimnames = list("a", c("lbound", "estimate", "ubound")))))
    },
    .package = "nlpsem"
  )
  stub <- OpenMx::mxModel("stub", OpenMx::mxMatrix("Full", 1, 1, free = TRUE, values = 1, labels = "a", name = "A"))
  est <- data.frame(Name = "a", Estimate = 1, SE = 0.1)
  out <- getEstimateStats(model = stub, est_in = est, CI_type = "likelihood", conf.level = 0.90)
  expect_equal(captured@intervals[[1]]@lowerdelta, qchisq(0.90, 1))
  expect_equal(out@likelihood$lik_lbound, 0.8)
  getEstimateStats(model = stub, est_in = est, CI_type = "likelihood")
  expect_equal(captured@intervals[[1]]@lowerdelta, qchisq(0.95, 1))
})

test_that("conf.level changes the width of bootstrap intervals", {
  skip_on_cran()
  dat <- .rebaseline_rms()
  fit <- suppressMessages(getLGCM(dat, t_var = "T", y_var = "M", curveFun = "LIN", intrinsic = FALSE,
                                  records = 1:9, paramOut = TRUE))
  set.seed(20191029)
  b80 <- suppressWarnings(suppressMessages(getEstimateStats(model = fit@mxOutput, est_in = fit@Estimates,
                                                            CI_type = "bootstrap", rep = 20, conf.level = 0.80)))
  set.seed(20191029)
  b99 <- suppressWarnings(suppressMessages(getEstimateStats(model = fit@mxOutput, est_in = fit@Estimates,
                                                            CI_type = "bootstrap", rep = 20, conf.level = 0.99)))
  w80 <- b80@bootstrap$boot_ubound - b80@bootstrap$boot_lbound
  w99 <- b99@bootstrap$boot_ubound - b99@bootstrap$boot_lbound
  expect_true(all(w99 >= w80))
  expect_true(any(w99 > w80))
})

test_that("getSummary returns a plain data frame for heterogeneous model lists", {
  eM <- .load_extdata("getMIX_examples.RData")
  eR <- .load_extdata("getMGroup_examples.RData")
  s <- suppressWarnings(getSummary(list(eM$Math_BLS_LGCM1@mxOutput, eM$Math_BLS_LGCM2@mxOutput,
                                        eM$Math_BLS_LGCM3@mxOutput), HetModels = TRUE))
  expect_s3_class(s, "data.frame")
  expect_false(inherits(s, "tbl_df"))
  expect_equal(nrow(s), 3)
  expect_true(all(c("Model", "No_Params", "-2ll", "AIC", "BIC", "Y_res_c1", "Y_res_c3", "%Class1", "%Class3") %in% names(s)))
  expect_equal(s[["%Class1"]][1], "100%")
  g <- suppressWarnings(getSummary(list(eR$MGroup_Math_BLS_LGCM_f@mxOutput), HetModels = TRUE))
  expect_false(inherits(g, "tbl_df"))
  expect_equal(nrow(g), 1)
  expect_true(all(c("%Class1", "%Class2") %in% names(g)))
})

test_that("getPosterior handles cluster-forming covariates", {
  eM <- .load_extdata("getMIX_examples.RData")
  post <- getPosterior(eM$RM_BLS_PLGCM3@mxOutput, nClass = 3, cluster_TIC = c("gx1", "gx2"))
  expect_equal(dim(post@prob), c(500, 3))
  expect_equal(unname(rowSums(post@prob)), rep(1, 500), tolerance = 1e-8)
  expect_true(all(post@membership %in% 1:3))
})
