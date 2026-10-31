# getSimData(): structure, generating distribution, and recovery of the true values by the matching fit.

test_that("simulated data have the columns and true values of the growth curve and change score fits", {
  specs <- list(c("LIN", FALSE), c("QUAD", FALSE), c("EXP", TRUE), c("JB", FALSE), c("JB", TRUE), c("BLS", TRUE), c("NonP", FALSE))
  for (s in specs) {
    intrinsic <- as.logical(s[2])
    dat <- getSimData(n = 30, curveFun = s[1], intrinsic = intrinsic, seed = 1)
    expect_identical(names(dat), c("id", paste0("Y", 1:10), paste0("T", 1:10)), info = paste(s, collapse = " "))
    truth <- attr(dat, "truth")
    # every free parameter of each fit that supports the form, without the derived change-score quantities
    if (s[1] != "NonP") expect_identical(names(truth), nlpsem:::.auto_names_LGCM(s[1], intrinsic, NULL), info = paste(s, collapse = " "))
    if (!s[1] %in% c("LIN", "BLS")) {
      lcsm <- nlpsem:::.auto_names_LCSM(s[1], intrinsic, NULL, 10)
      expect_identical(names(truth), lcsm[!grepl("_est[0-9]*$", lcsm)], info = paste(s, collapse = " "))
    }
    expect_equal(dat$T1, rep(0, 30)); expect_equal(dat$T10, rep(9, 30))
    expect_true(all(abs(as.matrix(dat[paste0("T", 2:9)]) - matrix(1:8, 30, 8, byrow = TRUE)) <= 0.25))
  }
})

test_that("a seed makes the data reproducible without changing the session's random numbers", {
  set.seed(5); before <- runif(1)
  set.seed(5); a <- getSimData(n = 20, curveFun = "LIN", seed = 11); after <- runif(1)
  expect_identical(before, after)
  expect_identical(a, getSimData(n = 20, curveFun = "LIN", seed = 11))
})

test_that("the generating distribution has the requested moments", {
  p <- list(mean = c(eta0 = 10, eta1 = 2, eta2 = -0.5, gamma = 3), cov = diag(c(4, 1, 0.25, 0.04)), residual = 0.5)
  dat <- getSimData(n = 20000, curveFun = "BLS", intrinsic = TRUE, params = p, jitter = 0, seed = 2)
  gf <- attr(dat, "growth_factors")
  expect_equal(colMeans(gf), p$mean, tolerance = 0.02)
  expect_equal(unname(diag(cov(gf))), diag(p$cov), tolerance = 0.05)
  curve <- gf[, "eta0"] + gf[, "eta1"] * pmin(as.matrix(dat[paste0("T", 1:10)]), gf[, "gamma"]) +
    gf[, "eta2"] * pmax(as.matrix(dat[paste0("T", 1:10)]) - gf[, "gamma"], 0)
  expect_equal(var(c(as.matrix(dat[paste0("Y", 1:10)]) - curve)), 0.5, tolerance = 0.02)
  # growth TICs: marginal means are the given means, conditional covariance the given cov
  q <- list(mean = c(eta0 = 10, eta1 = 2), cov = diag(c(4, 1)), residual = 1, tic_mean = 1, tic_cov = matrix(2),
            beta = matrix(c(1, 0.5), 2))
  dat <- getSimData(n = 20000, curveFun = "LIN", params = q, seed = 3)
  gf <- attr(dat, "growth_factors")
  expect_equal(colMeans(gf), q$mean, tolerance = 0.02)
  expect_equal(unname(coef(lm(gf[, "eta1"] ~ dat$ex1))[2]), 0.5, tolerance = 0.03)
  truth <- attr(dat, "truth")
  expect_equal(truth[c("alpha0", "alpha1", "beta10", "beta11", "mux1", "phi11", "mueta0", "mueta1")],
               c(alpha0 = 9, alpha1 = 1.5, beta10 = 1, beta11 = 0.5, mux1 = 1, phi11 = 2, mueta0 = 10, mueta1 = 2))
})

test_that("classes follow the proportions or the cluster-TIC logits", {
  p <- list(mean = c(eta0 = 10, eta1 = 2), cov = diag(2), residual = 1)
  dat <- getSimData(n = 20000, curveFun = "LIN", params = list(p, p, p), prop = c(0.2, 0.3, 0.5), seed = 4)
  expect_equal(as.numeric(table(dat$class)) / 20000, c(0.2, 0.3, 0.5), tolerance = 0.03)
  truth <- attr(dat, "truth")
  expect_equal(truth[c("p2", "p3")], c(p2 = log(1.5), p3 = log(2.5)))
  expect_true(all(c("c1mueta0", "c3residuals") %in% names(truth)))
  beta <- rbind(0, c(-0.5, 1, 0))
  dat <- getSimData(n = 20000, curveFun = "LIN", params = list(p, p), cluster = list(mean = c(0, 0), cov = diag(2), beta = beta), seed = 5)
  logit <- glm(I(class == 2) ~ gx1 + gx2, family = stats::binomial, data = dat)
  expect_equal(unname(coef(logit)), c(-0.5, 1, 0), tolerance = 0.06)
  expect_equal(attr(dat, "truth")[c("beta20", "beta21", "beta22")], c(beta20 = -0.5, beta21 = 1, beta22 = 0))
})

test_that("inputs that do not describe a model are rejected", {
  expect_error(getSimData(sub_Model = "LCSM", curveFun = "QUAD"), "latent change score model can also fit")
  expect_error(getSimData(curveFun = "LIN", intrinsic = TRUE), "intrinsic")
  expect_error(getSimData(curveFun = "LIN", params = list(mean = c(eta0 = 1, eta1 = 1), cov = diag(3), residual = 1)), "2 x 2")
  expect_error(getSimData(curveFun = "LIN", params = list(mean = c(eta0 = 1, eta1 = 1), cov = matrix(c(1, 2, 2, 1), 2), residual = 1)),
               "positive definite")
  expect_error(getSimData(curveFun = "LIN", jitter = 0.6), "jitter")
  p <- list(mean = c(eta0 = 10, eta1 = 2), cov = diag(2), residual = 1)
  expect_error(getSimData(curveFun = "LIN", params = list(p, p), cluster = list(mean = 0, cov = 1, beta = rbind(1, c(0, 1)))), "first row of zeros")
})

test_that("fits of the matching model recover the true values", {
  skip_on_cran()
  OpenMx::mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)
  recover <- function(dat, fit) {
    truth <- attr(dat, "truth"); est <- fit@Estimates
    i <- match(names(truth), est$Name)
    expect_false(anyNA(i))
    expect_equal(fit@mxOutput@output$status$code, 0)
    expect_length(getDiagnostics(fit)$issues, 0)
    expect_lt(max(abs((est$Estimate[i] - truth) / est$SE[i])), 4.5)
  }
  dat <- getSimData(n = 500, curveFun = "EXP", intrinsic = TRUE, seed = 20261027)
  set.seed(20191029)
  recover(dat, suppressMessages(getLGCM(dat = dat, t_var = "T", y_var = "Y", curveFun = "EXP", intrinsic = TRUE, records = 1:10,
                                        paramOut = TRUE)))
  # one data set, fitted as a growth curve and as change scores
  dat <- getSimData(n = 500, curveFun = "QUAD", seed = 20261027)
  set.seed(20191029)
  recover(dat, suppressMessages(getLGCM(dat = dat, t_var = "T", y_var = "Y", curveFun = "QUAD", intrinsic = FALSE, records = 1:10,
                                        paramOut = TRUE)))
  set.seed(20191029)
  recover(dat, suppressMessages(getLCSM(dat = dat, t_var = "T", y_var = "Y", curveFun = "QUAD", intrinsic = FALSE, records = 1:10,
                                        paramOut = TRUE)))
  dat <- getSimData(n = 500, curveFun = "NonP", seed = 20261027)
  set.seed(20191029)
  recover(dat, suppressMessages(getLCSM(dat = dat, t_var = "T", y_var = "Y", curveFun = "NonP", intrinsic = FALSE, records = 1:10,
                                        paramOut = TRUE)))
  p1 <- list(mean = c(eta0 = 40, eta1 = 20, eta2 = -1), cov = diag(c(25, 1, 0.09)), residual = 1)
  p2 <- p1; p2$mean[c("eta0", "eta2")] <- c(60, -1.5)
  dat <- getSimData(n = 500, curveFun = "QUAD", params = list(p1, p2), prop = c(0.4, 0.6), seed = 20261027)
  set.seed(20191029)
  fit <- suppressMessages(getMIX(dat = dat, prop_starts = c(0.5, 0.5), sub_Model = "LGCM", t_var = "T", y_var = "Y", curveFun = "QUAD",
                                 intrinsic = FALSE, records = 1:10, paramOut = TRUE))
  recover(dat, fit)
  membership <- getPosterior(fit@mxOutput, nClass = 2, cluster_TIC = NULL)@membership
  expect_gt(mean(membership == dat$class), 0.95)
  dat <- getSimData(n = 500, sub_Model = "MGM", curveFun = "LIN", seed = 20261027)
  set.seed(20191029)
  recover(dat, suppressMessages(getMGM(dat = dat, t_var = c("T", "T"), y_var = c("Y", "Z"), curveFun = "LIN", intrinsic = FALSE,
                                       records = list(1:10, 1:10), y_model = "LGCM", paramOut = TRUE)))
  dat <- getSimData(n = 500, sub_Model = "TVC", curveFun = "LIN", decompose = 1, nT = 5, seed = 20261027)
  set.seed(20191029)
  recover(dat, suppressMessages(getTVCmodel(dat = dat, t_var = "T", y_var = "Y", curveFun = "LIN", intrinsic = FALSE, records = 1:5,
                                            y_model = "LGCM", TVC = "TVC", decompose = 1, paramOut = TRUE)))
  dat <- getSimData(n = 500, sub_Model = "MED", curveFun = "LIN", x_type = "baseline", seed = 20261027)
  set.seed(20191029)
  recover(dat, suppressMessages(getMediation(dat = dat, t_var = c("T", "T"), y_var = "Y", m_var = "M", x_type = "baseline",
                                             x_var = "X", curveFun = "LIN", records = list(1:10, 1:10), paramOut = TRUE)))
})

test_that("multivariate, TVC and mediation data carry the true values of the matching fit", {
  configs <- list(
    list(sub_Model = "MGM", curveFun = "LIN"), list(sub_Model = "MGM", curveFun = "BLS", intrinsic = TRUE, y_var = c("Y", "Z", "W")),
    list(sub_Model = "MGM", curveFun = "NonP"), list(sub_Model = "MGM", curveFun = "JB", intrinsic = TRUE),
    list(sub_Model = "TVC", decompose = 0, curveFun = "EXP"), list(sub_Model = "TVC", decompose = 1, curveFun = "BLS", intrinsic = TRUE),
    list(sub_Model = "TVC", decompose = 2, curveFun = "QUAD"), list(sub_Model = "TVC", decompose = 3, curveFun = "NonP"),
    list(sub_Model = "MED", curveFun = "LIN", x_type = "baseline"), list(sub_Model = "MED", curveFun = "LIN", x_type = "longitudinal"),
    list(sub_Model = "MED", curveFun = "BLS", x_type = "baseline"), list(sub_Model = "MED", curveFun = "BLS", x_type = "longitudinal"))
  for (cfg in configs) {
    dat <- do.call(getSimData, c(list(n = 30, seed = 1), cfg))
    s <- attr(dat, "settings")
    ctx <- list(sub_Model = s$sub_Model, curveFun = s$curveFun, intrinsic = s$intrinsic, nT = s$nT,
                y_model = if (s$curveFun == "NonP") "LCSM" else "LGCM",
                decompose = s$decompose, x_type = s$x_type, y_var = s$y_var, growth_TIC = s$growth_TIC)
    label <- paste(unlist(cfg), collapse = " ")
    expect_identical(names(attr(dat, "truth")), nlpsem:::.sim_fit_names(ctx), info = label)
    expect_false(anyNA(dat), info = label)
  }
  # growth TICs in a TVC model with a decomposed TVC
  p <- nlpsem:::.sim_default_params(list(sub_Model = "TVC", form = "LIN", intrinsic = FALSE, nT = 10, decompose = 1))
  p$tic_mean <- c(0, 0); p$tic_cov <- diag(2); p$beta <- matrix(0.5, 2, 2); p$cov_bl <- c(0.3, 0.2)
  truth <- attr(getSimData(n = 30, sub_Model = "TVC", curveFun = "LIN", decompose = 1, params = p, seed = 1), "truth")
  expect_equal(truth[c("betaTIC10", "covBL1", "covBL2", "muTIC2", "phiTIC12")], c(betaTIC10 = 0.5, covBL1 = 0.3, covBL2 = 0.2, muTIC2 = 0, phiTIC12 = 0))
})

test_that("growth curve and change score fits share the names of the true values of every parametric form", {
  for (form in c("QUAD", "EXP", "JB")) for (intrinsic in c(FALSE, TRUE)) {
    if (form == "QUAD" && intrinsic) next
    names_of <- function(...) nlpsem:::.sim_fit_names(list(curveFun = form, intrinsic = intrinsic, nT = 10, y_var = c("Y", "Z"), ...))
    label <- paste(form, intrinsic)
    for (tic in list(NULL, c("ex1", "ex2"))) {
      expect_identical(names_of(sub_Model = "LGCM", y_model = "LGCM", growth_TIC = tic), names_of(sub_Model = "LGCM", y_model = "LCSM", growth_TIC = tic),
                       info = label)
      for (decompose in 0:3) {
        expect_identical(names_of(sub_Model = "TVC", y_model = "LGCM", decompose = decompose, growth_TIC = tic),
                         names_of(sub_Model = "TVC", y_model = "LCSM", decompose = decompose, growth_TIC = tic), info = label)
      }
    }
    expect_identical(names_of(sub_Model = "MGM", y_model = "LGCM"), names_of(sub_Model = "MGM", y_model = "LCSM"), info = label)
  }
})

test_that("multivariate data have the requested cross-process covariances", {
  dat <- getSimData(n = 20000, sub_Model = "MGM", curveFun = "LIN", jitter = 0, seed = 6)
  p <- attr(dat, "settings")$params[[1]]
  gf <- attr(dat, "growth_factors")
  expect_equal(unname(cov(gf)), unname(p$cov), tolerance = 0.05)
  eY <- as.matrix(dat[paste0("Y", 1:10)]) - (gf[, "Y_eta0"] + gf[, "Y_eta1"] * as.matrix(dat[paste0("T", 1:10)]))
  eZ <- as.matrix(dat[paste0("Z", 1:10)]) - (gf[, "Z_eta0"] + gf[, "Z_eta1"] * as.matrix(dat[paste0("T", 1:10)]))
  expect_equal(cov(c(eY), c(eZ)), p$residual[1, 2], tolerance = 0.05)
  expect_equal(attr(dat, "truth")[c("YZ_psi10", "YZ_res")], c(YZ_psi10 = p$cov[2, 3], YZ_res = 0.3))
})

test_that("a decomposed TVC acts through its slopes, changes or change from baseline", {
  for (decompose in 1:3) {
    dat <- getSimData(n = 20000, sub_Model = "TVC", curveFun = "LIN", decompose = decompose, seed = 7)
    p <- attr(dat, "settings")$params[[1]]
    gf <- attr(dat, "growth_factors"); tm <- as.matrix(dat[paste0("T", 1:10)])
    slope <- outer(gf[, "TVC_eta1"], p$tvc$rel_rate); lag <- t(apply(tm, 1, diff))
    true_x <- gf[, "TVC_eta0"] + cbind(0, t(apply(slope * lag, 1, cumsum)))
    state <- switch(decompose, cbind(0, slope), cbind(0, slope * lag), true_x - gf[, "TVC_eta0"])
    eY <- as.matrix(dat[paste0("Y", 1:10)]) - (gf[, "Y_eta0"] + gf[, "Y_eta1"] * tm) - p$kappa * state
    eX <- as.matrix(dat[paste0("TVC", 1:10)]) - true_x
    expect_equal(c(var(c(eY)), var(c(eX)), cov(c(eY), c(eX))), c(p$residual, p$tvc$residual, p$res_cov), tolerance = 0.05, info = paste("decompose", decompose))

    expect_equal(unname(coef(lm(gf[, "Y_eta0"] ~ gf[, "TVC_eta0"]))[2]), p$beta_tvc[1], tolerance = 0.1, info = paste("decompose", decompose))
  }
})

test_that("mediation paths are the requested effects", {
  dat <- getSimData(n = 20000, sub_Model = "MED", curveFun = "LIN", x_type = "baseline", seed = 8)
  p <- attr(dat, "settings")$params[[1]]; gf <- attr(dat, "growth_factors")
  expect_equal(unname(coef(lm(gf[, "M_eta1"] ~ dat$X))[2]), p$beta_xm[2, 1], tolerance = 0.05)
  fit <- lm(gf[, "Y_eta1"] ~ dat$X + gf[, "M_eta0"] + gf[, "M_eta1"])
  expect_equal(unname(coef(fit)[-1]), c(p$beta_xy[2, 1], p$beta_my[2, ]), tolerance = 0.05)
  truth <- attr(dat, "truth")
  expect_equal(truth[c("betaM0Y1", "alphaM0")], c(betaM0Y1 = p$beta_my[2, 1], alphaM0 = p$m$mean[["eta0"]] - p$beta_xm[1, 1] * p$x$mean))
})

test_that("multivariate, TVC and mediation inputs that do not describe a model are rejected", {
  expect_error(getSimData(sub_Model = "MED", curveFun = "QUAD", x_type = "baseline"), "linear or bilinear")
  expect_error(getSimData(sub_Model = "MED", curveFun = "LIN"), "x_type")
  expect_error(getSimData(sub_Model = "TVC", curveFun = "LIN"), "decompose")
  expect_error(getSimData(sub_Model = "MGM", curveFun = "LIN", y_var = "Y"), "two or more")
  p <- nlpsem:::.sim_default_params(list(sub_Model = "MED", form = "LIN", x_type = "longitudinal"))
  p$beta_my[1, 2] <- 0.1
  expect_error(getSimData(sub_Model = "MED", curveFun = "LIN", x_type = "longitudinal", params = p), "later ones")
})

test_that("named parameters keep their meaning in any order", {
  p <- nlpsem:::.sim_default_params(list(sub_Model = "TVC", form = "LIN", intrinsic = FALSE, nT = 10, decompose = 1))
  p$tvc$mean <- c(eta0 = 20, eta1 = 2)
  q <- p; q$tvc$mean <- rev(q$tvc$mean); q$mean <- rev(q$mean)
  a <- getSimData(n = 200, sub_Model = "TVC", curveFun = "LIN", decompose = 1, params = p, seed = 9)
  b <- getSimData(n = 200, sub_Model = "TVC", curveFun = "LIN", decompose = 1, params = q, seed = 9)
  expect_identical(as.data.frame(a), as.data.frame(b))
  expect_identical(attr(a, "growth_factors"), attr(b, "growth_factors"))
  m <- nlpsem:::.sim_default_params(list(sub_Model = "MED", form = "BLS", x_type = "longitudinal"))
  r <- m; r$x$mean <- rev(r$x$mean); r$m$mean <- rev(r$m$mean); r$y$mean <- rev(r$y$mean)
  expect_identical(as.data.frame(getSimData(n = 50, sub_Model = "MED", curveFun = "BLS", x_type = "longitudinal", params = m, seed = 2)),
                   as.data.frame(getSimData(n = 50, sub_Model = "MED", curveFun = "BLS", x_type = "longitudinal", params = r, seed = 2)))
})

test_that("column prefixes do not collide with the roles of the simulated processes", {
  stems <- function(d) unique(sub("[0-9]+$", "", names(d)))
  med <- getSimData(n = 20, sub_Model = "MED", curveFun = "LIN", x_type = "baseline", y_var = "M", m_var = "Mediator",
                    x_var = "Predictor", seed = 1)
  expect_identical(stems(med), c("id", "Mediator", "M", "T", "Predictor"))
  ref <- getSimData(n = 20, sub_Model = "MED", curveFun = "LIN", x_type = "baseline", seed = 1)
  expect_identical(unname(as.matrix(med[paste0("M", 1:10)])), unname(as.matrix(ref[paste0("Y", 1:10)])))
  expect_identical(unname(as.matrix(med[paste0("Mediator", 1:10)])), unname(as.matrix(ref[paste0("M", 1:10)])))
  tvc <- getSimData(n = 20, sub_Model = "TVC", curveFun = "LIN", decompose = 1, y_var = "TVC", TVC = "Covariate", seed = 1)
  expect_identical(stems(tvc), c("id", "TVC", "Covariate", "T"))
  expect_identical(stems(getSimData(n = 20, curveFun = "LIN", y_var = "TIC", seed = 1)), c("id", "TIC", "T"))
  expect_identical(stems(getSimData(n = 20, sub_Model = "MGM", curveFun = "LIN", y_var = c("Y", "TIC"), seed = 1)), c("id", "Y", "TIC", "T"))
  expect_error(getSimData(n = 20, curveFun = "LIN", y_var = "T"), "duplicated columns: T1")
  p <- nlpsem:::.sim_default_params(list(sub_Model = "LGCM", form = "LIN", intrinsic = FALSE, nT = 10))
  p$tic_mean <- 0; p$tic_cov <- 1; p$beta <- matrix(0.5, 2, 1)
  expect_error(getSimData(n = 20, curveFun = "LIN", params = p, growth_TIC = "Y3"), "duplicated columns: Y3")
})

test_that("malformed parameter values are rejected before any data are drawn", {
  p <- nlpsem:::.sim_default_params(list(sub_Model = "LGCM", form = "QUAD", intrinsic = FALSE, nT = 10))
  bad <- function(change, pattern) { q <- p; q <- change(q); expect_error(getSimData(n = 20, curveFun = "QUAD", params = q), pattern) }
  bad(function(q) { q$cov[1, 2] <- 3; q }, "'cov' must be symmetric")
  bad(function(q) { q$cov <- diag(c(1, -1, 1)); q }, "'cov' must be positive definite")
  bad(function(q) { q$mean[["eta1"]] <- Inf; q }, "'mean' must be finite")
  bad(function(q) { q$mean[["eta1"]] <- NA; q }, "'mean' must be finite")
  bad(function(q) { q$mean <- c(q$mean, eta0 = 1); q }, "once each")
  bad(function(q) { q$mean <- c(q$mean, eta3 = 1); q }, "once each")
  bad(function(q) { q$residual <- -1; q }, "'residual' must be positive")
  m <- nlpsem:::.sim_default_params(list(sub_Model = "MED", form = "LIN", x_type = "baseline"))
  m$beta_xm <- 0.7
  expect_error(getSimData(n = 20, sub_Model = "MED", curveFun = "LIN", x_type = "baseline", params = m), "'beta_xm' must give 2 effects")
  t <- nlpsem:::.sim_default_params(list(sub_Model = "TVC", form = "LIN", intrinsic = FALSE, nT = 10, decompose = 1))
  t$beta_tvc <- 0.2
  expect_error(getSimData(n = 20, sub_Model = "TVC", curveFun = "LIN", decompose = 1, params = t), "'beta_tvc' must be 2 finite")
  l <- nlpsem:::.sim_default_params(list(sub_Model = "LGCM", form = "LIN", intrinsic = FALSE, nT = 10))
  expect_error(getSimData(n = 20, curveFun = "LIN", params = list(l, l), cluster = list(mean = 0, cov = -1, beta = rbind(0, c(0, 1)))),
               "'cluster\\$cov' must be a symmetric positive-definite")
})

test_that("simulation names cannot be interpreted as data-frame column positions", {
  p <- nlpsem:::.sim_default_params(list(sub_Model = "LGCM", form = "LIN", intrinsic = FALSE, nT = 10))
  p$tic_mean <- 0; p$tic_cov <- 1; p$beta <- matrix(0.5, 2, 1)
  for (bad in list(1, NA_character_, "")) {
    expect_error(getSimData(n = 20, curveFun = "LIN", params = p, growth_TIC = bad), "growth_TIC.*character column names")
    expect_error(getSimData(n = 20, curveFun = "LIN", params = list(p, p),
                            cluster = list(mean = 0, cov = matrix(1), beta = rbind(c(0, 0), c(0, 1))),
                            cluster_TIC = bad), "cluster_TIC.*character column names")
  }
  d <- getSimData(n = 20, curveFun = "LIN", params = p, growth_TIC = "1", seed = 7)
  expect_identical(d$id, 1:20)
  expect_true("1" %in% names(d))
  ref <- getSimData(n = 20, sub_Model = "TVC", curveFun = "LIN", decompose = 1, seed = 7)
  named <- getSimData(n = 20, sub_Model = "TVC", curveFun = "LIN", decompose = 1,
                      y_var = c(outcome = "Y"), TVC = c(covariate = "TVC"), seed = 7)
  expect_identical(named, ref)
})

test_that("nonfinite simulation sizes and schedules are rejected", {
  for (bad in c(NA_real_, NaN, Inf, -Inf)) {
    expect_error(getSimData(curveFun = "LIN", n = bad), "'n'.*finite integer")
    expect_error(getSimData(curveFun = "LIN", nT = bad), "'nT'.*finite integer")
    expect_error(getSimData(curveFun = "LIN", times = c(0:8, bad), jitter = 0), "'times'.*finite increasing")
    expect_error(getSimData(curveFun = "LIN", jitter = bad), "'jitter'")
  }
})
