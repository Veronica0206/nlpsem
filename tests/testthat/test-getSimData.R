# getSimData(): structure, generating distribution, and recovery of the true values by the matching fit.

test_that("simulated data have the columns and true values of the matching fit", {
  specs <- list(c("LGCM", "LIN", FALSE), c("LGCM", "QUAD", FALSE), c("LGCM", "EXP", TRUE), c("LGCM", "JB", FALSE),
                c("LGCM", "BLS", TRUE), c("LCSM", "QUAD", FALSE), c("LCSM", "JB", TRUE), c("LCSM", "NonP", FALSE))
  for (s in specs) {
    intrinsic <- as.logical(s[3])
    dat <- getSimData(n = 30, sub_Model = s[1], curveFun = s[2], intrinsic = intrinsic, seed = 1)
    expect_identical(names(dat), c("id", paste0("Y", 1:10), paste0("T", 1:10)), info = paste(s, collapse = " "))
    fit_names <- if (s[1] == "LGCM") nlpsem:::.auto_names_LGCM(s[2], intrinsic, NULL) else
      nlpsem:::.auto_names_LCSM(s[2], intrinsic, NULL, 10)
    truth <- attr(dat, "truth")
    expect_true(all(names(truth) %in% fit_names), info = paste(s, collapse = " "))
    # every free parameter of the fit, without the derived change-score quantities
    expect_identical(names(truth), fit_names[!grepl("_est[0-9]*$", fit_names)], info = paste(s, collapse = " "))
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
  expect_error(getSimData(curveFun = "NonP"), "sub_Model = 'LCSM'")
  expect_error(getSimData(sub_Model = "LCSM", curveFun = "LIN"), "LCSM does not support")
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
  dat <- getSimData(n = 500, sub_Model = "LCSM", curveFun = "NonP", seed = 20261027)
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
})
