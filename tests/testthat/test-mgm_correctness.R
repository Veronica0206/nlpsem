# Regressions for covariance pairs, record alignment and class-specific output.
.mgm_correctness_data <- function(n = 240L, n_traj = 4L) {
  set.seed(914)
  b <- matrix(rnorm(n * n_traj * 2L), n, n_traj * 2L) + rnorm(n, sd = 0.3)
  dat <- as.data.frame(matrix(rep(0:3, each = n), n, 4L))
  names(dat) <- paste0("T", 1:4)
  for (i in seq_len(n_traj)) {
    for (j in 1:4) {
      dat[[paste0(LETTERS[i], j)]] <- 3 * i + b[, 2 * i - 1L] +
        (1 + i / 3 + b[, 2 * i] * 0.3) * (j - 1L) + rnorm(n, sd = 0.5)
    }
  }
  dat$group <- rep(1:2, length.out = n)
  dat
}

.mgm_test_optimizer <- function(expr) {
  old_optimizer <- OpenMx::mxOption(NULL, "Default optimizer")
  old_threads <- OpenMx::mxOption(NULL, "Number of Threads")
  on.exit({
    OpenMx::mxOption(NULL, "Default optimizer", old_optimizer)
    OpenMx::mxOption(NULL, "Number of Threads", old_threads)
  })
  OpenMx::mxOption(NULL, "Default optimizer", "CSOLNP")
  OpenMx::mxOption(NULL, "Number of Threads", 1L)
  eval.parent(substitute(expr))
}

test_that("four and five outcome models retain every distinct covariance pair", {
  for (n_traj in 4:5) {
    dat <- .mgm_correctness_data(80L, n_traj)
    outcomes <- LETTERS[seq_len(n_traj)]
    recs <- rep(list(1:4), n_traj)
    starts <- nlpsem:::getMULTI.initial(dat, rep("T", n_traj), outcomes,
                                       "LIN", recs, NULL, rep(0.2, choose(n_traj, 2)), FALSE)
    single <- nlpsem:::getMGM.mxModel(dat, rep("T", n_traj), outcomes,
                                     "LIN", FALSE, recs, "LGCM", starts)
    group <- nlpsem:::getsub.MGM(dat, 2L, rep("T", n_traj), outcomes,
                                "LIN", FALSE, recs, "LGCM", list(starts, starts),
                                "multigroup", "group")
    mixture <- nlpsem:::getsub.MGM(dat, 2L, rep("T", n_traj), outcomes,
                                  "LIN", FALSE, recs, "LGCM", list(starts, starts), "mixture")
    models <- c(list(single), group, mixture)
    prefixes <- c("", "c1", "c2", "c1", "c2")
    for (m_idx in seq_along(models)) {
      model <- models[[m_idx]]
      pfx <- prefixes[m_idx]
      # Two means and three within-process covariances per outcome, plus
      # four growth covariances and one residual covariance per outcome pair.
      expect_length(OpenMx::omxGetParameters(model), 6 * n_traj + 5 * choose(n_traj, 2))
      for (pair in combn(outcomes, 2, simplify = FALSE)) {
        from <- paste0(c("eta0", "eta1"), pair[1])
        to <- paste0(c("eta0", "eta1"), pair[2])
        expected <- matrix(paste0(pfx, paste0(pair, collapse = ""), "_psi", c("00", "01", "10", "11")),
                           2, 2, byrow = TRUE)
        expect_equal(unname(model$S$labels[from, to]), expected)
        expect_true(all(model$S$free[from, to]))
        expect_true(model$S$free[paste0(pair[1], "1"), paste0(pair[2], "1")])
        expect_equal(model$S$labels[paste0(pair[1], "1"), paste0(pair[2], "1")],
                     paste0(pfx, paste0(pair, collapse = ""), "_RES"))
        expect_true(paste0(pfx, paste0(pair, collapse = ""), "_psi") %in% names(model$algebras))
      }
    }
  }
})

test_that("public four-outcome MGM fits and reports independent between-process estimates", {
  .mgm_test_optimizer({
    dat <- .mgm_correctness_data()
    fit <- suppressMessages(getMGM(dat, rep("T", 4), LETTERS[1:4], "LIN", FALSE,
                                    rep(list(1:4), 4), "LGCM", paramOut = TRUE))
    expect_equal(fit@mxOutput@output$status$code, 0)
    expect_length(fit@mxOutput@output$estimate, 54)
    expect_equal(nrow(fit@Estimates), 54)
    expect_true(all(is.finite(fit@Estimates$Estimate)))
    expect_true(all(is.finite(fit@Estimates$SE)))
    for (pair in c("AB", "AC", "AD", "BC", "BD", "CD")) {
      actual <- fit@Estimates$Estimate[fit@Estimates$Name == paste0(pair, "_psi00")]
      expected <- OpenMx::mxEvalByName(paste0(pair, "_psi"), fit@mxOutput)[1, 1]
      expect_equal(unname(actual), unname(expected))
    }
    expect_false(isTRUE(all.equal(unname(fit@mxOutput@output$estimate["AD_psi00"]),
                                   unname(fit@mxOutput@output$estimate["BC_psi00"]))))
  })
})

test_that("MGM residual starts align common record IDs instead of column positions", {
  r1 <- cbind(c(1, -1, 1, -1), c(2, -2, 2, -2))
  r2 <- cbind(c(2, -2, 2, -2), c(-1, 1, -1, 1))
  expect_equal(nlpsem:::.residual_correlation_starts(list(r1, r2), records = list(1:2, 2:3)), 0.9)
  expect_warning(r <- nlpsem:::.residual_correlation_starts(list(r1, r2), records = list(1:2, 3:4)),
                 "falling back")
  expect_equal(r, 0.3)

  dat <- .mgm_correctness_data(80L, 2L)
  recs <- list(1:3, 2:4)
  residuals <- lapply(1:2, function(i) {
    gf <- nlpsem:::getUNI.GF(dat[, paste0(LETTERS[i], recs[[i]])],
                            dat[, paste0("T", recs[[i]])], 3L, "LIN")
    attr(gf, "residuals")
  })
  expected <- cor(as.vector(residuals[[1]][, 2:3]), as.vector(residuals[[2]][, 1:2]))
  expected <- max(-0.9, min(0.9, expected))
  starts <- nlpsem:::getMULTI.initial(dat, c("T", "T"), c("A", "B"), "LIN", recs,
                                     NULL, NULL, FALSE)
  actual <- starts$residuals[1, 2] / sqrt(prod(diag(starts$residuals)))
  expect_equal(actual, expected)
})

test_that("class-specific NonP MGM parameter tables contain correctly labeled estimates and SEs", {
  .mgm_test_optimizer({
    dat <- .mgm_correctness_data()
    fits <- list(group = suppressMessages(getMGroup(dat, "group", "MGM", rep("T", 2),
      rep(list(1:4), 2), LETTERS[1:2], "NonP", FALSE, y_model = "LCSM", paramOut = TRUE)))
    # Well separated classes make this a fast deterministic integration check of
    # the mixture output path rather than a difficult class-search benchmark.
    for (v in LETTERS[1:2]) for (j in 1:4) {
      dat[[paste0(v, j)]] <- dat[[paste0(v, j)]] + ifelse(dat$group == 2, 8, 0)
    }
    set.seed(914)
    fits$mixture <- suppressMessages(getMIX(dat, c(0.5, 0.5), "MGM", t_var = rep("T", 2),
      records = rep(list(1:4), 2), y_var = LETTERS[1:2], curveFun = "NonP", intrinsic = FALSE,
      y_model = "LCSM", paramOut = TRUE))
    for (fit in fits) {
      expect_equal(fit@mxOutput@output$status$code, 0)
      expect_gte(nrow(fit@Estimates), 114)
      expect_true(all(is.finite(fit@Estimates$Estimate)))
      expect_true(all(is.finite(fit@Estimates$SE)))
      for (k in 1:2) {
        prefix <- paste0("c", k)
        rate <- paste0(prefix, "A_rel_rate2")
        row <- fit@Estimates[fit@Estimates$Name == rate, ]
        expect_equal(nrow(row), 1L)
        expect_equal(row$Estimate, unname(fit@mxOutput@output$estimate[rate]))
        expect_equal(row$SE, unname(fit@mxOutput@output$standardErrors[rate, 1]))
        cross <- fit@Estimates$Estimate[fit@Estimates$Name == paste0(prefix, "AB_psi00")]
        expect_equal(cross, unname(OpenMx::mxEvalByName(paste0(prefix, "AB_psi"),
          fit@mxOutput@submodels[[k]])[1, 1]))
        slope <- fit@Estimates$Estimate[fit@Estimates$Name == paste0(prefix, "A_slp_val_est1")]
        expect_equal(slope, unname(OpenMx::mxEvalByName(paste0(prefix, "Aslp_m"),
          fit@mxOutput@submodels[[k]])[1, 1]))
      }
    }
  })
})
