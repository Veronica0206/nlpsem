# Independent construction and algebra regressions for review findings F4-F6.
.review_multi_data <- function(n_traj, records) {
  n <- 12L
  dat <- data.frame(group = rep(1:2, each = n / 2L))
  for (j in sort(unique(unlist(records)))) dat[[paste0("T", j)]] <- rep(j - 1, n)
  for (i in seq_len(n_traj)) for (j in records[[i]]) {
    dat[[paste0(LETTERS[i], j)]] <- 5 * i + j + sin(seq_len(n))
  }
  dat
}

.review_multi_starts <- function(n_traj, curve, identity = FALSE) {
  width <- switch(curve, LIN = 2L, QUAD = 3L, EXP = 3L, JB = 4L)
  mu <- switch(curve, LIN = c(3, 1), QUAD = c(3, 1, .1),
               EXP = c(3, 2, .3), JB = c(3, 1, -1, -.3))
  size <- width * n_traj
  covariance <- if (identity) diag(size) else
    diag(seq(1, 2, length.out = size)) + tcrossprod(seq_len(size) / 100)
  residuals <- diag(n_traj)
  residuals[upper.tri(residuals)] <- .2
  list(means = lapply(seq_len(n_traj), function(i) mu * (1 + (i - 1) / 10)),
       covMatrix = covariance, residuals = residuals)
}

# A mixture needs its full parent data context for algebra evaluation.
.review_multi_contexts <- function(container, n_traj, curve, intrinsic,
                                   records, y_model, starts) {
  dat <- .review_multi_data(n_traj, records)
  args <- list(dat = dat, t_var = rep("T", n_traj), y_var = LETTERS[seq_len(n_traj)],
               curveFun = curve, intrinsic = intrinsic, records = records, y_model = y_model)
  if (container == "single") {
    model <- do.call(nlpsem:::getMGM.mxModel, c(args, list(starts = starts)))
    return(list(list(model = model, parent = model, prefix = "", path = "")))
  }
  if (container == "multigroup") {
    models <- do.call(nlpsem:::getsub.MGM, c(args, list(nClass = 2L,
      starts = list(starts, starts), container = container, grp_var = "group")))
    return(lapply(1:2, function(k) list(model = models[[k]], parent = models[[k]],
      prefix = paste0("c", k), path = "")))
  }
  parent <- do.call(nlpsem:::getMIX.mxModel, c(args, list(nClass = 2L,
    sub_Model = "MGM", cluster_TIC = NULL, m_var = NULL, x_var = NULL,
    x_type = NULL, TVC = NULL, decompose = NULL, growth_TIC = NULL,
    starts = list(starts, starts, c(0, 0)), res_cor = NULL)))
  lapply(1:2, function(k) list(model = parent@submodels[[k]], parent = parent,
    prefix = paste0("c", k), path = paste0("Class", k, ".")))
}

test_that("EXP covariance starts retain every outcome in all model containers", {
  for (n in 2:4) for (y_model in c("LGCM", "LCSM")) {
    for (container in c("single", "multigroup", "mixture")) for (identity in c(TRUE, FALSE)) {
      starts <- .review_multi_starts(n, "EXP", identity)
      contexts <- .review_multi_contexts(container, n, "EXP", FALSE,
                                         rep(list(1:5), n), y_model, starts)
      # Independently select the intercept and amplitude from each 3-factor block.
      keep <- which(rep(c(TRUE, TRUE, FALSE), n))
      expected <- starts$covMatrix[keep, keep, drop = FALSE]
      latents <- paste0(rep(c("eta0", "eta1"), n), rep(LETTERS[seq_len(n)], each = 2))
      for (context in contexts) {
        actual <- unname(context$model$S$values[latents, latents])
        expect_equal(actual, expected)
        expect_gt(min(eigen(actual, symmetric = TRUE, only.values = TRUE)$values), 0)
        expect_true(all(context$model$S$free[latents, latents]))
      }
    }
  }
})

test_that("unequal multivariate wave counts retain independent derived quantities", {
  skip_on_cran()   # builds 18 multivariate configurations across three containers (about 6 s)
  for (curve in c("QUAD", "JB")) {
    for (intrinsic in if (curve == "QUAD") FALSE else c(FALSE, TRUE)) {
      for (container in c("single", "multigroup", "mixture")) for (unequal in c(FALSE, TRUE)) {
        records <- list(1:5, if (unequal) c(1, 2, 4, 5) else 1:5)
        starts <- .review_multi_starts(2L, curve)
        contexts <- .review_multi_contexts(container, 2L, curve, intrinsic,
                                           records, "LCSM", starts)
        for (context in contexts) for (i in 1:2) {
          # Independent derivative of a quadratic or Jenss-Bayley trajectory,
          # evaluated at interval midpoints, plus integration over observed lags.
          times <- records[[i]] - 1
          lags <- diff(times)
          mid <- (head(times, -1) + tail(times, -1)) / 2
          mu <- starts$means[[i]]
          width <- length(mu)
          if (curve == "QUAD") {
            rate <- cbind(1, 2 * mid)
            selected <- 2:3
            mean_factors <- mu[selected]
          } else {
            gamma <- mu[4]
            rate <- cbind(1, gamma * exp(gamma * mid))
            selected <- 2:3
            mean_factors <- mu[selected]
            if (intrinsic) {
              rate <- cbind(rate, mu[3] * exp(gamma * mid) * (1 + gamma * mid))
              selected <- 2:4
              mean_factors <- c(mu[2:3], 0)
            }
          }
          block <- (i - 1L) * width + selected
          psi <- starts$covMatrix[block, block, drop = FALSE]
          slope <- rate %*% mean_factors
          variance <- rate %*% psi %*% t(rate)
          omega <- matrix(0, length(lags), length(lags))
          for (j in seq_along(lags)) omega[j, seq_len(j)] <- lags[seq_len(j)]
          cumulative <- omega %*% rate
          expected <- list(r_loads = rate, slp_m = slope, slp_v = variance,
            chg_inv_m = slope * lags, chg_inv_v = matrix(diag(variance) * lags^2, 1),
            chg_bl_m = cumulative %*% mean_factors,
            chg_bl_v = cumulative %*% psi %*% t(cumulative))
          for (suffix in names(expected)) {
            name <- paste0(context$path, context$prefix, LETTERS[i], suffix)
            actual <- OpenMx::mxEvalByName(name, context$parent, compute = TRUE)
            expect_equal(unname(actual), unname(expected[[suffix]]),
                         info = paste(curve, intrinsic, container, unequal, i, suffix))
          }
        }
      }
    }
  }
})

test_that("multivariate residual covariances match common wave IDs", {
  schedules <- list(identical = list(1:4, 1:4),
                    reordered = list(1:4, 4:1),
                    partial = list(1:4, 3:6))
  for (y_model in c("LGCM", "LCSM")) {
    for (schedule in names(schedules)) for (container in c("single", "multigroup", "mixture")) {
      # LCSM intervals follow the user's chronological record order. Reordering
      # columns is an LGCM invariance check, not a request to reverse intervals.
      if (y_model == "LCSM" && schedule == "reordered") next
      records <- schedules[[schedule]]
      curve <- if (y_model == "LGCM") "LIN" else "QUAD"
      contexts <- .review_multi_contexts(container, 2L, curve, FALSE, records,
                                         y_model, .review_multi_starts(2L, curve))
      shared <- intersect(records[[1]], records[[2]])
      for (context in contexts) {
        model <- context$model
        from <- paste0("A", records[[1]])
        to <- paste0("B", records[[2]])
        expected <- outer(records[[1]], records[[2]], "==")
        expect_equal(unname(model$S$free[from, to, drop = FALSE]), expected)
        expect_equal(unname(model$S$values[from, to, drop = FALSE]), expected * .2)
        actual_labels <- model$S$labels[from, to, drop = FALSE][expected]
        expect_equal(unname(actual_labels), rep(paste0(context$prefix, "AB_RES"), length(shared)))
      }
    }
  }
})

test_that("disjoint multivariate records are rejected clearly before optimization", {
  records <- list(1:3, 4:6)
  starts <- .review_multi_starts(2L, "LIN")
  dat <- .review_multi_data(2L, records)
  pattern <- "Multivariate outcomes 'A' and 'B' must share at least one record ID"
  for (container in c("single", "multigroup", "mixture")) {
    expect_error(.review_multi_contexts(container, 2L, "LIN", FALSE, records,
                                       "LGCM", starts), pattern)
  }
  local_mocked_bindings(mxRun = function(...) stop("Optimization must not be reached"),
                       .package = "nlpsem")
  args <- list(dat = dat, t_var = c("T", "T"), records = records,
               y_var = c("A", "B"), curveFun = "LIN", intrinsic = FALSE, y_model = "LGCM")
  expect_error(do.call(getMGM, c(args, list(starts = starts))), pattern)
  expect_error(do.call(getMGroup, c(args, list(grp_var = "group", sub_Model = "MGM",
                                             starts = list(starts, starts)))), pattern)
  expect_error(do.call(getMIX, c(args, list(prop_starts = c(.5, .5), sub_Model = "MGM",
                                          starts = list(starts, starts, c(0, 0))))), pattern)
})
