.release_med_data <- function(prefixes = c("Y", "M", "X")) {
  dat <- data.frame(g = rep(1:2, each = 8), XB = seq_len(16) / 16)
  for (j in 1:8) {
    dat[[paste0("T", j)]] <- j - 1
    for (prefix in prefixes) dat[[paste0(prefix, j)]] <- j + seq_len(16) / 10
  }
  dat
}

.release_med_starts <- function(curve, xtype) {
  n <- if (curve == "LIN") 2L else 3L
  process <- list(means = if (curve == "LIN") c(0, 1) else c(1, 0, 1, 2.5),
                  covMatrix = diag(n), residuals = 1)
  Y <- M <- X <- process
  M$residuals <- 4
  X$residuals <- 9
  Y$beta_m <- matrix(0, n, n)
  Y$beta_x <- M$beta <- if (xtype == "baseline") rep(0, n) else matrix(0, n, n)
  if (xtype == "baseline") X <- list(means = 0, covMatrix = matrix(1))
  list(Y = Y, M = M, X = X)
}

.release_med_build <- function(container, curve, xtype, records,
                               prefixes = c("Y", "M", "X"), correlations = NULL) {
  nproc <- if (xtype == "baseline") 2L else 3L
  args <- list(dat = .release_med_data(prefixes), t_var = rep("T", nproc),
               y_var = prefixes[1], m_var = prefixes[2],
               x_var = if (xtype == "baseline") "XB" else prefixes[3],
               x_type = xtype, curveFun = curve, records = records,
               starts = .release_med_starts(curve, xtype), res_cor = correlations)
  if (container == "single") return(list(do.call(nlpsem:::getMED.mxModel, args)))
  args$starts <- rep(list(args$starts), 2)
  if (!is.null(correlations)) args$res_cor <- rep(list(correlations), 2)
  do.call(nlpsem:::getsub.MED,
          c(args, list(nClass = 2, container = container, grp_var = "g")))
}

test_that("mediation residuals match wave IDs independently of order and prefixes", {
  schedules <- list(ordered = rep(list(1:4), 3),
                    permuted = list(1:4, 4:1, c(3, 1, 4, 2)),
                    overlap = list(c(1, 3, 5), c(5, 3, 4), c(2, 5, 3)))
  for (container in c("single", "multigroup", "mixture")) {
    for (curve in c("LIN", "BLS")) for (xtype in c("baseline", "longitudinal")) {
      nproc <- if (xtype == "baseline") 2L else 3L
      for (prefixes in list(c("Y", "M", "X"), c("PHQ9_", "M2_", "X7_"))) {
        for (schedule in names(schedules)) {
          records <- schedules[[schedule]][seq_len(nproc)]
          models <- .release_med_build(container, curve, xtype, records, prefixes)
          for (k in seq_along(models)) for (i in seq_len(nproc - 1L)) {
            for (j in seq.int(i + 1L, nproc)) {
              from <- paste0(prefixes[i], records[[i]])
              to <- paste0(prefixes[j], records[[j]])
              expected <- outer(records[[i]], records[[j]], `==`)
              dimnames(expected) <- list(from, to)
              actual <- models[[k]]$S$free[from, to, drop = FALSE]
              expect_identical(actual, expected,
                               info = paste(container, curve, xtype, schedule, i, j, k))
              vals <- models[[k]]$S$values[from, to, drop = FALSE]
              expect_equal(unname(vals[expected]),
                           rep(.3 * sqrt(c(1, 4, 9)[i] * c(1, 4, 9)[j]), sum(expected)))
              expect_true(all(vals[!expected] == 0))
            }
          }
        }
      }
    }
  }
})

test_that("mediation uses each pair's correlation and both process variances", {
  # Distinct correlations prevent a correct variance paired with the wrong
  # correlation from passing. Pair order is YM, YX, MX.
  expected <- c(.1 * sqrt(1 * 4), .2 * sqrt(1 * 9), .3 * sqrt(4 * 9))
  for (container in c("single", "multigroup", "mixture")) for (curve in c("LIN", "BLS")) {
    models <- .release_med_build(container, curve, "longitudinal", rep(list(1:4), 3),
                                correlations = c(.1, .2, .3))
    for (model in models) {
      expect_equal(c(model$S$values["Y1", "M1"], model$S$values["Y1", "X1"],
                     model$S$values["M1", "X1"]), expected)
    }
  }
})

test_that("disjoint mediation schedules fail with an explicit model requirement", {
  for (container in c("single", "multigroup", "mixture")) for (curve in c("LIN", "BLS")) {
    expect_error(.release_med_build(container, curve, "baseline", list(1:3, 4:6)),
                 "must share at least one record ID")
    expect_error(.release_med_build(container, curve, "longitudinal", list(1:3, 1:3, 4:6)),
                 "must share at least one record ID")
  }
})
