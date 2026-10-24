# Distinct process times and zero covariance starts expose specification errors
# that same-time examples and nonzero starting values cannot detect.
review_spec_data <- function() {
  dat <- data.frame(g = rep(1:2, each = 6), Z1 = seq_len(12) / 10,
                    Z2 = (seq_len(12) %% 5) / 10, XB = seq_len(12) / 12)
  for (j in 1:4) {
    dat[[paste0("TY", j)]] <- j
    dat[[paste0("TM", j)]] <- j * 10
    dat[[paste0("TX", j)]] <- j * 100
    dat[[paste0("T", j)]] <- j - 1
    for (role in c("Y", "M", "X")) dat[[paste0(role, j)]] <- j + seq_len(12) / 10
  }
  dat
}

review_med_starts <- function(curve, xtype) {
  n <- if (curve == "LIN") 2L else 3L
  process <- list(means = if (curve == "LIN") c(0, 1) else c(1, 0, 1, 2.5),
                  covMatrix = diag(n), residuals = 1)
  Y <- M <- X <- process
  if (curve == "BLS") {
    M$means[4] <- 25
    X$means[4] <- 250
  }
  Y$beta_m <- matrix(0, n, n)
  Y$beta_x <- M$beta <- if (xtype == "baseline") rep(0, n) else matrix(0, n, n)
  if (xtype == "baseline") X <- list(means = 0, covMatrix = matrix(1))
  list(Y = Y, M = M, X = X)
}

review_med_models <- function(container, curve, xtype, dat, records) {
  # readr lazily caches its locale when the builder parses wave labels.
  # Keep that incidental cache from leaking into other test files.
  old_locale <- getOption("readr.default_locale")
  on.exit(options(readr.default_locale = old_locale), add = TRUE)
  nproc <- if (xtype == "baseline") 2L else 3L
  args <- list(dat = dat, t_var = c("TY", "TM", "TX")[seq_len(nproc)],
               y_var = "Y", m_var = "M", x_var = if (xtype == "baseline") "XB" else "X",
               x_type = xtype, curveFun = curve, records = records,
               starts = review_med_starts(curve, xtype), res_cor = NULL)
  if (container == "single") return(list(do.call(nlpsem:::getMED.mxModel, args)))
  args$starts <- rep(list(args$starts), 2)
  do.call(nlpsem:::getsub.MED,
          c(args, list(nClass = 2, container = container, grp_var = "g")))
}

test_that("mediation loadings use the time scale of each process in every container", {
  dat <- review_spec_data()
  scales <- c(Y = 1, M = 10, X = 100)
  for (container in c("single", "multigroup", "mixture")) {
    for (xtype in c("baseline", "longitudinal")) {
      roles <- if (xtype == "baseline") c("Y", "M") else c("Y", "M", "X")
      for (curve in c("LIN", "BLS")) {
        models <- review_med_models(container, curve, xtype, dat,
                                     rep(list(1:4), length(roles)))
        for (k in seq_along(models)) {
          model <- models[[k]]
          # Mixture children normally inherit the parent's data. Supply that
          # same data to this inspection copy so definition variables evaluate.
          if (container == "mixture") model <- OpenMx::mxModel(model, OpenMx::mxData(dat, type = "raw"))
          pfx <- if (container == "single") "" else paste0("c", k)
          context <- paste(container, xtype, curve, k)
          for (role in roles) for (wave in c(2, 4)) {
            expected_time <- unname(scales[role]) * wave
            expected_knot <- unname(scales[role]) * 2.5
            expect_identical(model@matrices[[paste0("t", role, wave)]]@labels[1, 1],
                             paste0("data.T", role, wave), info = context)
            for (loading in if (curve == "LIN") 1L else 1:2) {
              expected <- if (curve == "LIN") expected_time else {
                if (loading == 1L) min(0, expected_time - expected_knot) else max(0, expected_time - expected_knot)
              }
              actual <- OpenMx::mxEvalByName(paste0(pfx, "L", loading, wave, role),
                                             model, compute = TRUE)
              expect_equal(as.numeric(actual), expected, info = context)
            }
          }
        }
      }
    }
  }
})

test_that("unequal mediation schedules only reference their own selected time columns", {
  dat <- review_spec_data()
  dat$TM4 <- dat$M4 <- NULL
  dat$TX3 <- dat$X3 <- NULL
  for (container in c("single", "multigroup", "mixture")) {
    for (xtype in c("baseline", "longitudinal")) for (curve in c("LIN", "BLS")) {
      records <- if (xtype == "baseline") list(1:4, 1:3) else list(1:4, 1:3, c(1, 2, 4))
      models <- review_med_models(container, curve, xtype, dat, records)
      for (k in seq_along(models)) {
        model <- models[[k]]
        if (container == "mixture") model <- OpenMx::mxModel(model, OpenMx::mxData(dat, type = "raw"))
        pfx <- if (container == "single") "" else paste0("c", k)
        expected <- if (curve == "LIN") 4 else 1.5
        algebra <- paste0(pfx, if (curve == "LIN") "L14Y" else "L24Y")
        expect_equal(as.numeric(OpenMx::mxEvalByName(algebra, model, compute = TRUE)), expected)
        expect_false("tM4" %in% names(model@matrices))
        if (xtype == "longitudinal") expect_false("tX3" %in% names(model@matrices))
      }
    }
  }
})

review_tvc_starts <- function(eps, nTIC = 1L, decompose = 1L) {
  BL <- diag(nTIC + 1L)
  BL[row(BL) != col(BL)] <- eps
  list(
    Y_starts = list(means = c(0, 1), covMatrix = c(1, 0, 1), residuals = 1,
                    rel_rate = c(1, 1, 1)),
    TVC_starts = if (decompose) {
      list(means = c(0, 1), covMatrix = c(1, eps, 1), residuals = 1, rel_rate = c(1, 1, 1))
    } else list(means = rep(0, 4), covMatrix = rep(1, 4)),
    BL_starts = list(
      means = rep(0, if (decompose) nTIC + 1L else nTIC),
      covMatrix = if (decompose) BL else BL[seq_len(nTIC), seq_len(nTIC), drop = FALSE][lower.tri(diag(nTIC), diag = TRUE)]
    ),
    beta = matrix(0, 2, if (decompose) nTIC + 1L else nTIC), kappa = .1, XY_cov = .1
  )
}

review_tvc_models <- function(eps, container, decompose, y_model, nTIC = 1L) {
  args <- list(dat = review_spec_data(), t_var = "T", y_var = "Y",
               curveFun = if (y_model == "LCSM") "NonP" else "LIN", intrinsic = FALSE,
               records = 1:4, y_model = y_model, TVC = "X", decompose = decompose,
               growth_TIC = if (nTIC) paste0("Z", seq_len(nTIC)) else NULL,
               starts = review_tvc_starts(eps, nTIC, decompose))
  if (container == "single") return(list(do.call(nlpsem:::getTVC.mxModel, args)))
  args$starts <- rep(list(args$starts), 2)
  do.call(nlpsem:::getsub.TVC,
          c(args, list(nClass = 2, container = container, grp_var = "g")))
}

test_that("zero TVC covariance starts preserve free parameters and structural constraints", {
  for (container in c("single", "multigroup", "mixture")) for (decompose in 1:3) {
    for (y_model in c("LGCM", "LCSM")) for (nTIC in 1:2) {
      zero <- review_tvc_models(0, container, decompose, y_model, nTIC)
      small <- review_tvc_models(.1, container, decompose, y_model, nTIC)
      for (k in seq_along(zero)) {
        pfx <- if (container == "single") "" else paste0("c", k)
        zero_free <- names(OpenMx::omxGetParameters(zero[[k]], free = TRUE))
        small_free <- names(OpenMx::omxGetParameters(small[[k]], free = TRUE))
        context <- paste(container, decompose, y_model, nTIC, k)
        expect_identical(sort(zero_free), sort(small_free), info = context)
        expected_free <- paste0(pfx, c("X_psi01", paste0("covBL", seq_len(nTIC)),
                                       if (nTIC == 2) "phi12"))
        expect_true(all(expected_free %in% zero_free), info = context)
        expect_identical(zero[[k]]$S@free, small[[k]]$S@free, info = context)
        # The expected mask is based on the model: TICs covary with each other
        # and eta0x; only their covariances with eta1x are structural zeros.
        vars <- c(paste0("Z", seq_len(nTIC)), "eta0x", "eta1x")
        expected_mask <- matrix(TRUE, length(vars), length(vars), dimnames = list(vars, vars))
        expected_mask[seq_len(nTIC), "eta1x"] <- FALSE
        expected_mask["eta1x", seq_len(nTIC)] <- FALSE
        expect_identical(zero[[k]]$S@free[vars, vars], expected_mask, info = context)
        expect_equal(unname(zero[[k]]$S@values[paste0("Z", seq_len(nTIC)), "eta1x"]),
                     rep(0, nTIC), info = context)
        expect_equal(unname(small[[k]]$S@values[paste0("Z", seq_len(nTIC)), "eta1x"]),
                     rep(0, nTIC), info = context)
      }
    }
  }
})

test_that("TVC zero-start invariance also holds without TICs and without decomposition", {
  cases <- rbind(data.frame(decompose = 1:3, nTIC = 0L),
                 data.frame(decompose = 0L, nTIC = 1:2))
  for (i in seq_len(nrow(cases))) {
    zero <- review_tvc_models(0, "single", cases$decompose[i], "LGCM", cases$nTIC[i])[[1]]
    small <- review_tvc_models(.1, "single", cases$decompose[i], "LGCM", cases$nTIC[i])[[1]]
    expect_setequal(names(OpenMx::omxGetParameters(zero, free = TRUE)),
                    names(OpenMx::omxGetParameters(small, free = TRUE)))
    expect_identical(zero$S@free, small$S@free)
  }
})
