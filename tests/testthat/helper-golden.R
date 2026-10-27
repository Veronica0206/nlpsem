# Shared definitions for the fitting-layer reference ("golden") tests. The same case list drives
# data-raw/make_golden_fits.R, which writes tests/testthat/golden/fit_golden.rds.

# A 200-person subsample of RMS_dat (or the first n persons) with re-baselined occasions, one
# standardised growth TIC (ex1), one standardised cluster TIC (gx1) and a standardised time-varying
# covariate (Rs).
.golden_data <- function(n = 200) {
  env <- new.env()
  data("RMS_dat", package = "nlpsem", envir = env)
  dat <- env$RMS_dat[seq_len(n), ]
  baseT <- dat$T1
  for (j in 1:9) dat[[paste0("T", j)]] <- (dat[[paste0("T", j)]] - baseT) / 12
  dat$ex1 <- as.numeric(scale(dat$Approach_to_Learning))
  dat$gx1 <- as.numeric(scale(dat$INCOME))
  BL_mean <- mean(dat$R1); BL_sd <- sd(dat$R1)
  for (j in 1:9) dat[[paste0("Rs", j)]] <- (dat[[paste0("R", j)]] - BL_mean) / BL_sd
  dat
}

.golden_cases <- function() {
  uni <- function(fn, curveFun, intrinsic, growth_TIC = NULL)
    list(fn = fn, args = list(t_var = "T", y_var = "M", curveFun = curveFun, intrinsic = intrinsic,
                              records = 1:9, growth_TIC = growth_TIC))
  tvc <- function(curveFun, intrinsic, y_model, decompose, growth_TIC = NULL)
    list(fn = "getTVCmodel", args = list(t_var = "T", y_var = "M", curveFun = curveFun, intrinsic = intrinsic,
                                         records = 1:9, y_model = y_model, TVC = "Rs", decompose = decompose,
                                         growth_TIC = growth_TIC))
  mgm <- function(curveFun, intrinsic, y_model)
    list(fn = "getMGM", args = list(t_var = c("T", "T"), y_var = c("R", "M"), curveFun = curveFun,
                                    intrinsic = intrinsic, records = list(1:9, 1:9), y_model = y_model))
  med <- function(curveFun, x_type)
    list(fn = "getMediation",
         args = if (x_type == "baseline")
           list(t_var = rep("T", 2), y_var = "M", m_var = "R", x_type = "baseline", x_var = "ex1",
                curveFun = curveFun, records = list(1:9, 1:9))
         else
           list(t_var = rep("T", 3), y_var = "S", m_var = "M", x_type = "longitudinal", x_var = "R",
                curveFun = curveFun, records = list(2:9, 1:9, 1:9)))
  grp <- function(sub_Model, extra)
    list(fn = "getMGroup", args = c(list(grp_var = "SEX", sub_Model = sub_Model), extra))
  mix <- function(sub_Model, extra, cluster_TIC = NULL)
    list(fn = "getMIX", args = c(list(prop_starts = c(0.5, 0.5), sub_Model = sub_Model, cluster_TIC = cluster_TIC), extra))
  lgcm_args <- function(curveFun, intrinsic) list(t_var = "T", y_var = "M", curveFun = curveFun, intrinsic = intrinsic, records = 1:9)
  cases <- list(
    LGCM_LIN = uni("getLGCM", "LIN", FALSE), LGCM_QUAD = uni("getLGCM", "QUAD", FALSE),
    LGCM_EXP_f = uni("getLGCM", "EXP", FALSE), LGCM_EXP_r = uni("getLGCM", "EXP", TRUE),
    LGCM_JB_f = uni("getLGCM", "JB", FALSE), LGCM_JB_r = uni("getLGCM", "JB", TRUE),
    LGCM_BLS_f = uni("getLGCM", "BLS", FALSE), LGCM_BLS_r = uni("getLGCM", "BLS", TRUE),
    LGCM_BLS_f_TIC = uni("getLGCM", "BLS", FALSE, "ex1"),
    LCSM_NonP = uni("getLCSM", "NonP", FALSE), LCSM_QUAD = uni("getLCSM", "QUAD", FALSE),
    LCSM_EXP_f = uni("getLCSM", "EXP", FALSE), LCSM_EXP_r = uni("getLCSM", "EXP", TRUE),
    LCSM_JB_f = uni("getLCSM", "JB", FALSE), LCSM_JB_r = uni("getLCSM", "JB", TRUE),
    LCSM_NonP_TIC = uni("getLCSM", "NonP", FALSE, "ex1"),
    TVC_LGCM_LIN_d0 = tvc("LIN", FALSE, "LGCM", 0),
    TVC_LGCM_BLS_f_d1_TIC = tvc("BLS", FALSE, "LGCM", 1, "ex1"),
    TVC_LGCM_BLS_r_d1 = tvc("BLS", TRUE, "LGCM", 1),
    TVC_LGCM_EXP_r_d2 = tvc("EXP", TRUE, "LGCM", 2),
    TVC_LCSM_NonP_d3 = tvc("NonP", FALSE, "LCSM", 3),
    MGM_LGCM_LIN = mgm("LIN", FALSE, "LGCM"), MGM_LGCM_BLS_f = mgm("BLS", FALSE, "LGCM"),
    MGM_LGCM_BLS_r = mgm("BLS", TRUE, "LGCM"), MGM_LCSM_NonP = mgm("NonP", FALSE, "LCSM"),
    MED_LIN_baseline = med("LIN", "baseline"), MED_BLS_baseline = med("BLS", "baseline"),
    MED_LIN_longitudinal = med("LIN", "longitudinal"), MED_BLS_longitudinal = med("BLS", "longitudinal"),
    MGroup_LGCM_BLS_f = grp("LGCM", lgcm_args("BLS", FALSE)),
    MGroup_LCSM_NonP = grp("LCSM", lgcm_args("NonP", FALSE)),
    MGroup_TVC_LIN_d0 = grp("TVC", list(t_var = "T", y_var = "M", curveFun = "LIN", intrinsic = FALSE, records = 1:9,
                                        y_model = "LGCM", TVC = "Rs", decompose = 0)),
    MGroup_MGM_LIN = grp("MGM", list(t_var = c("T", "T"), y_var = c("R", "M"), curveFun = "LIN", intrinsic = FALSE,
                                     records = list(1:9, 1:9), y_model = "LGCM")),
    MGroup_MED_LIN = grp("MED", list(t_var = rep("T", 2), y_var = "M", m_var = "R", x_type = "baseline", x_var = "ex1",
                                     curveFun = "LIN", records = list(1:9, 1:9))),
    MIX_LGCM_LIN = mix("LGCM", lgcm_args("LIN", FALSE)),
    MIX_LGCM_BLS_f = mix("LGCM", lgcm_args("BLS", FALSE)),
    MIX_LGCM_BLS_f_clusterTIC = mix("LGCM", lgcm_args("BLS", FALSE), cluster_TIC = "gx1"),
    MIX_LCSM_NonP = mix("LCSM", lgcm_args("NonP", FALSE)),
    MIX_TVC_LIN_d0 = mix("TVC", list(t_var = "T", y_var = "M", curveFun = "LIN", intrinsic = FALSE, records = 1:9,
                                     y_model = "LGCM", TVC = "Rs", decompose = 0)),
    MIX_MGM_LIN = mix("MGM", list(t_var = c("T", "T"), y_var = c("R", "M"), curveFun = "LIN", intrinsic = FALSE,
                                  records = list(1:9, 1:9), y_model = "LGCM")),
    MIX_MED_LIN = mix("MED", list(t_var = rep("T", 2), y_var = "M", m_var = "R", x_type = "baseline", x_var = "ex1",
                                  curveFun = "LIN", records = list(1:9, 1:9))),
    # Jenss-Bayley models on the reading scores of all 500 children, which accelerate clearly over the
    # study (the math-score JB cases above do not converge)
    LGCM_JB_f_read = c(list(n = 500L), uni("getLGCM", "JB", FALSE)),
    LCSM_JB_f_read = c(list(n = 500L), uni("getLCSM", "JB", FALSE))
  )
  for (nm in c("LGCM_JB_f_read", "LCSM_JB_f_read")) cases[[nm]]$args$y_var <- "R"
  # Cases that do not reach status 0 from the default starts on the subsample get a few jittered
  # extra attempts (mxTryHard); the seed set in .fit_golden_case() keeps them reproducible.
  light <- c("LGCM_EXP_f", "LGCM_JB_f", "LGCM_JB_r", "LCSM_EXP_f", "LCSM_JB_f", "LCSM_JB_r", "MED_LIN_baseline")
  heavy <- c("TVC_LGCM_BLS_f_d1_TIC", "TVC_LGCM_BLS_r_d1", "MGM_LGCM_BLS_f", "MGM_LGCM_BLS_r", "MGM_LCSM_NonP",
             "MIX_LGCM_BLS_f", "MIX_LGCM_BLS_f_clusterTIC", "MIX_LCSM_NonP", "MIX_MGM_LIN")
  for (nm in light) cases[[nm]]$tries <- 5L
  for (nm in heavy) cases[[nm]]$tries <- 3L
  for (nm in c("LGCM_JB_f_read", "LCSM_JB_f_read")) cases[[nm]]$tries <- 10L
  cases
}

# The data of a case: the 200-person subsample unless the case names a larger sample.
.golden_case_data <- function(case, dat = .golden_data()) {
  if (is.null(case$n) || case$n == nrow(dat)) dat else .golden_data(case$n)
}

# Fit one case deterministically (fixed optimizer, fixed seed for the k-means step of mixtures) and
# return the quantities that are compared: -2 log-likelihood, free-parameter estimates, status code,
# and the problems getDiagnostics() reports (covariance admissibility, Hessian, standard errors).
.fit_golden_case <- function(case, dat) {
  OpenMx::mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)
  set.seed(20191029)
  fn <- get(case$fn, envir = asNamespace("nlpsem"))
  args <- case$args
  if (!is.null(case$tries)) args <- c(args, list(tries = case$tries))
  fit <- suppressMessages(suppressWarnings(do.call(fn, c(list(dat = dat), args))))
  model <- fit@mxOutput
  list(minus2LL = unname(model@output$Minus2LogLikelihood),
       estimates = OpenMx::omxGetParameters(model),
       status = unname(model@output$status$code),
       issues = getDiagnostics(fit)$issues,
       seconds = unname(model@output$wallTime))
}
