# Numerical checks of the inferential output, rather than just ggplot class checks.

.plot_test_data <- function(fit, prefix = "Y", curve = "BLS", family = "LGCM", y_model = "LGCM", times = 20,
                            group = NULL){
  eval_model <- if (is.null(group)) fit else fit@submodels[[group]]
  eval_prefix <- if (is.null(group)) prefix else paste0("c", group, prefix)
  se_prefix <- if (is.null(group)) prefix else paste0("Class", group, ".c", group, prefix)
  nlpsem:::.trajectory_plot_data(fit, eval_model, eval_prefix, se_prefix, family, curve, y_model, times)
}

.expect_plot_se <- function(plot_data, column, expression, fit, tolerance = 3e-4){
  # mxSE differentiates the complete nonlinear expression numerically; allow its
  # finite-difference error. Analytic gradients are independently checked below.
  expected <- c(OpenMx::mxSE(expression, fit, silent = TRUE))
  actual <- (plot_data[[paste0(column, "_u")]] - plot_data[[paste0(column, "_l")]]) / (2 * qnorm(.975))
  expect_equal(actual, expected, tolerance = tolerance)
}

test_that("pointwise intervals retain non-diagonal parameter covariance and untruncated bounds", {
  dist <- list(est = c(10, 2), cov = matrix(c(4, -1, -1, 1), 2))
  out <- nlpsem:::.plot_curve_delta(dist, times = 2, curveFun = "LIN")
  # Var(a + 2b) = 4 + 4*1 + 4*(-1) = 4.
  expect_equal(out$est, 14)
  expect_equal(out$lower, 14 - qnorm(.975) * 2)
  expect_equal(out$upper, 14 + qnorm(.975) * 2)
  wide <- nlpsem:::.plot_curve_delta(list(est = c(10, 0), cov = diag(c(10000, 1))), 0, "LIN")
  expect_equal(wide$lower, 10 - qnorm(.975) * 100)
  expect_equal(wide$upper, 10 + qnorm(.975) * 100)
})

test_that("curve gradients agree with an independent finite-difference calculation", {
  cases <- list(
    LIN = list(p = c(10, 2), f = function(p,t) p[1]+p[2]*t, r = function(p,t) rep(p[2],length(t))),
    QUAD = list(p = c(10, 2, -.05), f = function(p,t) p[1]+p[2]*t+p[3]*t^2, r = function(p,t) p[2]+2*p[3]*t),
    EXP = list(p = c(10, 12, .15), f = function(p,t) p[1]+p[2]*(1-exp(-p[3]*t)), r = function(p,t) p[2]*p[3]*exp(-p[3]*t)),
    JB = list(p = c(10, 2, -6, -.15), f = function(p,t) p[1]+p[2]*t+p[3]*(exp(p[4]*t)-1), r = function(p,t) p[2]+p[3]*p[4]*exp(p[4]*t)),
    BLS = list(p = c(10, 2, .5, 3), f = function(p,t) ifelse(t<=p[4],p[1]+p[2]*t,p[1]+p[2]*p[4]+p[3]*(t-p[4])))
  )
  tt <- c(1, 5)
  for (curve in names(cases)){
    spec <- cases[[curve]]; p <- spec$p
    # A positive definite matrix with nonzero off-diagonals of both signs.
    a <- diag(length(p)); a[1, 2] <- -.4
    covariance <- crossprod(a) * .01
    for (quantity in if (curve == "BLS") "status" else c("status", "rate", "change")){
      fn <- if (quantity == "rate") spec$r else spec$f
      if (quantity == "change") fn <- local({ff <- spec$f; function(p,t) ff(p,t)-ff(p,0)})
      jac <- sapply(seq_along(p), function(j){
        plus <- minus <- p; plus[j] <- plus[j]+1e-5; minus[j] <- minus[j]-1e-5
        (fn(plus,tt)-fn(minus,tt))/(2e-5)
      })
      expected_se <- sqrt(diag(jac %*% covariance %*% t(jac)))
      out <- nlpsem:::.plot_curve_delta(list(est=p,cov=covariance),tt,curve,quantity)
      expect_equal(out$est,fn(p,tt),info=paste(curve,quantity))
      expect_equal((out$upper-out$lower)/(2*qnorm(.975)),expected_se,tolerance=1e-6,info=paste(curve,quantity))
    }
  }
})

test_that("the reported BLS example has covariance-aware public plot intervals", {
  e <- .load_extdata("getLGCM_examples.RData")
  fit <- e$Math_LGCM_BLS_f@mxOutput
  out <- getFigure(fit, sub_Model="LGCM", y_var="M", curveFun="BLS", t_var="T", records=1:9, xstarts=0)
  dat <- out@figures[[1]][[1]]$data
  at20 <- dat[which.min(abs(dat$time-20)),]
  expect_equal(at20$time,20)
  mean0 <- c(OpenMx::mxEvalByName("Y_mean0", fit))
  expect_equal(at20$Y_hat, mean0[1] + 20 * mean0[2], tolerance=1e-10)
  .expect_plot_se(at20,"Y_hat","Y_mean0[1,1] + 20 * Y_mean0[2,1]",fit)
  # After the knot, uncertainty in the estimated knot contributes to the interval.
  after <- .plot_test_data(fit,times=40)
  .expect_plot_se(after,"Y_hat",paste0("Y_mean0[1,1] + Y_mean0[2,1]*Y_mean0[4,1] + ",
                                      "Y_mean0[3,1]*(40-Y_mean0[4,1])"),fit)
  expect_match(as.character(out@figures[[1]][[1]]$layers[[2]]$mapping$linetype),"Pointwise")
})

test_that("fixed knots and covariate-adjusted mean uncertainty are propagated", {
  e <- .load_extdata("getLGCM_examples.RData")
  for (nm in c("Math_LGCM_BLS_r","Math_LGCM_TIC_BLS_f")){
    fit <- e[[nm]]@mxOutput
    out <- .plot_test_data(fit,times=40)
    knot <- if(length(OpenMx::mxEvalByName("Y_mean0",fit))==4) "Y_mean0[4,1]" else "Y_alpha0[4,1]"
    expression <- paste0("Y_mean0[1,1]+Y_mean0[2,1]*",knot,"+Y_mean0[3,1]*(40-",knot,")")
    .expect_plot_se(out,"Y_hat",expression,fit)
  }
})

test_that("parametric change-score intervals include random-rate uncertainty", {
  e <- .load_extdata("getLCSM_examples.RData")
  cases <- list(
    list(nm="Read_LCSM_QUAD",curve="QUAD"),
    list(nm="Read_LCSM_EXP_r",curve="EXP"),
    list(nm="Read_LCSM_JB_r",curve="JB")
  )
  for (case in cases){
    fit <- e[[case$nm]]@mxOutput
    times <- c(OpenMx::mxEvalByName("M_TIME",fit))[-1]
    out <- .plot_test_data(fit,curve=case$curve,family="LCSM",y_model="LCSM",times=times)
    .expect_plot_se(out,"cY_hat","Ychg_bl_m",fit)
    .expect_plot_se(out,"dY_hat","Yslp_m",fit)
  }
})

test_that("nonparametric rate intervals include uncertainty in estimated relative rates", {
  e <- .load_extdata("getLCSM_examples.RData")
  for (nm in c("Read_LCSM_NonP","Read_LCSM_NonP_TIC")){
    fit <- e[[nm]]@mxOutput
    times <- c(OpenMx::mxEvalByName("M_TIME",fit))[-1]
    out <- .plot_test_data(fit,curve="NonP",family="LCSM",y_model="LCSM",times=times)
    .expect_plot_se(out,"dY_hat","Yslp_m",fit)
    .expect_plot_se(out,"cY_hat","Ychg_bl_m",fit)
  }
})

test_that("TVC, multivariate and mediation curves propagate their joint covariance", {
  tvc <- .load_extdata("getTVCmodel_examples.RData")$Math_TVCslp_BLS_f@mxOutput
  knot <- c(OpenMx::mxEvalByName("Y_mean0",tvc))[4]
  tt <- knot + 1
  out <- .plot_test_data(tvc,family="TVC",times=tt)
  .expect_plot_se(out,"Y_hat",paste0("Y_mean0[1,1]+Y_mean0[2,1]*Y_mean0[4,1]+Y_mean0[3,1]*(",tt,"-Y_mean0[4,1])"),tvc)
  mgm <- .load_extdata("getMGM_examples.RData")$RM_PLGCM.f@mxOutput
  for(prefix in c("R","M")){
    out <- .plot_test_data(mgm,prefix=prefix,family="MGM",times=20)
    .expect_plot_se(out,"Y_hat",paste0(prefix,"_mean0[1,1]+20*",prefix,"_mean0[2,1]"),mgm)
  }
  med <- .load_extdata("getMediation_examples.RData")$Med2_LGCM_BLS@mxOutput
  for(prefix in c("Y","M")){
    out <- .plot_test_data(med,prefix=prefix,family="MED",times=20)
    .expect_plot_se(out,"Y_hat",paste0(prefix,"_mean[2,1]+",prefix,"_mean[1,1]*(20-",prefix,"_mug)"),med)
  }
})

test_that("multiple-group and mixture curves use the parent fitted covariance", {
  cases <- list(
    .load_extdata("getMGroup_examples.RData")$MGroup_Math_BLS_LGCM_f@mxOutput,
    .load_extdata("getMIX_examples.RData")$Math_BLS_LGCM3@mxOutput
  )
  for(fit in cases){
    for(k in seq_along(fit@submodels)){
      prefix <- paste0("Class",k,".c",k,"Y")
      mean0 <- c(OpenMx::mxEvalByName(paste0("c",k,"Y_mean0"),fit@submodels[[k]]))
      knot_expr <- if(length(mean0)==4) paste0(prefix,"_mean0[4,1]") else paste0(prefix,"_alpha0[4,1]")
      knot <- if(length(mean0)==4) mean0[4] else c(OpenMx::mxEvalByName(paste0("c",k,"Y_alpha0"),fit@submodels[[k]]))[4]
      tt <- knot+1
      out <- .plot_test_data(fit,group=k,times=tt)
      expression <- paste0(prefix,"_mean0[1,1]+",prefix,"_mean0[2,1]*",knot_expr,"+",prefix,"_mean0[3,1]*(",tt,"-",knot_expr,")")
      .expect_plot_se(out,"Y_hat",expression,fit)
    }
  }
})

test_that("an uncertain spline knot omits only the undefined pointwise interval", {
  dist <- list(est=c(10,2,.5,3),cov=diag(4))
  out <- nlpsem:::.plot_curve_delta(dist,c(2,3,4),"BLS")
  expect_true(is.na(out$lower[2]) && is.na(out$upper[2]))
  expect_true(all(is.finite(out$lower[c(1,3)])))
  expect_true(all(is.finite(out$est)))
  dist$cov[4,4] <- 0
  expect_true(all(is.finite(nlpsem:::.plot_curve_delta(dist,c(2,3,4),"BLS")$lower)))
})

test_that("multi-character prefixes preserve observed growth and change layers", {
  cases <- list(
    list(fit=.load_extdata("getLGCM_examples.RData")$Math_LGCM_BLS_f@mxOutput,family="LGCM",curve="BLS",old="M"),
    list(fit=.load_extdata("getLCSM_examples.RData")$Read_LCSM_NonP@mxOutput,family="LCSM",curve="NonP",old="R")
  )
  for(case in cases){
    fit <- case$fit
    for(j in 1:9){
      fit$data$observed[[paste0("Outcome",j)]] <- fit$data$observed[[paste0(case$old,j)]]
      fit$data$observed[[paste0("Time",j)]] <- fit$data$observed[[paste0("T",j)]]
    }
    fig <- suppressWarnings(getFigure(fit,sub_Model=case$family,y_var="Outcome",curveFun=case$curve,
                                       t_var="Time",records=1:9,xstarts=0))
    observed <- fig@figures[[1]][[1]]$layers[[4]]$data
    expected_waves <- if(case$family=="LGCM") 9 else 8
    expect_equal(nrow(observed),nrow(fit$data$observed)*expected_waves)
    expect_equal(length(unique(observed$Wave)),expected_waves)
  }
})


test_that("parametric change intervals use the first selected time as baseline", {
  dist <- list(est = c(10, 2, 1), cov = diag(3))
  out <- nlpsem:::.plot_curve_delta(dist, times = 3, curveFun = "QUAD",
                                   quantity = "change", baseline_time = 2)
  # (10 + 2*3 + 3^2) - (10 + 2*2 + 2^2) = 7; gradient = (0, 1, 5).
  expect_equal(out$est, 7)
  expect_equal(out$upper, 7 + qnorm(.975) * sqrt(26))
  expect_equal(out$lower, 7 - qnorm(.975) * sqrt(26))
})

test_that("the public change plot uses a nonzero first selected observation", {
  fit <- .load_extdata("getLCSM_examples.RData")$Read_LCSM_QUAD@mxOutput
  figure <- getFigure(fit,sub_Model="LCSM",y_var="R",curveFun="QUAD",t_var="T",records=2:9,xstarts=0)
  out <- figure@figures[[1]][[1]]$data
  t0 <- mean(fit$data$observed$T2)
  tt <- mean(fit$data$observed$T3)
  gf <- c(OpenMx::mxEvalByName("Y_mean0",fit))
  expect_equal(out$cY_hat[1],gf[2]*(tt-t0)+gf[3]*(tt^2-t0^2))
  expression <- paste0("Y_mean0[2,1]*",tt-t0,"+Y_mean0[3,1]*",tt^2-t0^2)
  .expect_plot_se(out[1,],"cY_hat",expression,fit)
})

test_that("class-specific nonlinear TIC trajectories use adjusted mean and rate covariance", {
  fit <- .load_extdata("getMGroup_examples.RData")$MGroup_Read_EXP_LGCM_TIC_r@mxOutput
  for(k in seq_along(fit@submodels)){
    out <- .plot_test_data(fit,curve="EXP",group=k,times=2)
    prefix <- paste0("Class",k,".c",k,"Y")
    gf <- c(OpenMx::mxEvalByName(paste0("c",k,"Y_mean0"),fit@submodels[[k]]))
    rate <- if(length(gf)==3) paste0(prefix,"_mean0[3,1]") else paste0(prefix,"_alpha0[3,1]")
    expression <- paste0(prefix,"_mean0[1,1]+",prefix,"_mean0[2,1]*(1-exp(-",rate,"*2))")
    .expect_plot_se(out,"Y_hat",expression,fit)
  }
})

test_that("TVC plots label the growth component and add the full model-implied mean", {
  fit <- .load_extdata("getTVCmodel_examples.RData")$Math_TVC_BLS_f@mxOutput
  fig <- suppressWarnings(getFigure(fit, sub_Model = "TVC", y_var = "M", curveFun = "BLS", y_model = "LGCM", t_var = "T",
                                    records = 1:9, xstarts = 0))
  p <- fig@figures[[1]][[1]]
  labels <- ggplot2::ggplot_build(p)$plot$scales$get_scales("linetype")$get_labels()
  expect_true("Model Implied Growth Component (TVC Effect Excluded)" %in% labels)
  expect_true("Model Implied Mean Outcome (TVC Effect Included)" %in% labels)
  expect_false("Model Implied Growth Status" %in% labels)
  # With a direct TVC, the full mean is the growth component plus kappa times the TVC mean. The model's
  # loadings smooth |t - knot| as sqrt((t - knot)^2 + 1e-6) and the plotted curve does not, so at an occasion
  # near the knot the two differ by up to the half-difference of the slopes times 1e-3.
  full <- p$layers[[5]]$data
  times <- colMeans(fit$data$observed[, paste0("T", 1:9)])
  component <- nlpsem:::.plot_curve_delta(nlpsem:::.plot_gf_distribution(fit, fit, "Y", "Y", "TVC", "BLS"), times, "BLS")$est
  par <- OpenMx::omxGetParameters(fit)
  expect_equal(full$time, unname(times))
  gap <- abs(full$value - component - unname(par["kappa"] * par[paste0("TVC_m", 1:9)]))
  knot <- c(OpenMx::mxEvalByName("Y_mean0", fit))[4]
  expect_lt(max(gap[abs(times - knot) > 0.1]), 1e-5)
  expect_lt(max(gap), 1e-3)
})

test_that("a fit with the records in descending order can be plotted", {
  skip_on_cran()
  dat <- .golden_data(n = 200)
  OpenMx::mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)
  fit_records <- function(records){
    set.seed(20191029)
    suppressMessages(getLGCM(dat, t_var = "T", y_var = "M", curveFun = "LIN", intrinsic = FALSE, records = records))
  }
  descending <- plot(fit_records(9:1), draw = FALSE)@figures[[1]][[1]]$data
  ascending <- plot(fit_records(1:9), draw = FALSE)@figures[[1]][[1]]$data
  expect_equal(descending$time, ascending$time)
  expect_equal(descending$Y_hat, ascending$Y_hat, tolerance = 1e-6)
})

test_that("mixture and multiple-group TVC plots add each class's full model-implied mean", {
  skip_on_cran()
  dat <- .golden_data(n = 200)
  OpenMx::mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)
  args <- list(t_var = "T", y_var = "M", curveFun = "LIN", intrinsic = FALSE, records = 1:9, y_model = "LGCM", TVC = "Rs", decompose = 0)
  set.seed(20191029)
  fit <- suppressWarnings(suppressMessages(do.call(getMGroup, c(list(dat = dat, grp_var = "SEX", sub_Model = "TVC"), args))))@mxOutput
  fig <- suppressWarnings(getFigure(fit, nClass = 2, grp_var = "SEX", sub_Model = "TVC", y_var = "M", curveFun = "LIN", y_model = "LGCM",
                                    t_var = "T", records = 1:9, xstarts = 0))
  p <- fig@figures[[1]][[1]]
  labels <- ggplot2::ggplot_build(p)$plot$scales$get_scales("linetype")$get_labels()
  expect_true(all(c("Model Implied Growth Component (TVC Effect Excluded)", "Model Implied Mean Outcome (TVC Effect Included)") %in% labels))
  full <- p$layers[[5]]$data
  par <- OpenMx::omxGetParameters(fit)
  # With a direct TVC, each group's full mean is its growth component plus its kappa times its TVC means,
  # at the group's own mean occasions.
  for (k in 1:2) {
    times <- colMeans(fit@submodels[[k]]@data$observed[, paste0("T", 1:9)])
    expect_equal(full$time[full$Class == k], unname(times))
    mean0 <- c(OpenMx::mxEvalByName(paste0("c", k, "Y_mean0"), fit@submodels[[k]]))
    expect_equal(full$value[full$Class == k],
                 unname(mean0[1] + mean0[2] * times + par[paste0("c", k, "kappa")] * par[paste0("c", k, "TVC_m", 1:9)]), tolerance = 1e-8)
  }
})

test_that("groups with different schedules each get their own full means and an observed curve", {
  skip_on_cran()
  OpenMx::mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)
  set.seed(1027); n <- 400; g <- rep(1:2, each = n / 2)
  b0 <- 10 + rnorm(n, sd = 2); b1 <- 1 + rnorm(n, sd = 0.3)
  dat <- data.frame(id = 1:n, group = g)
  for (j in 1:6) {
    t <- (j - 1) * ifelse(g == 1, 1, 2); x <- rnorm(n, mean = 0.7 * t, sd = 1)
    dat[[paste0("T", j)]] <- t; dat[[paste0("X", j)]] <- x; dat[[paste0("Y", j)]] <- b0 + b1 * t + 0.5 * x + rnorm(n, sd = 0.5)
  }
  fit <- suppressWarnings(suppressMessages(getMGroup(dat = dat, grp_var = "group", sub_Model = "TVC", t_var = "T", y_var = "Y",
    records = 1:6, curveFun = "LIN", intrinsic = FALSE, TVC = "X", y_model = "LGCM", decompose = 0)))
  expect_identical(getDiagnostics(fit)$issues, character(0))
  p <- plot(fit, draw = FALSE)@figures[[1]][[1]]
  full <- Filter(function(l) inherits(l$geom, "GeomPoint"), p$layers)[[1]]$data
  expect_equal(full$time[full$Class == 1], 0:5)
  expect_equal(full$time[full$Class == 2], seq(0, 10, 2))
  for (k in 1:2) {
    expect_equal(full$value[full$Class == k], nlpsem:::.tvc_full_mean(fit@mxOutput@submodels[[k]], paste0("Y", 1:6)), tolerance = 1e-10)
  }
  # 1,200 observations per group at six occasions: the observed curve is drawn, without warnings
  observed <- which(vapply(p$layers, function(l) inherits(l$stat, "StatSummary"), logical(1)))
  expect_length(observed, 1)
  built <- expect_no_warning(suppressMessages(ggplot2::ggplot_build(p)))
  expect_equal(nrow(built$data[[observed]]), 12)
})

test_that("the observed curve falls back to occasion means only where the default smoother cannot fit", {
  build <- function(time) {
    d <- data.frame(time = time, value = time + rnorm(length(time)), Class = rep(1:2, length.out = length(time)))
    p <- ggplot2::ggplot() + nlpsem:::.observed_curve(d, ggplot2::aes(x = time, y = value, group = Class))
    list(stat = class(p$layers[[1]]$stat)[1], rows = nrow(suppressMessages(ggplot2::ggplot_build(p))$data[[1]]))
  }
  set.seed(3)
  fixed <- build(rep(0:4, length.out = 2400))
  expect_identical(fixed$stat, "StatSummary"); expect_equal(fixed$rows, 10)
  expect_identical(build(rep(0:4, length.out = 400))$stat, "StatSmooth")
  varying <- build(runif(2400, 0, 4))
  expect_identical(varying$stat, "StatSmooth"); expect_gt(varying$rows, 0)
})

test_that("observed curves count only occasions with usable outcomes", {
  d <- data.frame(time = rep(0:9, each = 200), value = rep(0:9, each = 200))
  d$value[d$time == 9] <- NA_real_
  p <- ggplot2::ggplot() + nlpsem:::.observed_curve(d, ggplot2::aes(x = time, y = value))
  built <- expect_no_warning(suppressMessages(ggplot2::ggplot_build(p)))
  expect_equal(built$data[[1]]$x, 0:8)
  expect_equal(built$data[[1]]$y, 0:8)
  d$value[] <- NA_real_
  p <- ggplot2::ggplot() + nlpsem:::.observed_curve(d, ggplot2::aes(x = time, y = value))
  built <- expect_no_warning(suppressMessages(ggplot2::ggplot_build(p)))
  expect_equal(nrow(built$data[[1]]), 0)
})

test_that("a grouping column named ID retains its group identities in plots", {
  fit <- .load_extdata("getMGroup_examples.RData")$MGroup_Math_BLS_LGCM_f@mxOutput
  figure <- function(model, grp_var) getFigure(model, nClass = 2, grp_var = grp_var,
    sub_Model = "LGCM", y_var = "M", curveFun = "BLS", t_var = "T", records = 1:9,
    xstarts = 0, outcome = "Math")@figures[[1]][[1]]
  original <- figure(fit, "SEX")
  for (k in seq_along(fit@submodels)) {
    dat <- fit@submodels[[k]]@data@observed
    dat$ID <- dat$SEX
    dat$SEX <- NULL
    fit@submodels[[k]]@data@observed <- dat
  }
  renamed <- figure(fit, "ID")
  expect_equal(renamed$layers[[4]]$data, original$layers[[4]]$data)
  expect_equal(renamed$scales$get_scales("colour")$labels, original$scales$get_scales("colour")$labels)
  built <- expect_no_warning(suppressMessages(ggplot2::ggplot_build(renamed)))
  expect_equal(length(unique(built$data[[4]]$group)), 2)
})
