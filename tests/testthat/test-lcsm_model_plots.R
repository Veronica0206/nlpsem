# Check the quantities displayed by the plots against the RAM expectation, not
# another implementation of the continuous curve or its analytic derivatives.
.lcsm_plot_case <- function(curve, random = FALSE, tic = FALSE, records = 1:4){
  means <- switch(curve, NonP = c(50, 2), QUAD = c(50, 2, -.01),
                  EXP = c(50, 10, .2), JB = c(50, 2, -10, -.2))
  covariance <- diag(rep(.01, length(means)))
  starts <- list(Y_starts = list(means = means, covMatrix = covariance[lower.tri(covariance, diag = TRUE)],
                                 residuals = 1, rel_rate = c(1, 2, .5)),
                 TIC_starts = list(means = 2, covMatrix = 1),
                 beta = matrix(seq(.05, by = .03, length.out = length(means)), ncol = 1))
  dat <- data.frame(ID = 1:8, x = seq(1, 3, length.out = 8), grp = rep(1:2, each = 4))
  for (j in seq_along(records)){
    dat[[paste0("T", records[j])]] <- c(0, 10, 20, 35)[j]
    dat[[paste0("Y", records[j])]] <- 50 + j + seq_len(nrow(dat)) / 10
    dat[[paste0("Z", records[j])]] <- 10 + j + seq_len(nrow(dat)) / 10
    dat[[paste0("X", records[j])]] <- j + seq_len(nrow(dat)) / 20
  }
  args <- list(dat = dat, t_var = "T", y_var = "Y", curveFun = curve, intrinsic = random,
                records = records, growth_TIC = if (tic) "x" else NULL, starts = starts)
  list(model = do.call(nlpsem:::getLCSM.mxModel, args), starts = starts, dat = dat, args = args)
}

.lcsm_ram_change <- function(model, columns){
  expected <- c(OpenMx::mxGetExpected(model, "means", defvar.row = 1)[1, columns])
  unname(expected - expected[1])
}

test_that("LCSM plotted recurrence agrees with RAM for coarse schedules and shifted rate means", {
  for (curve in c("NonP", "QUAD", "EXP", "JB")){
    for (random in if (curve %in% c("EXP", "JB")) c(FALSE, TRUE) else FALSE){
      for (tic in c(FALSE, TRUE)){
        case <- .lcsm_plot_case(curve, random, tic, records = 2:5)
        out <- nlpsem:::.lcsm_plot_components(case$model, "Y", "Y", curve, "LCSM")
        expected <- .lcsm_ram_change(case$model, paste0("Y", 2:5))
        expect_equal(out$est, expected, tolerance = 1e-10, info = paste(curve, random, tic))
        expect_equal(c(OpenMx::mxEvalByName("Ychg_bl_m",case$model,compute=TRUE)),expected[-1],tolerance=1e-10)
        expect_equal(c(OpenMx::mxEvalByName("Yslp_m",case$model,compute=TRUE)),diff(expected)/diff(c(0,10,20,35)),tolerance=1e-10)
        evaluated <- OpenMx::mxModel(case$model, OpenMx::mxAlgebraFromString(out$expression, name = "plot_check"))
        expect_equal(c(OpenMx::mxEvalByName("plot_check", evaluated, compute = TRUE)), expected, tolerance = 1e-10)
        if (curve == "EXP" && !tic) expect_equal(out$est[2], 7.35758882342885, tolerance = 1e-10)
      }
    }
  }
})

test_that("multivariate and TVC growth-component recurrences agree with RAM", {
  for (curve in c("NonP", "QUAD", "EXP", "JB")){
    for (random in if (curve %in% c("EXP", "JB")) c(FALSE, TRUE) else FALSE){
      case <- .lcsm_plot_case(curve, random, tic = TRUE)
      multi_starts <- list(means = rep(list(case$starts$Y_starts$means), 2),
                           covMatrix = diag(rep(.01, length(case$starts$Y_starts$means)*2)),
                           residuals = diag(2), rel_rate = rep(list(c(1,2,.5)),2))
      multi <- nlpsem:::getMGM.mxModel(case$dat, c("T","T"), c("Y","Z"), curve, random,
                                     list(1:4,1:4), "LCSM", multi_starts)
      for (prefix in c("Y", "Z")){
        out <- nlpsem:::.lcsm_plot_components(multi, prefix, prefix, curve, "MGM")
        expect_equal(out$est, .lcsm_ram_change(multi, paste0(prefix,1:4)), tolerance = 1e-10)
      }
      tvc_starts <- case$starts
      tvc_starts$BL_starts <- case$starts$TIC_starts
      tvc_starts$TVC_starts <- list(means = 1:4, covMatrix = rep(1,4))
      tvc_starts$kappa <- .5
      tvc <- nlpsem:::getTVC.mxModel(case$dat,"T","Y",curve,random,1:4,"LCSM","X",0,"x",tvc_starts)
      out <- nlpsem:::.lcsm_plot_components(tvc,"Y","Y",curve,"TVC")
      # Remove the direct time-varying contribution to isolate the documented
      # growth component. The covariate paths to its growth factors remain.
      component <- OpenMx::omxSetParameters(tvc, labels = "kappa", values = 0)
      expect_equal(out$est, .lcsm_ram_change(component,paste0("Y",1:4)), tolerance = 1e-10)
      expect_equal(c(OpenMx::mxEvalByName("Ychg_bl_m",tvc,compute=TRUE)),out$est[-1],tolerance=1e-10)
    }
  }
})

test_that("decomposed TVC mean rate deviations enter reported rates and changes", {
  for (curve in c("EXP","JB")) for (decompose in 1:3) for (tic in c(FALSE,TRUE)){
    case <- .lcsm_plot_case(curve,random=TRUE,tic=tic)
    starts <- case$starts
    starts$TVC_starts <- list(means=c(1,.3),covMatrix=c(1,0,.1),residuals=1,rel_rate=c(1,2,.5))
    starts$BL_starts <- list(means=if(tic) c(2,1) else 1,covMatrix=if(tic) diag(2) else 1)
    starts$beta <- matrix(seq(.03,by=.01,length.out=length(starts$Y_starts$means)*(1+tic)),ncol=1+tic)
    starts$kappa <- .5; starts$XY_cov <- .1
    model <- nlpsem:::getTVC.mxModel(case$dat,"T","Y",curve,TRUE,1:4,"LCSM","X",decompose,if(tic) "x" else NULL,starts)
    out <- nlpsem:::.lcsm_plot_components(model,"Y","Y",curve,"TVC")
    component <- OpenMx::omxSetParameters(model,labels="kappa",values=0)
    expected <- .lcsm_ram_change(component,paste0("Y",1:4))
    expect_equal(out$est,expected,tolerance=1e-10)
    expect_equal(c(OpenMx::mxEvalByName("Ychg_bl_m",model,compute=TRUE)),expected[-1],tolerance=1e-10)
    expect_equal(c(OpenMx::mxEvalByName("Yslp_m",model,compute=TRUE)),diff(expected)/diff(c(0,10,20,35)),tolerance=1e-10)
  }
})

test_that("class-specific fitted recurrences retain their own grid and rate shifts", {
  for (container in c("multigroup", "mixture")){
    case <- .lcsm_plot_case("EXP", random = TRUE, tic = TRUE)
    dat <- case$dat
    for (j in 1:4) dat[dat$grp == 2, paste0("T",j)] <- dat[dat$grp == 2,paste0("T",j)] + 2
    classes <- nlpsem:::getsub.LCSM(dat,2,"T","Y","EXP",TRUE,1:4,"x",rep(list(case$starts),2),
                                  container,if(container == "multigroup") "grp" else NULL)
    for (k in 1:2){
      observed <- if (container == "multigroup") dat[dat$grp == k,] else dat
      for (j in 1:4) observed[[paste0("T",j)]] <- mean(observed[[paste0("T",j)]])
      eval_model <- OpenMx::mxModel(classes[[k]], OpenMx::mxData(observed, type="raw"))
      out <- nlpsem:::.lcsm_plot_components(eval_model,paste0("c",k,"Y"),paste0("Class",k,".c",k,"Y"),"EXP","LCSM")
      expect_equal(out$est,.lcsm_ram_change(eval_model,paste0("Y",1:4)),tolerance=1e-10)
    }
  }
})

test_that("public LCSM subsets preserve endpoints, baseline and joint contrast uncertainty", {
  skip_on_cran()   # 24 public subsets with joint contrast standard errors (about 4 s)
  e <- .load_extdata("getLCSM_examples.RData")
  cases <- c(NonP="Read_LCSM_NonP", QUAD="Read_LCSM_QUAD", EXP="Read_LCSM_EXP_r", JB="Read_LCSM_JB_r")
  for (curve in names(cases)){
    fit <- e[[cases[[curve]]]]@mxOutput
    full <- c(0, OpenMx::mxEvalByName("Ychg_bl_m",fit))
    for (records in list(1:9,2:9,c(1,3,5,7,9),c(2,4,6,8,9),c(2,9),9:2)){
      figure <- getFigure(fit,sub_Model="LCSM",y_var="R",curveFun=curve,t_var="T",records=records,xstarts=0)
      actual <- figure@figures[[1]][[1]]$data
      times <- unname(colMeans(fit$data$observed[,paste0("T",records),drop=FALSE]))
      expect_equal(nrow(actual),length(records)-1L)
      expect_equal(unname(actual$time),times[-1])
      expect_equal(actual$cY_hat,full[records[-1]]-full[records[1]],tolerance=1e-10)
      expect_equal(actual$dY_hat,diff(full[records])/diff(times),tolerance=1e-10)
      # Joint contrast uncertainty: the covariance with a selected baseline
      # contributes and cannot be recovered from separate marginal intervals.
      term <- function(j) if(j == 1) "0" else paste0("Ychg_bl_m[",j-1,",1]")
      expressions <- vapply(records[-1],function(j) paste0("(",term(j),"-",term(records[1]),")"),character(1))
      expected_se <- c(OpenMx::mxSE(paste0("rbind(",paste(expressions,collapse=","),")"),fit,silent=TRUE))
      actual_se <- (actual$cY_hat_u-actual$cY_hat_l)/(2*qnorm(.975))
      expect_equal(actual_se,expected_se,tolerance=3e-4)
    }
  }
})

test_that("legacy aliases and invalid selections cannot recycle fitted changes", {
  fit <- .load_extdata("getLCSM_examples.RData")$Read_LCSM_NonP@mxOutput
  for (j in 1:9){
    fit$data$observed[[paste0("Outcome",j)]] <- fit$data$observed[[paste0("R",j)]]
    fit$data$observed[[paste0("Time",j)]] <- fit$data$observed[[paste0("T",j)]]
  }
  figure <- suppressWarnings(getFigure(fit,sub_Model="LCSM",y_var="Outcome",curveFun="NonP",t_var="Time",records=c(2,9),xstarts=0))
  expect_equal(nrow(figure@figures[[1]][[1]]$data),1L)
  expect_equal(figure@figures[[1]][[1]]$data$cY_hat,
               c(OpenMx::mxEvalByName("Ychg_bl_m",fit))[8]-c(OpenMx::mxEvalByName("Ychg_bl_m",fit))[1])
  call <- function(records) getFigure(fit,sub_Model="LCSM",y_var="R",curveFun="NonP",t_var="T",records=records,xstarts=0)
  expect_error(call(c(2,2,9)),"distinct")
  expect_error(call(2),"at least two")
  expect_error(call(c(1,5,3)),"fitted order")
  fit$data$observed$R10 <- fit$data$observed$R9
  fit$data$observed$T10 <- fit$data$observed$T9+1
  expect_error(call(c(1,10)),"included in the fitted model")
})

test_that("random-rate TIC contrast covariance agrees with numerical RAM derivatives", {
  case <- .lcsm_plot_case("EXP",random=TRUE,tic=TRUE)
  model <- case$model
  components <- nlpsem:::.lcsm_plot_components(model,"Y","Y","EXP","LCSM")
  parameters <- OpenMx::omxGetParameters(model)
  pcov <- diag(length(parameters)) * .001
  pcov[1,2] <- pcov[2,1] <- -.0003
  dimnames(pcov) <- list(names(parameters),names(parameters))
  delta_cov <- OpenMx::mxSE(components$expression,model,cov=pcov,details=TRUE,silent=TRUE)$Cov
  jacobian <- vapply(seq_along(parameters),function(j){
    step <- 1e-5 * max(1,abs(parameters[j]))
    plus <- OpenMx::omxSetParameters(model,labels=names(parameters)[j],values=parameters[j]+step)
    minus <- OpenMx::omxSetParameters(model,labels=names(parameters)[j],values=parameters[j]-step)
    (.lcsm_ram_change(plus,paste0("Y",1:4))-.lcsm_ram_change(minus,paste0("Y",1:4)))/(2*step)
  },numeric(4))
  expect_equal(unname(delta_cov),unname(jacobian %*% pcov %*% t(jacobian)),tolerance=3e-4)
})

test_that("nonparametric rate segments occupy the intervals they summarize", {
  fit <- .load_extdata("getLCSM_examples.RData")$Read_LCSM_NonP@mxOutput
  records <- c(2,4,6,9)
  p <- getFigure(fit,sub_Model="LCSM",y_var="R",curveFun="NonP",t_var="T",records=records,xstarts=10)@figures[[1]][[2]]
  built <- ggplot2::ggplot_build(p)
  times <- unname(colMeans(fit$data$observed[,paste0("T",records)]))+10
  full <- c(0,OpenMx::mxEvalByName("Ychg_bl_m",fit))
  rate <- diff(full[records])/diff(times)
  expect_equal(built$data[[1]]$x,times[-length(times)])
  expect_equal(built$data[[1]]$xend,times[-1])
  expect_equal(built$data[[1]]$y,rate)
  expect_equal(built$data[[1]]$yend,rate)
  expect_equal(built$data[[2]]$x,times[-length(times)])
  expect_equal(built$data[[3]]$xend,times[-1])
})

test_that("public group and mixture subsets use the parent fitted covariance", {
  skip_on_cran()
  OpenMx::mxOption(NULL,"Default optimizer","CSOLNP")
  dat <- getSimData(n=160,nT=4,curveFun="LIN",times=c(0,1,3,5),jitter=0,seed=416)
  dat$grp <- rep(1:2,each=80)
  group <- suppressMessages(getMGroup(dat,grp_var="grp",sub_Model="LCSM",t_var="T",y_var="Y",curveFun="NonP",intrinsic=FALSE,records=1:4))
  mixture <- suppressMessages(getMIX(dat,prop_starts=1,sub_Model="LCSM",t_var="T",y_var="Y",curveFun="NonP",intrinsic=FALSE,records=1:4))
  for (object in list(group,mixture)){
    fit <- object@mxOutput
    n <- length(fit@submodels)
    grp <- if (n == 2) "grp" else NULL
    out <- getFigure(fit,nClass=n,grp_var=grp,sub_Model="LCSM",t_var="T",y_var="Y",curveFun="NonP",records=c(2,4),xstarts=0)@figures[[1]][[1]]$data
    expect_equal(nrow(out),n)
    for (k in seq_len(n)){
      pfx <- paste0("Class",k,".c",k,"Ychg_bl_m")
      expression <- paste0(pfx,"[3,1]-",pfx,"[1,1]")
      value <- OpenMx::mxEvalByName(paste0("c",k,"Ychg_bl_m"),fit@submodels[[k]])
      se <- c(OpenMx::mxSE(expression,fit,silent=TRUE))
      expect_equal(out$cY_hat[k],c(value)[3]-c(value)[1],tolerance=1e-10)
      expect_equal((out$cY_hat_u[k]-out$cY_hat_l[k])/(2*qnorm(.975)),se,tolerance=3e-4)
    }
    figures <- getFigure(fit,nClass=n,grp_var=grp,sub_Model="LCSM",t_var="T",y_var="Y",curveFun="NonP",records=c(2,4),xstarts=0)@figures[[1]]
    for (j in seq_along(figures)){
      p <- figures[[j]]
      expect_s3_class(p$layers[[1]]$geom,"GeomPointrange")
      built <- suppressMessages(ggplot2::ggplot_build(p))
      prefix <- if(j == 1L) "cY_hat" else "dY_hat"
      expect_equal(nrow(built$data[[1]]),n)
      expect_equal(built$data[[1]]$y,p$data[[prefix]])
      expect_equal(built$data[[1]]$ymin,p$data[[paste0(prefix,"_l")]])
      expect_equal(built$data[[1]]$ymax,p$data[[paste0(prefix,"_u")]])
      expect_true(all(p$data$time %in% built$layout$panel_params[[1]]$x$breaks))
    }
  }
})
