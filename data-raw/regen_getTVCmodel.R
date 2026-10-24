# Regenerates inst/extdata/getTVCmodel_examples.RData from the package data with the current code.
# Run through data-raw/regenerate_extdata.R (which loads the package and sets `out_dir`),
# or interactively after `pkgload::load_all()` with `out_dir` set to a writable folder.
if (!exists("out_dir")) out_dir <- "inst/extdata"
mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)

data("RMS_dat")
RMS_dat0 <- RMS_dat
baseT <- RMS_dat0$T1
RMS_dat0$T1 <- (RMS_dat0$T1 - baseT)/12
RMS_dat0$T2 <- (RMS_dat0$T2 - baseT)/12
RMS_dat0$T3 <- (RMS_dat0$T3 - baseT)/12
RMS_dat0$T4 <- (RMS_dat0$T4 - baseT)/12
RMS_dat0$T5 <- (RMS_dat0$T5 - baseT)/12
RMS_dat0$T6 <- (RMS_dat0$T6 - baseT)/12
RMS_dat0$T7 <- (RMS_dat0$T7 - baseT)/12
RMS_dat0$T8 <- (RMS_dat0$T8 - baseT)/12
RMS_dat0$T9 <- (RMS_dat0$T9 - baseT)/12
RMS_dat0$ex1 <- scale(RMS_dat0$Approach_to_Learning)
BL_mean <- mean(RMS_dat0[, "R1"])
BL_var <- var(RMS_dat0[, "R1"])
RMS_dat0$Rs1 <- (RMS_dat0$R1 - BL_mean)/sqrt(BL_var)
RMS_dat0$Rs2 <- (RMS_dat0$R2 - BL_mean)/sqrt(BL_var)
RMS_dat0$Rs3 <- (RMS_dat0$R3 - BL_mean)/sqrt(BL_var)
RMS_dat0$Rs4 <- (RMS_dat0$R4 - BL_mean)/sqrt(BL_var)
RMS_dat0$Rs5 <- (RMS_dat0$R5 - BL_mean)/sqrt(BL_var)
RMS_dat0$Rs6 <- (RMS_dat0$R6 - BL_mean)/sqrt(BL_var)
RMS_dat0$Rs7 <- (RMS_dat0$R7 - BL_mean)/sqrt(BL_var)
RMS_dat0$Rs8 <- (RMS_dat0$R8 - BL_mean)/sqrt(BL_var)
RMS_dat0$Rs9 <- (RMS_dat0$R9 - BL_mean)/sqrt(BL_var)
xstarts <- mean(baseT)

set.seed(20191029)
Math_TVC_BLS_f <- getTVCmodel(
  dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS", intrinsic = TRUE,
  records = 1:9, y_model = "LGCM", TVC = "Rs", decompose = 0, growth_TIC = "ex1",
  tries = 10
)
set.seed(20191029)
Math_TVCslp_BLS_f <- getTVCmodel(
  dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS", intrinsic = TRUE,
  records = 1:9, y_model = "LGCM", TVC = "Rs", decompose = 1, growth_TIC = "ex1",
  tries = 10, paramOut = TRUE
)
set.seed(20191029)
Math_TVCslp_BLS_r <- getTVCmodel(
  dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS", intrinsic = FALSE,
  records = 1:9, y_model = "LGCM", TVC = "R", decompose = 1, growth_TIC = "ex1",
  tries = 10, paramOut = TRUE)
set.seed(20191029)
Math_TVCchg_BLS_r <- getTVCmodel(
  dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS", intrinsic = FALSE,
  records = 1:9, y_model = "LGCM", TVC = "R", decompose = 2, growth_TIC = "ex1",
  tries = 10, paramOut = TRUE)
set.seed(20191029)
Math_TVCchgBL_BLS_r <- getTVCmodel(
  dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS", intrinsic = FALSE,
  records = 1:9, y_model = "LGCM", TVC = "R", decompose = 3, growth_TIC = "ex1",
  tries = 10, paramOut = TRUE)

save(baseT, BL_mean, BL_var, Math_TVC_BLS_f, Math_TVCchg_BLS_r,
     Math_TVCchgBL_BLS_r, Math_TVCslp_BLS_f, Math_TVCslp_BLS_r,
     RMS_dat, RMS_dat0, xstarts,
     file = file.path(out_dir, "getTVCmodel_examples.RData"), compress = "xz")
message("getTVCmodel_examples.RData written to ", out_dir)
