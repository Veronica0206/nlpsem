# Regenerates inst/extdata/getMIX_examples.RData from the package data with the current code.
# Run through data-raw/regenerate_extdata.R (which loads the package and sets `out_dir`),
# or interactively after `pkgload::load_all()` with `out_dir` set to a writable folder.
if (!exists("out_dir")) out_dir <- "inst/extdata"
mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)

data("RMS_dat")
RMS_dat0 <- RMS_dat
baseT <- RMS_dat0$T1
RMS_dat0$T1 <- RMS_dat0$T1 - baseT
RMS_dat0$T2 <- RMS_dat0$T2 - baseT
RMS_dat0$T3 <- RMS_dat0$T3 - baseT
RMS_dat0$T4 <- RMS_dat0$T4 - baseT
RMS_dat0$T5 <- RMS_dat0$T5 - baseT
RMS_dat0$T6 <- RMS_dat0$T6 - baseT
RMS_dat0$T7 <- RMS_dat0$T7 - baseT
RMS_dat0$T8 <- RMS_dat0$T8 - baseT
RMS_dat0$T9 <- RMS_dat0$T9 - baseT
RMS_dat0$ex1 <- scale(RMS_dat0$Approach_to_Learning)
RMS_dat0$ex2 <- scale(RMS_dat0$Attention_focus)
RMS_dat0$gx1 <- scale(RMS_dat0$INCOME)
RMS_dat0$gx2 <- scale(RMS_dat0$EDU)
xstarts <- mean(baseT)

Math_BLS_LGCM1 <- getLGCM(
  dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "BLS", intrinsic = FALSE,
  records = 1:9, tries = 10
)
Math_BLS_LGCM2 <- getMIX(
  dat = RMS_dat0, prop_starts = c(0.45, 0.55), sub_Model = "LGCM", y_var = "M",
  t_var = "T", records = 1:9, curveFun = "BLS", intrinsic = FALSE, tries = 10
)
set.seed(20191029)
Math_BLS_LGCM3 <- getMIX(
  dat = RMS_dat0, prop_starts = c(0.30, 0.40, 0.30), sub_Model = "LGCM", y_var = "M",
  t_var = "T", records = 1:9, curveFun = "BLS", intrinsic = FALSE, tries = 10
)
set.seed(20191029)
RM_BLS_PLGCM3 <- getMIX(
  dat = RMS_dat0, prop_starts = c(0.30, 0.40, 0.30), sub_Model = "MGM",
  cluster_TIC = c("gx1", "gx2"), t_var = c("T", "T"), y_var = c("R", "M"),
  curveFun = "BLS", intrinsic = FALSE, records = list(1:9, 1:9),
  y_model = "LGCM", tries = 10, paramOut = TRUE
)

save(baseT, Math_BLS_LGCM1, Math_BLS_LGCM2, Math_BLS_LGCM3,
     RM_BLS_PLGCM3, RMS_dat, RMS_dat0, xstarts,
     file = file.path(out_dir, "getMIX_examples.RData"), compress = "xz")
message("getMIX_examples.RData written to ", out_dir)
