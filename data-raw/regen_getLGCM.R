# Regenerates inst/extdata/getLGCM_examples.RData from the package data with the current code.
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
xstarts <- mean(baseT)

Math_LGCM_BLS_f <- getLGCM(
  dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "bilinear spline",
  intrinsic = TRUE, records = 1:9, growth_TIC = NULL, tries = 10
)
Math_LGCM_BLS_r <- getLGCM(
  dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "bilinear spline",
  intrinsic = FALSE, records = 1:9, growth_TIC = NULL, tries = 10
)
Math_LGCM_TIC_BLS_f <- getLGCM(
  dat = RMS_dat0, t_var = "T", y_var = "M", curveFun = "bilinear spline",
  intrinsic = TRUE, records = 1:9, growth_TIC = c("ex1", "ex2"),
  tries = 10, paramOut = TRUE
)

save(baseT, Math_LGCM_BLS_f, Math_LGCM_BLS_r, Math_LGCM_TIC_BLS_f,
     RMS_dat, RMS_dat0, xstarts,
     file = file.path(out_dir, "getLGCM_examples.RData"), compress = "xz")
message("getLGCM_examples.RData written to ", out_dir)
