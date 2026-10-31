# Regenerates inst/extdata/getMediation_examples.RData from the package data with the current code.
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
xstarts <- mean(baseT)

Med2_LGCM_BLS <- getMediation(
  dat = RMS_dat0, t_var = rep("T", 2), y_var = "M", m_var = "R",
  x_type = "baseline", x_var = "ex1", curveFun = "bilinear spline",
  records = list(1:9, 1:9), tries = 10,
  paramOut = TRUE
)
set.seed(20191029)

Med3_LGCM_BLS <- getMediation(
  dat = RMS_dat0, t_var = rep("T", 3), y_var = "S", m_var = "M", x_type = "longitudinal",
  x_var = "R", curveFun = "bilinear spline", records = list(2:9, 1:9, 1:9),
  tries = 10, paramOut = TRUE
)

save(baseT, Med2_LGCM_BLS, Med3_LGCM_BLS,
     RMS_dat, RMS_dat0, xstarts,
     file = file.path(out_dir, "getMediation_examples.RData"), compress = "xz")
message("getMediation_examples.RData written to ", out_dir)
