# Regenerates inst/extdata/getMGM_examples.RData from the package data with the current code.
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
xstarts <- mean(baseT)

RM_PLGCM.r <- getMGM(
  dat = RMS_dat0, t_var = c("T", "T"), y_var = c("R", "M"), curveFun = "BLS",
  intrinsic = FALSE, records = list(1:9, 1:9), y_model = "LGCM",
  tries = 10, paramOut = TRUE
)
RM_PLGCM.f <- getMGM(
  dat = RMS_dat0, t_var = c("T", "T"), y_var = c("R", "M"), curveFun = "BLS",
  intrinsic = TRUE, records = list(1:9, 1:9), y_model = "LGCM",
  tries = 10, paramOut = TRUE
)

save(baseT, RM_PLGCM.f, RM_PLGCM.r, RMS_dat, RMS_dat0, xstarts,
     file = file.path(out_dir, "getMGM_examples.RData"), compress = "xz")
message("getMGM_examples.RData written to ", out_dir)
