# Regenerates inst/extdata/getLCSM_examples.RData from the package data with the current code.
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
RMS_dat0$ex2 <- scale(RMS_dat0$Attention_focus)
xstarts <- mean(baseT)/12

Read_LCSM_NonP <- getLCSM(
  dat = RMS_dat0, t_var = "T", y_var = "R", curveFun = "nonparametric",
  intrinsic = FALSE, records = 1:9, growth_TIC = NULL,
  tries = 10, paramOut = TRUE
)
Read_LCSM_NonP_TIC <- getLCSM(
  dat = RMS_dat0, t_var = "T", y_var = "R", curveFun = "nonparametric",
  intrinsic = FALSE, records = 1:9, growth_TIC = c("ex1", "ex2"),
  tries = 10, paramOut = TRUE
)
Read_LCSM_QUAD <- getLCSM(
  dat = RMS_dat0, t_var = "T", y_var = "R", curveFun = "quadratic", intrinsic = FALSE,
  records = 1:9, tries = 10, paramOut = TRUE
)
set.seed(20191029)
Read_LCSM_EXP_r <- getLCSM(
  dat = RMS_dat0, t_var = "T", y_var = "R", curveFun = "negative exponential",
  intrinsic = FALSE, records = 1:9, tries = 10, paramOut = TRUE
)
set.seed(20191029)
Read_LCSM_JB_r <- getLCSM(
  dat = RMS_dat0, t_var = "T", y_var = "R", curveFun = "Jenss-Bayley",
  intrinsic = FALSE, records = 1:9, tries = 10, paramOut = TRUE
)

save(baseT, Read_LCSM_EXP_r, Read_LCSM_JB_r, Read_LCSM_NonP,
     Read_LCSM_NonP_TIC, Read_LCSM_QUAD, RMS_dat, RMS_dat0, xstarts,
     file = file.path(out_dir, "getLCSM_examples.RData"), compress = "xz")
message("getLCSM_examples.RData written to ", out_dir)
