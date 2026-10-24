# Regenerates inst/extdata/getMGroup_examples.RData from the package data with the current code.
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

set.seed(20191029)
MGroup_Math_BLS_LGCM_f <- getMGroup(
  dat = RMS_dat0, grp_var = "SEX", sub_Model = "LGCM", y_var = "M", t_var = "T",
  records = 1:9, curveFun = "BLS", intrinsic = TRUE, tries = 20
)
set.seed(20191029)
MGroup_Read_EXP_LGCM_TIC_r <- getMGroup(
  dat = RMS_dat0, grp_var = "SEX", sub_Model = "LGCM", y_var = "R", t_var = "T",
  records = 1:9, curveFun = "EXP", intrinsic = FALSE,
  growth_TIC = c("ex1", "ex2"), paramOut = TRUE
)

save(MGroup_Math_BLS_LGCM_f, MGroup_Read_EXP_LGCM_TIC_r, xstarts,
     file = file.path(out_dir, "getMGroup_examples.RData"), compress = "xz")
message("getMGroup_examples.RData written to ", out_dir)
