# Regenerate every fitted example object shipped in inst/extdata with the current package code.
#
#   Rscript data-raw/regenerate_extdata.R            # writes into inst/extdata
#   Rscript data-raw/regenerate_extdata.R /tmp/out   # writes into another folder for comparison first
#
# Each script fits several models with mxTryHard (10 to 50 extra attempts), so a full run takes a long
# time; the scripts are independent and can be run in parallel. The seeds inside the scripts make the
# jittered attempts reproducible for a given OpenMx version.
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else "inst/extdata"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
pkgload::load_all(".", quiet = TRUE)
suppressPackageStartupMessages(library(OpenMx))
scripts <- c("regen_getLGCM.R", "regen_getLCSM.R", "regen_getTVCmodel.R", "regen_getMGM.R",
             "regen_getMediation.R", "regen_getMGroup.R", "regen_getMIX.R")
for (s in scripts) {
  message("== ", s, " (", format(Sys.time(), "%H:%M:%S"), ")")
  env <- new.env(); env$out_dir <- out_dir
  source(file.path("data-raw", s), local = env)
}
writeLines(c(paste("Regenerated:", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
             paste("R:", R.version.string), paste("OpenMx:", as.character(packageVersion("OpenMx"))),
             paste("nlpsem:", as.character(packageVersion("nlpsem")))),
           file.path(out_dir, "REGENERATED.txt"))
message("done; see ", file.path(out_dir, "REGENERATED.txt"))
