# Reconstruct from the documented local public-use exports without overwriting data/.
# Rscript data-raw/reconstruct_RMS_dat.R <input-directory> <new-output-directory>
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
if (!requireNamespace("haven", quietly = TRUE)) stop("Install haven for SAS import.")
if (!requireNamespace("digest", quietly = TRUE)) stop("Install digest for SHA256 verification.")
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install jsonlite for provenance.")
mapping <- read.csv("data-raw/RMS_dat-variable-map.csv", stringsAsFactors = FALSE)
provenance <- jsonlite::read_json("data-raw/RMS_dat-provenance.json", simplifyVector = TRUE)
if (dir.exists(args[2])) stop("Use a new output directory to preserve prior artifacts.")
frames <- lapply(unique(mapping$file), function(file) {
  path <- file.path(args[1], file)
  expected <- provenance$files$sha256[match(file, provenance$files$file)]
  if (!identical(digest::digest(file = path, algo = "sha256"), expected)) {
    stop("Input checksum differs from the documented export: ", file)
  }
  dat <- as.data.frame(haven::read_sas(path))
  variables <- unique(c("CHILDID", mapping$source[mapping$file == file]))
  if (!all(variables %in% names(dat))) stop("Missing documented columns in ", file)
  if (anyNA(dat$CHILDID) || anyDuplicated(dat$CHILDID)) stop("Invalid subject IDs in ", file)
  dat[, variables, drop = FALSE]
})
joined <- Reduce(function(x, y) merge(x, y, by = "CHILDID", all.x = TRUE, sort = TRUE), frames)
selected <- joined[, mapping$source, drop = FALSE]
names(selected) <- mapping$output
selected[] <- lapply(selected, function(x) {
  if (is.numeric(x)) {
    x <- as.numeric(x)
    x[is.nan(x) | x %in% provenance$reconstruction$missing_codes] <- NA_real_
  }
  x
})
complete <- selected[stats::complete.cases(selected), , drop = FALSE]
stopifnot(nrow(complete) == provenance$reconstruction$complete_cases)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(provenance$reconstruction$seed)
ids <- sample(unique(complete$ID), provenance$reconstruction$sample_n, replace = FALSE)
RMS_dat <- complete[complete$ID %in% ids, , drop = FALSE]
RMS_dat$EDU <- pmax(RMS_dat$P1EDU, RMS_dat$P2EDU)
RMS_dat <- RMS_dat[, !names(RMS_dat) %in% c("P1EDU", "P2EDU"), drop = FALSE]
original <- new.env()
load("data/RMS_dat.RData", envir = original)
comparison <- data.frame(variable = names(RMS_dat), equal = vapply(names(RMS_dat), function(nm) {
  isTRUE(all.equal(RMS_dat[[nm]], original$RMS_dat[[nm]], check.attributes = FALSE))
}, logical(1)))
stopifnot(identical(names(RMS_dat), names(original$RMS_dat)), all(comparison$equal))
dir.create(args[2], recursive = TRUE)
save(RMS_dat, file = file.path(args[2], "RMS_dat.RData"))
write.csv(comparison, file.path(args[2], "column-comparison.csv"), row.names = FALSE)
writeLines(c("All 49 variables and 500 rows reproduce the shipped values.",
             capture.output(sessionInfo())), file.path(args[2], "reconstruction.log"))
cat("Verified reconstruction in", args[2], "\n")
