# Self-tests for the structural comparison gate. Run from the package root:
# Rscript data-raw/test_structure_comparator.R
# These fixtures check the comparator itself, not model fitting.
source("data-raw/compare_model_structures.R")

constant <- function(rows = 4L) list(
  class = "FullMatrix", values = matrix(1, rows, 1L),
  free = matrix(FALSE, rows, 1L), labels = matrix(NA_character_, rows, 1L),
  lbound = matrix(NA_real_, rows, 1L), ubound = matrix(NA_real_, rows, 1L),
  dimnames = NULL
)
fixture <- function(prefix = "") {
  old <- list(name = "example", matrices = list(ADD1 = constant(), S = constant(2L)),
              algebras = list(), submodels = list())
  for (outcome in c("A", "B")) {
    old$algebras[[paste0(prefix, outcome, "r_loads")]] <- list(
      formula = paste0("cbind(ADD1, ", outcome, "M_MID_TIME[-1, ] * 2)"))
  }
  old$algebras$unchanged <- list(formula = "cbind(1, 2)")
  new <- old
  new$matrices$ADD1 <- NULL
  for (outcome in c("A", "B")) {
    key <- paste0(prefix, outcome, "ADD1")
    new$matrices[[key]] <- constant()
    new$algebras[[paste0(prefix, outcome, "r_loads")]]$formula <-
      paste0("cbind(", key, ", ", outcome, "M_MID_TIME[-1, ] * 2)")
  }
  new$matrices <- new$matrices[order(names(new$matrices))]
  list(old = old, new = new)
}
compare_fixture <- function(old, new, expected_count = 1L, numeric_tolerance = 1e-10) {
  invisible(capture.output(result <- compare_model_structures(
    list(example = old), list(example = new), expected_count = expected_count,
    numeric_tolerance = numeric_tolerance)))
  result
}
checks <- 0L
pass <- function(label, code) {
  force(code)
  checks <<- checks + 1L
  cat("PASS:", label, "\n")
}
reject <- function(label, old, new, pattern) {
  message <- tryCatch({compare_fixture(old, new); NA_character_},
                      error = function(e) conditionMessage(e))
  if (is.na(message) || !grepl(pattern, message, fixed = TRUE)) {
    stop(label, ": expected rejection containing '", pattern, "'; got ", message)
  }
  pass(label, TRUE)
}

f <- fixture()
pass("identical structures are unchanged", {
  x <- compare_fixture(f$old, f$old)
  stopifnot(x$exact == 1L, x$approved_renames == 0L)
})
pass("both correct outcome constants are accepted", {
  x <- compare_fixture(f$old, f$new)
  stopifnot(x$exact == 0L, x$approved_renames == 1L)
})
pass("class prefixes are preserved", {
  g <- fixture("c12")
  stopifnot(compare_fixture(g$old, g$new)$approved_renames == 1L)
})
pass("renames recurse into independently named submodels", {
  a <- fixture("c1"); b <- fixture("c2")
  old <- list(matrices = list(), algebras = list(),
              submodels = list(Class1 = a$old, Class2 = b$old))
  new <- old
  new$submodels <- list(Class1 = a$new, Class2 = b$new)
  x <- compare_fixture(old, new)
  stopifnot(length(x$renamed_scopes) == 2L, x$approved_renames == 1L)
})

bad <- f$new
bad$matrices$BADD1 <- NULL
bad$algebras$Br_loads$formula <- "cbind(AADD1, BM_MID_TIME[-1, ] * 2)"
reject("shared replacement constant cannot masquerade as split", f$old, bad,
       "Outcome-specific ADD1 constant names differ")

bad <- f$new
bad$algebras$Br_loads$formula <- "cbind(AADD1, BM_MID_TIME[-1, ] * 2)"
reject("wrong outcome reference is rejected even when both constants exist", f$old, bad,
       "Incorrect ADD1 outcome mapping or formula")

bad <- f$new
bad$matrices$CADD1 <- constant()
reject("extra renamed constant is rejected", f$old, bad,
       "Outcome-specific ADD1 constant names differ")

bad <- f$new
bad$matrices$AADD1$values[1, 1] <- 2
reject("changed constant value is rejected", f$old, bad, "ADD1 constants differ in content")

bad <- f$new
bad$matrices$AADD1$free[1, 1] <- TRUE
reject("changed constant free mask is rejected", f$old, bad, "ADD1 constants differ in content")

bad <- f$new
bad$matrices$AADD1$labels[1, 1] <- "unexpected"
reject("changed constant parameter label is rejected", f$old, bad, "ADD1 constants differ in content")

bad <- f$new
bad$matrices$AADD1 <- constant(3L)
reject("changed constant dimensions are rejected", f$old, bad, "ADD1 constants differ in content")

bad <- f$new
bad$algebras$Ar_loads$formula <- "cbind(AADD1, AM_MID_TIME[-1, ] * 3)"
reject("formula changes beyond the renamed symbol are rejected", f$old, bad,
       "Incorrect ADD1 outcome mapping or formula")

bad <- f$new
bad$algebras$Ar_loads$formula <- "cbind(AADD1suffix, AM_MID_TIME[-1, ] * 2)"
reject("symbol substrings cannot be normalized", f$old, bad,
       "Incorrect ADD1 outcome mapping or formula")

bad <- f$new
bad$algebras$unchanged$formula <- "cbind(AADD1, 2)"
reject("unaffected formulas are not normalized", f$old, bad,
       "Unexpected structural differences")

bad <- f$new
bad$matrices$S$free[1, 1] <- TRUE
reject("unrelated model constraints are preserved exactly", f$old, bad,
       "Unexpected structural differences")

bad <- f$new
bad$algebras$Ar_loads$extra <- "unexpected"
reject("other algebra fields are preserved exactly", f$old, bad,
       "Unexpected structural differences")

bad <- f$new
bad$submodels <- list(extra = f$new)
reject("unexpected submodels are rejected", f$old, bad, "Submodel names differ")

bad <- f$new
bad$error <- "construction failed"
reject("model construction errors are rejected", f$old, bad, "Construction error")

pass("default gate requires the complete 688-model grid", {
  message <- tryCatch({compare_model_structures(list(example = f$old), list(example = f$new)); ""},
                      error = function(e) conditionMessage(e))
  stopifnot(grepl("Expected 688 model structures", message, fixed = TRUE))
})
pass("model identifiers cannot silently change", {
  message <- tryCatch({compare_model_structures(list(example = f$old), list(other = f$new),
                                                expected_count = 1L); ""},
                      error = function(e) conditionMessage(e))
  stopifnot(grepl("Model names differ", message, fixed = TRUE))
})

pass("tiny value roundoff is reported separately from exact agreement", {
  new <- f$old
  new$matrices$S$values[1, 1] <- new$matrices$S$values[1, 1] + 5e-11
  x <- compare_fixture(f$old, new)
  stopifnot(x$exact == 0L, x$approved_renames == 0L, x$roundoff_only == 1L,
            x$roundoff_equivalent == 1L, x$max_abs_roundoff > 0,
            x$max_abs_roundoff <= 1e-10, x$max_scaled_roundoff <= 1e-10,
            length(x$roundoff_fields) == 1L)
})
pass("combined renames and roundoff have their own count", {
  new <- f$new
  new$matrices$S$values[1, 1] <- new$matrices$S$values[1, 1] + 5e-11
  x <- compare_fixture(f$old, new)
  stopifnot(x$exact == 0L, x$approved_renames == 1L, x$roundoff_equivalent == 1L,
            x$roundoff_only == 0L, x$renames_only == 0L, x$renames_and_roundoff == 1L)
})
for (field in c("lbound", "ubound")) {
  pass(paste("tiny roundoff is allowed in", field), {
    old <- f$old
    old$matrices$S[[field]] <- matrix(10, 2L, 1L)
    new <- old
    new$matrices$S[[field]][1, 1] <- 10 + 5e-11
    stopifnot(compare_fixture(old, new)$roundoff_only == 1L)
  })
}
pass("tiny checksum roundoff is allowed", {
  old <- f$old
  old$data <- list(checksum = 1000)
  new <- old
  new$data$checksum <- 1000 + 5e-11
  stopifnot(compare_fixture(old, new)$roundoff_only == 1L)
})
pass("large checksum roundoff uses the documented relative scale", {
  old <- f$old
  old$data <- list(checksum = 1e9)
  new <- old
  new$data$checksum <- 1e9 + .05
  x <- compare_fixture(old, new)
  stopifnot(x$roundoff_only == 1L, x$max_abs_roundoff > 1e-10,
            x$max_scaled_roundoff > 0, x$max_scaled_roundoff <= 1e-10)
})

old <- f$old
old$data <- list(checksum = 1e9)
bad <- old
bad$data$checksum <- 1e9 + .2
reject("large checksum changes beyond the relative bound fail", old, bad,
       "Unexpected structural differences")

bad <- f$old
bad$matrices$S$values[1, 1] <- 1 + 1e-8
reject("meaningful numerical changes still fail", f$old, bad,
       "Unexpected structural differences")

bad <- f$old
storage.mode(bad$matrices$S$values) <- "integer"
reject("equal numbers with different storage types still fail", f$old, bad,
       "Unexpected structural differences")

bad <- f$old
dim(bad$matrices$S$values) <- c(1L, 2L)
reject("numeric dimensions remain exact", f$old, bad,
       "Unexpected structural differences")

bad <- f$old
dimnames(bad$matrices$S$values) <- list(c("x", "y"), "value")
reject("numeric dimension names remain exact", f$old, bad,
       "Unexpected structural differences")

old <- f$old
old$fitfunction <- list(offset = 1)
bad <- old
bad$fitfunction$offset <- 1 + 5e-11
reject("numeric fields outside the allowlist remain exact", old, bad,
       "Unexpected structural differences")

bad <- f$old
bad$matrices$S$lbound[1, 1] <- NaN
reject("NA cannot silently become NaN", f$old, bad, "Unexpected structural differences")

old <- f$old
old$matrices$S$ubound[1, 1] <- Inf
bad <- old
bad$matrices$S$ubound[1, 1] <- -Inf
reject("infinite bound signs remain exact", old, bad, "Unexpected structural differences")

bad <- f$new
bad$matrices$AADD1$values[1, 1] <- 1 + 5e-11
reject("ADD1 constants retain exact content checks", f$old, bad,
       "ADD1 constants differ in content")

pass("zero tolerance can require exact numeric agreement", {
  new <- f$old
  new$matrices$S$values[1, 1] <- 1 + 5e-11
  message <- tryCatch({compare_fixture(f$old, new, numeric_tolerance = 0); ""},
                      error = function(e) conditionMessage(e))
  stopifnot(grepl("Unexpected structural differences", message, fixed = TRUE))
})
pass("numeric allowance cannot exceed the documented bound", {
  message <- tryCatch({compare_fixture(f$old, f$old, numeric_tolerance = 1e-6); ""},
                      error = function(e) conditionMessage(e))
  stopifnot(grepl("numeric_tolerance must be between 0 and 1e-10", message, fixed = TRUE))
})
# Build a small structural fixture carrying the same reviewed S/M/R contract
# as the immutable longitudinal mediation cases, without loading/replacing them.
med_fixture <- function(prefix = "") {
  vars <- c("S1", "M1", "R1")
  values <- matrix(c(1, .6, .9, .6, 4, .9, .9, .9, 9), 3L, dimnames = list(vars, vars))
  labels <- matrix(paste0(prefix, c("S_residuals", "SM_RES", "SR_RES", "SM_RES",
                                    "M_residuals", "MR_RES", "SR_RES", "MR_RES", "R_residuals")), 3L)
  old <- list(matrices = list(S = list(values = values, labels = labels,
                                       free = matrix(TRUE, 3L, 3L))), submodels = list())
  new <- old
  new$matrices$S$values[cbind(c(2L, 3L), c(3L, 2L))] <- 1.8
  list(old = old, new = new)
}
compare_med <- function(old, new, id = "MED_LIN_longitudinal") {
  invisible(capture.output(result <- compare_model_structures(setNames(list(old), id),
    setNames(list(new), id), expected_count = 1L)))
  result
}
reject_med <- function(label, old, new, id = "MED_LIN_longitudinal") {
  message <- tryCatch({compare_med(old, new, id); NA_character_}, error = conditionMessage)
  if (is.na(message)) stop(label, ": unexpected acceptance")
  pass(label, TRUE)
}
g <- med_fixture()
pass("MED M-X start correction is derived from unchanged residual variances", {
  x <- compare_med(g$old, g$new)
  stopifnot(x$exact == 0L, x$approved_mediation_starts == 1L,
            x$mediation_starts_only == 1L, x$approved_renames == 0L)
})
pass("MED start correction recurses through class models", {
  a <- med_fixture("c1"); b <- med_fixture("c2")
  old <- list(submodels = list(Class1 = a$old, Class2 = b$old))
  new <- list(submodels = list(Class1 = a$new, Class2 = b$new))
  for (container in c("MGroup", "MIX", "MIXc")) {
    x <- compare_med(old, new, paste0(container, "_MED_BLS_longitudinal"))
    stopifnot(x$approved_mediation_starts == 1L, length(x$mediation_start_scopes) == 2L)
  }
})
reject_med("MED correction is restricted to approved model identifiers", g$old, g$new,
           "MED_LIN_baseline")
for (bad_value in c(.9 + 1e-6, 1.8 + 1e-6, -1.8, NA_real_)) {
  bad <- g$new
  bad$matrices$S$values[2, 3] <- bad_value
  reject_med("MED wrong/asymmetric covariance correction is rejected", g$old, bad)
}
bad <- g$new
bad$matrices$S$free[2, 3] <- FALSE
reject_med("MED correction does not permit changed constraints", g$old, bad)
bad <- g$new
bad$matrices$S$labels[2, 3] <- "different"
reject_med("MED correction does not permit changed labels", g$old, bad)
bad <- g$new
bad$matrices$S$values[2, 2] <- 4.1
reject_med("MED correction does not permit changed residual variance", g$old, bad)
bad <- g$new
bad$matrices$S$values[1, 2] <- .61
reject_med("MED correction does not permit another changed covariance", g$old, bad)
bad <- g$new
bad$matrices$S$values[2, 2] <- 9
bad$matrices$S$values[cbind(c(2L, 3L), c(3L, 2L))] <- 2.7
reject_med("MED changing variance and covariance consistently still fails", g$old, bad)
cat("All", checks, "structural comparator self-tests passed.\n")
