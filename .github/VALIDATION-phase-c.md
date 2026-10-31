# nlpsem 1.0.0.9002 — local validation record

Validated on 2026-10-24. This development update includes the preceding unified
builder refactor and correctness fixes, followed by Phase C's fitted-model API,
provenance, documentation, and CI infrastructure. The published manuscript and
existing releases are not modified by this branch.

## Results

| Check | Local result |
|---|---|
| Regular regression suite (`NOT_CRAN=true`) | 1,553 passing expectations; no failures, errors, warnings, or skips |
| Numerical reference suite (`NOT_CRAN=true`) | All 36 eligible fitted models reproduced their stored references; 147 passing expectations; five explicit historical-reference exclusions |
| Combined tests | 1,700 passing expectations; no failures, errors, or warnings |
| New Phase C coverage | 447 expectations for provenance, standard methods, serialization, legacy fits, lazy evaluation, RNG preservation, and canonical curve labels |
| Structural grid | 688 models constructed without errors; canonical structures exactly match the 1.0.0.9001 baseline |
| Full `R CMD check --as-cran` | 0 errors, 0 warnings, 2 notes |
| Installed-package example | Passed; fitted model status 0; automatic plot verified |
| Additional automatic-plot checks | Five newly fitted families all status 0; eight figures match explicit `getFigure()` calls |
| pkgdown | 47 HTML pages; 1,292 local link/asset targets exist; all 27 article plots have descriptions |

The package check included ordinary and `donttest` examples, tests, rebuilt
vignettes, and the PDF manual. Its two notes concerned the development version
number and the local system's outdated HTML Tidy. HTML manual validation is not
claimed. The pkgdown build also reported an unset hosted URL and Pandoc
syntax-highlighting deprecation messages; no website was deployed.

Environment: Apple Silicon, macOS 27.0.1, R 4.5.3, OpenMx 2.22.11,
testthat 3.3.2, roxygen2 7.3.3, pkgdown 2.2.0, and rlang 1.1.7. Numerical
regressions used one OpenMx thread, the existing fixed seeds/optimizer, and
unchanged comparison tolerances (1e-6 likelihood; 1e-4 free parameters).

## Known numerical-reference limitations

The five excluded references were already nonconverged. They remain explicit
exclusions, and a coverage test prevents that set from silently growing. Final
diagnostic reruns produced the following statuses without changing references:

| Case | Status |
|---|---:|
| LGCM_JB_f | 6 |
| TVC_LGCM_BLS_f_d1_TIC | 0 |
| MIX_LGCM_BLS_f | 6 |
| MIX_LGCM_BLS_f_clusterTIC | 6 |
| MIX_MED_LIN | 5 |

Four remain nonconverged. A converged diagnostic rerun does not make an unavailable
historical reference valid. These local results do not establish convergence or
statistical recovery for every supported model.

## Reproduce

From the package repository root, with the declared dependencies and `pkgload`
installed:

```r
Sys.setenv(NOT_CRAN = "true")
OpenMx::mxOption(NULL, "Number of Threads", 1L)
testthat::test_local(".", reporter = "summary", stop_on_failure = TRUE)
```

The broad construction-only grid is available separately:

```sh
Rscript data-raw/snapshot_models.R . /path/to/model-structures.rds
```

Build the package, then run `R CMD check --as-cran` on the resulting source
archive. Build documentation with `pkgdown::build_site()` when pkgdown is installed.

The published code, tests, data, and installed documentation match the validated
Phase C source. Repository preparation preserves the newer public README's
portfolio context, adds Git ignores and this validation note, and removes trailing
whitespace in a historical data-preparation script. It does not alter fitting code
or reference fixtures.

## Hosted checks

Local validation is distinct from GitHub Actions results. This branch adds routine
checks on Windows, macOS, and Linux, a manual full-fitting workflow, and a
website-artifact workflow. Their live status is shown on the pull request; no
hosted success is implied by the local results above. The workflows do not publish
releases or deploy documentation.
