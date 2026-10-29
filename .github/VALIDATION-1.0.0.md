# nlpsem 1.0.0 validation

The release candidate starts from development commit e9a631e5f3eb6c8067107f9c358b667b0eb9ad87.
The final review corrects simulator argument validation and plotting edge cases, qualifies
single-sample approximation results, and prepares release metadata. NLS initialization,
the SEM optimizer and retry settings, reference values and numerical tolerances are unchanged.

## Checks and retained evidence

Use the workflow results for the commit being assessed. A configured workflow does not
establish that a run passed.

- The submission-archive workflow builds with current R release, checks that archive, then
  checks the identical archive on R-devel, Windows and macOS. It retains the archive,
  SHA-256, source revision, build environment and check logs.
- The package workflow separately checks release R on Linux, Windows and macOS,
  R-devel and R-oldrel-1.
- The full-fitting workflow runs the 688-model structural replay, all numerical
  references on macOS, and the regular suite with testthat 3.2.0 on R 4.5.3.
- The numerical workflow runs seven independently generated basic model-family
  cases and a 50-replicate linear LGCM mean-recovery screen on three platforms.
- The documentation workflow builds the reference pages and all articles.

Targeted local review on R 4.5.3/OpenMx 2.22.11 passed 234 simulator assertions,
108 plotting assertions, and 219 diagnostics/standardization assertions. All 23
shipped example fits rendered 34 figures without plotting warnings. Additional
checks reproduced rejection of missing SEs and a warning for an inadmissible
mixture solution despite optimizer status zero.

The maintainer's submission bundle retains the full local test/check logs,
review findings, archive manifest and final hosted run identifiers.

## Numerical scope

The reference inventory has 43 configurations: 33 positive comparisons, five
explicitly nonconverged exclusions and five diagnostic-only inadmissible fixtures.
Positive comparisons require convergence, finite likelihood, covariance
availability, a positive-definite Hessian and complete finite positive SEs.

Reference agreement is regression evidence. It does not prove global optimality,
identification, or calibrated uncertainty across all data-generating settings.
The portable recovery experiment concerns linear LGCM means only. The
approximation examples in data-raw/approximation-bias retain their one-sample
design and do not estimate repeated-sample bias or coverage.

Four stored group/mixture example solutions fail covariance admissibility;
the guides identify them. They remain illustrations of diagnostics, not accepted
scientific solutions. Missing uncertainty information is distinct from a check
that was performed and passed. Users should inspect diagnostics and compare
defensible alternative starts before interpreting model results.

## Release boundary

Building or checking the archive does not publish the package. The maintainer
submits the prepared source tarball through the CRAN submission form and confirms
the resulting email. GitHub main, tags, releases and documentation deployment are
separate from preparation of the submission candidate.

