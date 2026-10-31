# nlpsem 1.0.0 validation

This review extends candidate de5769206a61fbbdfbd0daba19d755b1dceab642. It corrects
mediation disturbance covariance starts, fitted change-score plotting and selected
waves, and simulator truth metadata for shifted time origins. The per-person NLS
method, SEM optimizer and retry controls, numerical reference values and comparison
tolerances are preserved.

The follow-up to merged revision 8519d00 repairs intrinsic EXP/JB LGCM plotting
when covariates shift the mean rate, and restores explicit mediation `res_scale`
precedence. Plot regressions compare the mean and numerical derivatives against
the RAM expectation, including class-specific and TVC components, expansion-point
amplitudes, and parent-model parameter covariance. Mediation regressions check
per-process and per-class scales for linear and bilinear-spline starts, plus the
unchanged `NULL` automatic estimate and fallback. Specified covariance matrices
in derivative tests check propagation, not empirical interval coverage.

The follow-up to merged revision 548cd83 extends RAM-referenced plotting tests
to intrinsic bilinear splines. Curves and intervals retain the expansion knot,
smoothed loadings, adjusted internal means and mean knot deviation. Tests cover
simultaneous slope and knot effects, knot-only effects, class-local evaluation
with parent covariance, TVC decompositions and multivariate processes. Stored
BLS plot tests now use RAM expectations and numerical derivatives rather than
the reconstructed sharp-spline formula. Nonintrinsic spline controls remain.
The NLS method, optimizer, model equations, references and their tolerances
are unchanged.

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

Targeted review on R 4.5.3/OpenMx 2.22.11 includes 104 mediation-start/NLS
assertions, 387 plotting assertions, and 85 new simulator-origin assertions.
All 23 shipped example fits rendered 34 figures without plotting warnings;
eight additional selected-wave and coarse nonlinear figures were inspected.
Sixty-two deterministic simulator comparisons retained exactly the generated
data, original truth, growth factors, settings and RNG behavior. The added
truth attributes are metadata, not changes to random generation.

The 688-model construction replay has 596 exact matches, 12 approved constant
renames, eight approved historical mediation-start corrections and 72 approved
derived LCSM mean corrections, with no unexplained changes locally. The
comparator's 76 self-tests include rejection of incorrect or partial formula
changes. Replay supplies historical starts; separate tests and fitted references
exercise the corrected automatic starts.

The maintainer's submission bundle retains the full local test/check logs,
review findings, archive manifest and final hosted run identifiers.

## Numerical scope

The reference inventory has 43 configurations: 33 positive comparisons, five
explicitly nonconverged exclusions and five diagnostic-only inadmissible fixtures.
Positive comparisons require convergence, agreement with the finite reference
likelihood, no reported covariance issues, a positive-definite Hessian and
complete finite positive SEs. The independent numerical smoke workflow additionally
requires available covariance checks.

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
