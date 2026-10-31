# nlpsem 1.0.0.9005 validation

This candidate restores per-person stats::nls initialization for EXP, JB and BLS
using the 0.3 formulas and default Gauss-Newton controls. The inspected 0.4 and
earlier 1.0 development sources actually used grid/OLS initialization. The 9004
correctness fixes are retained. This remains a development candidate.

## Matched numerical study

Both platforms fit exactly the same serialized starting models and data using
CSOLNP, one computation thread, seed 20191029 and identical retry budgets.
All 14 fits completed without timeout. macOS used R 4.5.3; Linux ARM64 used
R 4.5.2; both used OpenMx 2.22.11.

The fitted values are compared using numerical tolerances, not exact equality:
1e-6 for scalar -2 log likelihood and 1e-4 for the parameter vector. The latter
is an aggregate relative/absolute comparison, not a per-parameter maximum.
The standalone comparison now explicitly uses testthat edition 3; postanalysis
under editions 2 and 3 produced the same decisions for every completed case.

| Case | Basic diagnostic screen, both platforms | Paired point-estimate comparison |
|---|---|---|
| LGCM_JB_r | Fail | Fail |
| LCSM_JB_f | Fail | Fail; likelihood alone agrees |
| LCSM_JB_r | Fail | Fail; likelihood alone agrees |
| TVC_LGCM_BLS_r_d1 | Pass | Pass |
| MGM_LGCM_BLS_r | Pass after correcting diagnostic cross-block handling | Pass |
| MED_BLS_longitudinal | Pass | Pass |
| LGCM_EXP_r | Fail | Fail |

The diagnostic checker previously treated an asymmetric cross-process covariance
block as a complete covariance matrix. It now recognizes such blocks from RAM
variable mappings, reports them separately, and checks the assembled joint
covariance. Marginal/joint covariance, Hessian and negative-variance checks are
retained; no thresholds were relaxed. Saved fits were reassessed without refitting.

The screen checks convergence, finite estimates, stored covariance admissibility,
a positive-definite stored Hessian and finite positive SEs. It does not establish
identification, global optimality, interval coverage or SE reproducibility.

MED standard errors remain sensitive across platforms. Y_alpha1 SE is 1.182703
versus 1.561111. The naive 95% Wald interval for M_alpha2 includes zero on Linux
but excludes it on macOS. This remains a qualification concern despite close
point estimates and passing the basic diagnostic screen.

NLS failures remain explicit. For 200 math trajectories, NLS converged for
EXP 191/200, JB 49/200 and BLS 135/200; other rows used audited heuristic
fallbacks. Missing/reordered occasions and the historical BLS seed-indexing
error are handled explicitly. No grid values are silently substituted.

## Validation and remaining gates

The macOS regular suite passed 3,955 assertions before adding the 49 new NLS
assertions; the Linux regular suite passed all 4,004 then-current assertions.
The corrected diagnostic suite passes 86 assertions, and final combined
diagnostic/NLS tests pass 135 on Linux. The 688-model structural replay has
zero unexplained differences when supplied the frozen historical starts.
The basic seven-family synthetic smoke and limited 50-replicate linear recovery
screen also pass; they are not nonlinear coverage qualification.

The final macOS archive check has zero errors, zero warnings and two notes
(development version and local HTML Tidy); 4,011 CRAN-mode assertions pass and
11 tests intentionally skip. Examples, donttest, vignettes and PDF manual pass.
The final Linux archive check has zero errors/warnings/notes, with tests/manual
omitted and tests run separately. Artifact hashes are recorded in the sibling
development validation folder and NLS_VALIDATION.md. Earlier check attempts,
source snapshots, original fit classifications and reassessments are preserved.

V1.0.0 remains on hold pending the JB/EXP failures, mediation uncertainty
reproducibility, broader recovery/coverage qualification, the five historical
reference gaps, refreshed tutorial fits and remaining platform/runtime checks.
No reference values were replaced, and nothing was pushed or published.
