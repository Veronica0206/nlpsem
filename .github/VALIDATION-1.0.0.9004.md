# nlpsem 1.0.0.9004 validation

This development update addresses the October 25, 2026 publication-readiness
review, starting from 1.0.0.9003 commit
ca795d15c5daba5c2ecfcd67e4cc1d4bae21561b. **V1.0.0 publication remains on hold.**
The seven reproduced code defects are repaired; difficult nonlinear fits still
fail numerical qualification.

## Fixes

Mediation paths match explicit wave IDs and use the correct M-X covariance start.
Posterior output supports one class, validates class counts/covariates, and uses
stable log-space normalization. Summaries align residual schemas and retain one
proportion column per class. Custom parameter names are validated. Kappa displays
the requested confidence level. Factor-score support and invisible-return
contracts are explicit. Opt-in getDiagnostics checks stored covariance, Hessian,
SEs and component sizes without refitting or claiming statistical acceptance.

Focused regressions fail on unchanged 9003 and pass on 9004. Reference files,
historical serialized examples and generic numerical tolerances are unchanged.
The structural comparator allows only mathematically checked starting-value
corrections in addition to its previously approved constant renames.

## Local evidence

| Check | Result and scope |
|---|---|
| macOS regular suite, R 4.5.3 / OpenMx 2.22.11 | 3,955 pass; zero failures, errors, warnings or skips |
| Minimum testthat 3.2.0 on macOS R 4.5.3 | 3,955 pass; zero failures, errors, warnings or skips; not minimum R/OpenMx qualification |
| Linux ARM64 regular suite, R 4.5.2 / OpenMx 2.22.11 | 3,955 pass; zero failures, errors, warnings or skips; Ubuntu 24.04.4 container |
| Final macOS archive check --as-cran | 0 errors, 0 warnings, 2 notes; examples, donttest, tests, vignettes and PDF manual pass |
| Final Linux package check | 0 errors, 0 warnings, 0 notes; --no-tests --no-manual, with regular tests run separately |
| Structural replay | 688 constructions: 668 exact, 12 approved ADD1 renames, 8 MED start corrections; 0 unexplained differences; 49 comparator tests pass |
| Independent synthetic smoke, macOS and Linux | 7/7 basic family cases accepted on each; not all nonlinear curves/containers |
| Linear mean recovery, macOS and Linux | 50/50 replicates accepted on each; intercept/slope coverage 47/50 and 46/50 |
| Full macOS reference suite | 146 passing expectations, 1 failure, 5 historical exclusions; MED_BLS_longitudinal status 5, with reference-like likelihood and estimates |
| Six formerly divergent Linux cases | 1/6 matches unchanged references; 5/6 fail; no timeouts |
| Historical-object API smoke | 23 diagnostic reports and 46 summaries pass; not requalification of the fits |
| Documentation | 51 HTML pages; 1,420 local file targets checked, zero broken targets; local build only |

The two macOS notes concern the development version number and the local HTML
Tidy installation. The 50-replicate linear screen has nominal coverage Monte
Carlo SE about 0.031: it does not establish precise coverage or qualify nonlinear
or mixture inference. All claims are bounded by the recorded platforms and
versions. Prepared CI is not a completed hosted run.

## Remaining V1.0.0 gates

1. Resolve MED_BLS_longitudinal on macOS: the corrected start exposes an indefinite
   Hessian and problematic SEs despite reference-like likelihood and estimates.
2. Resolve Linux mismatches for LGCM_JB_r, LCSM_JB_r, TVC_LGCM_BLS_r_d1,
   MGM_LGCM_BLS_r and MED_BLS_longitudinal. LCSM_JB_f now matches. Use matched
   environments, repeated starts, admissibility and identification diagnostics;
   do not widen tolerances to erase differences.
3. Add accepted cases or explicitly limit support for the five original exclusions:
   LGCM_JB_f, TVC_LGCM_BLS_f_d1_TIC, MIX_LGCM_BLS_f,
   MIX_LGCM_BLS_f_clusterTIC and MIX_MED_LIN.
4. Extend recovery/uncertainty calibration to nonlinear, weak-mixture, missingness,
   scaling, irregular-schedule and boundary scenarios. The simple nested-normal
   bootstrap smoke establishes execution only; random/fixed-knot LRT calibration
   remains unresolved.
5. Refit and reconcile the 23 tutorial fits in a separate release area. The
   historical examples and inadmissible TVC example are now identified explicitly,
   and were not silently constrained or replaced.
6. Run Windows, R-devel, R-oldrel and a defensible minimum R/OpenMx combination on
   the final candidate. Complete reverse-dependency checks and final 0.4-to-1.0
   migration/release prose.
7. Resolve the original data source release/export identifier. The exact local
   exports, input hashes, named variables, seed and reconstruction are verified:
   all 500 rows and 49 columns match. This is a complete-case teaching sample.
8. Build/check the actual 1.0.0 archive, then perform separately authorized
   integration, tagging, release and CRAN submission.

All 4,193 pre-existing project files match their before-work hashes. Detailed
logs, probes and checksums are in the versioned development folder alongside this
checkout. No hosted CI, merge, release or CRAN submission was performed for 9004.
