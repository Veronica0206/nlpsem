# Validation of nlpsem 1.0.0.9003

This development update addresses the 11 findings from the independent review of
Phase C commit `31f3d61bb256b12b9ec1b345dc364c8680925fe6`. It retains the preceding
Phase B refactor, correctness work, and Phase C specification/inspection APIs.
The manuscript, release artifacts, and historical fitted examples are unchanged.

## Corrections and regression coverage

| Finding | Change | Regression evidence |
|---|---|---|
| F1 | Mediation uses outcome/mediator/predictor time prefixes in that order. | Distinct time values and unequal schedules; LIN/BLS, both predictor types, all containers. |
| F2 | Estimates, SEs, and numeric confidence bounds retain full precision. | All seven output entry points; unit invariance, small SEs, exact OpenMx comparisons. |
| F3 | TVC constraints use an explicit structural mask. | Zero/nonzero start invariance across decomposition methods, TIC counts, outcome models, and containers; structural-zero controls. |
| F4 | EXP multivariate starts retain every outcome's free growth factors. | Actual covariance blocks for two, three, and four outcomes and both PD start patterns across LGCM/LCSM and containers. |
| F5 | QUAD/JB constants have outcome-specific names. | Actual rate and change algebras with unequal waves; independent derivative/integration expectations. |
| F6 | Residual covariance paths match common wave IDs. | Identical, reordered, partial-overlap, and disjoint schedules. Disjoint pairs receive a clear error before optimization. |
| F7 | Kappa variance includes uncertainty in empirical marginals. | Independent numerical derivatives for imbalanced, multi-category, negative-agreement, and sparse tables at three confidence levels. |
| F8 | Summary parsing consumes the complete class/group number. | Distinct residual values for 9, 10, and 11 groups/classes. |
| F9 | One-class output omits nonexistent free mixture weights. | Real one-/two-class public fits with parameter output, plus weight-output checks. |
| F10 | The LRT example uses the correct full/reduced fixture order. | The actual documented example asserts df=4 and a finite p-value. |
| F11 | The declared testthat minimum is 3.2.0. | Regular suite run against an isolated 3.2.0 installation and a dedicated hosted job. |

The four new regression files contain **1,437 expectations**, included in the
regular-suite total below. Their tests reproduce failures against the unchanged
Phase C source. The corrected LRT example was separately executed against both
versions: Phase C gives df=-4 and no p-value; this version gives df=4 and
p=6.401156e-06. Historical output-table tests compare display rounding only when
reading old saved fixtures; new outputs are compared with full-precision values.

## Local results

Validated on October 24, 2026 using R 4.5.3, OpenMx 2.22.11, and macOS arm64.

| Check | Result |
|---|---|
| Regular suite, testthat 3.3.2 | 2,990 passed; zero failures, errors, warnings, or skips |
| Regular suite, declared minimum testthat 3.2.0 | 2,990 passed; zero failures, errors, warnings, or skips |
| Full numerical reference suite | 36 eligible fits matched unchanged references; 147 passed, five explicit historical skips; zero failures, errors, or warnings |
| Full source-package check with PDF manual | 0 errors, 0 warnings, 2 notes (incoming checks and local HTML validator) |
| Historical structural grid | 688 constructed; 676 exact matches, 12 approved ADD1 renames; no unexplained differences or numeric drift locally |
| Comparator self-tests | 37 passed, including deliberately missing/misassigned constants and unexpected structural changes |
| Documentation | 48 HTML pages; 1,339 local targets checked, none missing; 27 article plots have descriptions |

The pre-push incoming note includes the development version and a 404 for this new
validation-record URL before publication. The other note is an outdated local HTML
Tidy; the PDF manual builds successfully. Hosted checks verify the published head.

The minimum-dependency run clears the inherited LC_ALL environment variable:
testthat 3.2.0 otherwise emits locale-change warnings under the local R 4.5
environment. No package warnings were suppressed.
The hosted minimum-dependency job uses R 4.5.3 because testthat 3.2.0 cannot
load with R 4.6.1's removed C APIs. Current-R coverage uses modern testthat in
the main package and full-fitting jobs; the historical minimum is tested as
an explicit R 4.5.3 combination.

Numeric tables and interval endpoints now carry more digits. Existing serialized
`@Estimates` tables retain their historical rounding; refitting or recomputing the
table is necessary to recover precision. The existing character p-value format
and S4 slots are preserved.

Every multivariate outcome pair must share at least one record ID because the
model estimates a residual covariance for each pair. Disjoint schedules are
rejected clearly. Kappa intervals treat supplied paired labels as observed and
subjects as independent; they do not include uncertainty from latent-model
estimation or class assignment.

## Reference limits and reproducibility

The five unchanged unavailable Phase B references remain explicitly excluded:
`LGCM_JB_f`, `TVC_LGCM_BLS_f_d1_TIC`, `MIX_LGCM_BLS_f`,
`MIX_LGCM_BLS_f_clusterTIC`, and `MIX_MED_LIN`. No fixture was regenerated
or promoted to make this update pass. Numerical fitting tolerances remain
1e-6 for likelihoods and 1e-4 for free-parameter comparisons.

The historical structural snapshot, hashes, source revision, environment, and
comparison policy are committed under `data-raw/structural-baseline/`.
The hosted structural gate replays the exact saved Phase C data, arguments, and
starting values through the public fitting wrappers. Capture is accepted only
when all 688 old-source structures exactly match the immutable historical
snapshot. This avoids conflating platform-dependent automatic initialization
with changes in model construction. The full numerical fitting suite continues
to exercise automatic initialization.
Only the reviewed ADD1 renames are normalized. A bounded 1e-10 scaled numeric
allowance applies to matrix values/bounds and data checksums for cross-platform
roundoff; types, dimensions, names, free masks, labels, formulas, and ADD1
contents remain exact. Statistical correctness is also checked through independent
numerical expectations and invariants in the new tests.

From the repository root:

```sh
mkdir -p test-results
Rscript -e 'source("data-raw/run_validation_tests.R")' . . test-results/full-tests
Rscript -e 'source("data-raw/snapshot_models.R")' . test-results/model-structures.rds . replay data-raw/structural-baseline/phase-c-1.0.0.9002-inputs.rds
Rscript data-raw/test_structure_comparator.R
Rscript data-raw/compare_model_structures.R data-raw/structural-baseline/phase-c-1.0.0.9002.rds test-results/model-structures.rds
R CMD build .
NOT_CRAN=false R CMD check --as-cran nlpsem_1.0.0.9003.tar.gz
```

The [full fitting workflow](https://github.com/Veronica0206/nlpsem/blob/ca795d15c5daba5c2ecfcd67e4cc1d4bae21561b/.github/workflows/full-fitting-tests.yaml) now runs on pull
requests and checks the PR head, retains numerical/structural evidence, and tests
the minimum dependency. Three-platform package checks and the documentation
workflow run separately. Hosted results must be read for the final commit on
[PR 1](https://github.com/Veronica0206/nlpsem/pull/1); local results do not imply
hosted success. This is a development update, with no merge or release publication.
