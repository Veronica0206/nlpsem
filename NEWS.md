# nlpsem 1.0.0.9006 (standardized fitting development update)

* New `getSimData()` simulates wide-format longitudinal data from a latent growth curve model (linear,
  quadratic, negative exponential, Jenss-Bayley, bilinear spline) or a latent change score model
  (quadratic, negative exponential, Jenss-Bayley, nonparametric), with fixed or random growth ratios and
  knots, individually varying measurement occasions, growth TICs, and several latent classes with fixed
  proportions or cluster-TIC multinomial logits. The true parameter values are attached under the names of
  the matching fit's `Estimates`, and each person's growth factors are attached too, so a fit can be
  compared with the values that generated the data. It also generates multivariate growth models (two or
  more outcomes, with growth factors and residuals correlated across outcomes), models with a time-varying
  covariate (a direct effect, or decompositions 1, 2 or 3, in which the outcome responds to the TVC's slopes,
  interval changes or change from baseline and the TVC's baseline predicts the outcome's growth factors),
  and longitudinal mediation models (linear or bilinear spline, with a baseline or longitudinal predictor).
  A new section of the workflow vignette generates data, fits the matching model and compares the
  estimates with the truth.
* Quadratic, negative exponential and Jenss-Bayley models, with fixed or random ratios, are
  fitted to standardized data and converted back exactly in `getLGCM()`, `getLCSM()`,
  `getTVCmodel()`, `getMGM()`, `getMGroup()` and `getMIX()`. Each outcome is centred and scaled
  at its first recorded occasion, a TVC is scaled, and time is divided by the mean of its last
  recorded occasion. The model starts from the raw-scale starting values, mapped exactly, and the
  optimizer and retry settings are unchanged.
* The estimates, Hessian, parameter covariance matrix, standard errors and -2 log-likelihood are
  converted back exactly. The unit factor of every free parameter was checked against the
  likelihood on all 480 QUAD/EXP/JB specifications of the structural grid; the largest relative
  mismatch was 1.5e-15. `getSpecification()` records the centres, scales and condition numbers
  under `fitting$standardization`. `options(nlpsem.standardize = FALSE)` fits on the raw scale.
* On the 200-person reference sample with CSOLNP on macOS and Linux, all six QUAD and seven of
  eight EXP cases pass the basic diagnostic screen on both platforms, with -2LL agreeing within 4e-9.
* The fitting-layer reference values (`tests/testthat/golden/fit_golden.rds`) were regenerated with
  1.0.0.9006 on macOS. The QUAD and EXP references are unchanged. The four JB references fit the
  nearly linear math scores and stop with status 6; they are excluded explicitly, and three
  bilinear spline references that did not converge before now do.
  The random-ratio EXP fits reach the 1.0.0.9000 reference likelihoods that the 9005 raw-scale
  fits missed. JB fits improve but stop with status 6: their random slope and acceleration terms
  are nearly collinear on this sample.
* Mediation, linear, nonparametric and bilinear spline models are fitted as before.
* In multiple-group and mixture models, the standardized fit scales each outcome (and a TVC) by the
  standard deviation of its first occasion pooled within the groups, or within the provisional classes
  of the automatic starting values; the overall standard deviation is used with user-supplied starts.
  With well-separated groups or classes the overall standard deviation is several times the
  within-group one, so every within-group variance became small on the standardized scale and the
  optimizer stopped with status 6 short of the optimum that the raw-scale fit reached with status 0.
  `getSpecification()` records the scale basis.
* A growth or class-membership covariate that is also a time, outcome or TVC column (for example
  baseline age, `growth_TIC = "T1"`) is read from an unscaled copy in the standardized fit. It was
  rescaled with the time column before, which made the converted solution and likelihood wrong.
* Mixture starting values: each provisional class must have enough persons for its class-specific
  starts (one more than the number of growth factors, and two more than the number of growth TICs).
  When k-means on the individual growth factors leaves a smaller class (it can isolate a few persons
  with outlying negative exponential fits), the classes are formed on the normal scores of the growth
  factors instead; if a class is still too small, the k-means classes are kept with a warning. A covariate
  effect that the persons used cannot determine (growth-factor effects of TICs, mediation paths and the
  TVC effect) starts at 0 with a warning. Before, a two-TIC negative exponential mixture on 200 persons
  failed with "infinite or missing values in 'x'".
* Nonparametric (latent change score) starting values: each person's residuals now come from that
  person's own intercept and scale along the average interval-specific slopes. Measured against the
  average curve itself, they also contained the differences between persons, which made the residual
  variance start about ten times too large (twenty in mixtures, which double it); two-class
  nonparametric mixtures then collapsed into one class. The change applies to every model family.
* The fitting-layer reference tests check convergence, covariance admissibility, the Hessian and the
  standard errors with `getDiagnostics()` before comparing numbers. Five references reach status 0
  with inadmissible solutions on the 200-person subsample (four two-class mixtures and the
  multiple-group bilinear spline model); they are kept as fixtures for `getDiagnostics()` rather than
  as model coverage, which leaves 33 positive references. Two fixed-ratio JB references fitted to the
  reading scores of all 500 children were added.
* Bilinear spline starting values: the mean knot now starts at the common knot profiled over all persons
  (each with their own intercept and slopes; the least-squares knot of the fixed-knot model) instead of
  the average of the per-person knots. Persons whose own fit fails keep the midpoint seed, which pulled
  that average towards the middle of the time range (for example 3.9 for a true knot of 3.5). The
  per-person knots are shifted to the profiled knot, so the knot variance rule is unchanged. From the
  new start, two random-knot references reach better solutions, and their stored values were updated:
  `TVC_LGCM_BLS_r_d1` (-2LL 16093.60 to 16090.76) and `MGM_LGCM_BLS_r` (25005.40 to 24975.52), both
  passing `getDiagnostics()`. The other references are unchanged.
* Plots of fits whose `records` are not in ascending order no longer fail. Plots of TVC models label
  their curves as the growth component, which excludes the effect of the TVC, and add the full
  model-implied mean of the outcome (or its change from baseline) at each wave, for each class or
  group in mixture and multiple-group plots.

# nlpsem 1.0.0.9005 (NLS initialization development update)

* Restores per-person nonlinear least-squares initialization with `stats::nls()`
  for EXP, JB and BLS trajectories, following the 0.3 formulas. The shared
  initializer supplies the same method across the seven model families.
* Retains the current residual-data interface and records individual NLS failures
  and heuristic fallbacks. A fallback is not reported as NLS convergence.
* Corrects the historical BLS second-slope seed indexing and handles missing or
  reordered observations explicitly. The final SEM optimizer is unchanged.
* Matched numerical comparisons use CSOLNP on both platforms and identical
  serialized unfitted models. Earlier grid-initialized reference results remain
  historical evidence and are not replaced automatically.
* Fit diagnostics now distinguish cross-process covariance blocks from complete
  covariance matrices. Cross blocks are reported separately and their joint RAM
  covariance is checked; legitimate asymmetry is no longer flagged as invalid.

# nlpsem 1.0.0.9004 (development review fixes)

* Mediation residual covariance paths now match explicit wave IDs, including
  reordered schedules and outcome prefixes containing digits. Longitudinal M-X
  covariance starts use the M and X residual variances. Disjoint process schedules
  are rejected clearly because the current mediation output requires a covariance
  for every process pair.
* One-class mixture posterior calculations return numeric missing entropy; their
  summaries and plots work. Posterior normalization uses log space, validates the
  fitted class count, and rejects invalid likelihood rows or mismatched covariates.
* Multivariate summaries have one proportion column per class, combine residual
  schemas by name, and retain absent values as NA. Displaying models together does
  not establish that their information criteria are statistically comparable.
* Parameter names must have valid, unique entries and the correct length instead
  of being silently recycled or padded. The documented LCSM core-name convention
  is retained.
* Kappa output retains and prints the requested confidence level. Legacy objects
  without that metadata do not assume 95%. Table/show methods honor their
  documented invisible return values without changing S4 slots.
* Factor scoring rejects unsupported mixture/multiple-group parent models clearly
  and preserves matrix dimensions for a single person or latent factor.
* New opt-in `getDiagnostics()` inspects stored covariance admissibility, Hessian
  and standard-error information, and mixture component sizes. It does not refit,
  impose constraints, or establish global optimality or repeated-start stability.
* The guides identify historical fitted objects and the inadmissible stored TVC
  example. LRT documentation states nesting, same-data and regularity requirements
  and removes an unqualified random-knot significance claim.
* Example-data reconstruction now has verified input hashes and named-variable
  mappings. R-devel and R-oldrel jobs are prepared for hosted checking. This remains
  a development candidate; scientific qualification and V1 publication are pending.

# nlpsem 1.0.0.9003 (development correctness update)

* Corrected mediation time prefixes to follow outcome, mediator, and predictor
  order, including distinct observation schedules in all model containers.
* Estimate and standard-error tables, including derived parameters and numeric
  confidence bounds, now retain full precision. Display p-values retain their
  existing format. Older saved estimate tables remain unchanged and may contain
  rounded values.
* Decomposed TVC models with growth TICs now determine free covariance parameters
  from model structure, independently of zero or nonzero starting values.
* Non-intrinsic EXP multivariate models now select covariance starts for every
  outcome. QUAD/JB multivariate LCSMs use outcome-specific constant matrices,
  allowing different numbers of waves.
* Multivariate residual covariance paths match shared wave IDs, including
  reordered selections. Outcome pairs without shared waves receive a clear
  error before optimization.
* Kappa intervals now include empirical marginal uncertainty using the
  multinomial delta variance for paired observed labels.
* Corrected class/group summary parsing for numbers above nine, and one-class
  mixture parameter output when there are no free mixture weights.
* Corrected the LRT example's full/reduced fixture order and added executable
  assertions for its degrees of freedom and finite p-value.
* Raised the minimum testthat version to 3.2.0, added targeted regression tests,
  and enabled the full fitting workflow on pull requests. Preserved the Phase C
  structural baseline with provenance and an explicit comparison tool.

# nlpsem 1.0.0.9002 (development version, Phase C)

* New fits from all seven fitting functions carry versioned metadata containing
  their structural specification, evaluated fitting settings, resolved starting
  values, random-number state, software versions, and optimizer diagnostics.
  `getSpecification()` exposes the metadata without duplicating input data or
  retaining the calling environment. Manifest-group metadata preserve original
  group labels and retained/dropped row indices.
* Added standard `print()`, `summary()`, `coef()`, `logLik()`, and `plot()` methods.
  Coefficients retain full precision; likelihood attributes support standard
  `AIC()` and `BIC()` calculations. Automatic plotting uses saved specifications
  and permits display changes without silently changing the fitted structure.
* Preserved the existing S4 slots and compatibility with older saved fits.
  Legacy fits return `NULL` from `getSpecification()` and require explicit
  structural plotting arguments. The concise S4 display uses the new print method.
* Added regression coverage for every fitting wrapper, serialization, original
  group labels, omitted rows, random-number-state preservation, standard methods,
  and legacy saved objects. Existing fitted examples and golden references are
  unchanged.
* Added a workflow vignette, a pkgdown website configuration, and GitHub Actions
  workflows for three-platform package checks, manual full fitting tests, and
  documentation artifacts. No website deployment or release publication is enabled.
* Added the full GPL version 3 license text for the source repository. The package
  license remains GPL version 3 or later. All 1.0.0.9001 correctness fixes are retained.

# nlpsem 1.0.0.9001 (development correctness update)

* Corrected multivariate outcome-pair indexing across builders and output helpers.
  Four or more outcomes no longer overwrite cross-process parameters or fail with
  conflicting starting values. Residual-correlation starts align shared wave IDs.
* Fixed class-specific nonparametric multivariate change-score output extraction
  and aligned automatic output labels with the numeric estimate/SE order.
* Plot confidence intervals now propagate the joint covariance of growth means,
  rates, knots and relevant derived parameters using the delta method. Intervals
  are pointwise 95% intervals for model-implied means, are not prediction intervals,
  and are no longer numerically clipped. At an estimated bilinear knot where the
  derivative is undefined, the pointwise interval is reported as missing.
* Longer outcome/time prefixes retain their observed-data plotting layers.
* Bootstrap confidence intervals exclude failed or non-finite replications and
  expose replication counts in the bootstrap table's `bootstrap_diagnostics`
  attribute. Fewer than two usable replications fail informatively. Confidence
  levels and integer replication counts are validated before estimation.
* `getLatentKappa()` matches categories by label values instead of factor positions.
  Different observed category sets use their union; class-label alignment remains
  an explicit user responsibility. Invalid confidence levels are rejected.
* Manifest-group derived changes use group-specific observation-time summaries;
  latent mixture classes retain a pooled reference grid. Selected nonconsecutive
  waves use the preceding selected occasion and positional loading indices.
* Added regression tests for the corrected calculations, wider input configurations,
  real four-outcome and class-specific fits, and independent uncertainty calculations.
  Five nonconverged Phase B golden references are now explicit exclusions, and a
  coverage check prevents the exclusion set from silently growing. Existing fitted
  examples and numerical reference files are unchanged.
* Preserved the 688-specification structural snapshot harness in `data-raw/`.
* Corrected the full/reduced-model order in the `getLRT()` documentation example.

# nlpsem 1.0.0 (development version, Phase B)

## Structure

* One builder per sub-model. `.build_LGCM()`, `.build_LCSM()`, `.build_TVC()`, `.build_MGM()` and `.build_MED()`
  build one `mxModel` per prefix: the single-group model with an empty prefix, the class-specific models of a
  multiple-group or mixture model with `c1`, `c2`, .... `get<Model>.mxModel()` and `getsub.<Model>(container =)`
  are thin wrappers; the container contributes only its data rows (multiple-group) or the row-wise fit function
  (mixture). The loading, derived-parameter and TVC-component helpers are parameterised the same way, and the
  multiple-group and mixture output helpers share `.class_output()`. The ten `MGroup.SUBMODEL.*`/`MIX.SUBMODEL.*`
  builders, six `MIX_*` helpers and the second output helper are gone: `R/` shrinks from 27,525 to 14,692 lines.
  A structural snapshot of 688 model specifications (every sub-model, curve form, intrinsic option, TIC and
  decomposition setting, single-group, multiple-group and mixture with and without cluster TICs) against the
  0.4.1 builders shows identical matrices, free parameters and algebras apart from the items listed below.

## Bug fixes

* `getMGroup()`/`getMIX()` with `sub_Model = "TVC"`, `y_model = "LCSM"`, `curveFun = "EXP"`, `intrinsic = TRUE`,
  `decompose != 0` and no growth TICs could not be fitted ("Unknown reference 'c1Y_alpha0'"): the class-specific
  builder branch was a copy of the wrong case. The unified builder uses the single-group model for this case.
* The mixture weight rows of `getMIX()` tables (`p2`, ...) now carry the standard errors of the softmax weights;
  0.4 repeated the estimates in the SE column.
* `getLRT()` returns numeric p-values; 0.4 replaced values below 0.0001 by the string `"<0.0001"`, which turned the
  column into text.
* `getFigure()` for multiple-group Jenss-Bayley change-score models used the lower-bound acceleration terms in the
  upper bound of the growth-rate interval and vice versa.

* Data-driven residual-correlation starts for multivariate and decomposed-TVC models are now computed from the
  per-person curve residuals (`getUNI.GF()` attaches them as the `"residuals"` attribute) instead of
  person-mean-centred scores. The 0.4 heuristic left the growth trend in the deviations, so two processes that grew
  together started from a residual correlation above 0.9 (the fitted value on the example data is 0.24); with that
  start `getMIX()` could not recover the three-class bivariate example. Mediation models keep the builders' default
  start of 0.3 (an unused copy of the heuristic was removed from the mediation start helper).

* `getMIX()` with `sub_Model = "MGM"` now damps the class-specific fixed-knot bilinear-spline covariance starts as
  `getMGroup()` already did (each provisional class is a smaller sample whose data-driven transformed covariance is
  too large to start from). With both start changes the three-class bivariate example converges at the first
  optimisation attempt.

* `getMGroup()` and `getMIX()` with `sub_Model = "MED"` and a baseline predictor (`x_type = "baseline"`)
  could not be fitted unless the predictor column was literally named `X`: the class-specific mediation
  builders hard-coded `"X"` in the paths for the predictor's mean, variance and effects. They now use
  `x_var`, as the single-group builder always did. Found by the new fitting-layer tests.


## Changes in labels and names

* Single-group bilinear-spline LGCMs label the knot `Y_knot` (was `knot`) and the reparameterised TIC coefficients
  `beta01`, ... (was `beta01s`, ...), matching the multiple-group and mixture models. This affects
  `omxGetParameters()` names and `mxCI()` labels only; estimate tables are unchanged.
* Mediation models name the latent growth factors `eta1Y`, `etarY`, `eta2Y` (the single-group convention) in all
  containers; the class-specific models used `etaYr`, `etaY2`.
* The TVC absolute-rate algebras of class-specific decomposed-TVC models are named `c<k>X_abs_rate<j>` (was
  `c<k>abs_rate<j>`), as in the single-group models.
* Single-group TVC model names state the decomposition ("LGCM with decomposed TVC, ...") whenever `decompose != 0`,
  and the two multivariate negative-exponential change-score models state the ratio type.

## Reproducibility and testing

* The fitted example objects in `inst/extdata` were regenerated with the current code (R 4.5.3, OpenMx 2.22.11) by
  the scripts now kept in `data-raw/`, which follow the vignette calls. Twenty of the twenty-three fits reproduce the
  previous log-likelihoods exactly. `Math_TVCslp_BLS_f` now reaches a lower -2LL (42269.2 instead of 42457.3), while
  `RM_PLGCM.f` (62609.2 instead of 62589.5) and `RM_BLS_PLGCM3` (61903.6 instead of 61865.3, and now labelled with the
  process names R/M instead of the pre-0.4 Y/Z) converge to nearby local optima of these multimodal likelihoods from the
  new start values; `inst/extdata/REGENERATED.txt` records the comparison. As with any random-knot or mixture model,
  compare several start sets (`prop_starts`, `res_cor`, `starts`) before interpreting a single fit.

* The scripts that regenerate the fitted example objects in `inst/extdata`, the derivation of
  `RMS_dat` from the public-use ECLS-K:2011 files, and the generator of the fitting-layer reference
  values now live in `data-raw/` (excluded from the built package).
* Added fitting-layer regression tests (`tests/testthat/test-fitting_golden.R`): every sub-model and
  functional form is fitted on a subsample and compared with stored reference values. They are
  skipped on CRAN and take several minutes locally.
* Examples run in a few minutes instead of 37: post-processing examples load the shipped fits, the
  fitting examples keep one fast model each, and the long intrinsically nonlinear and mixture fits are
  shown but not run.

---

# nlpsem 0.4.1

## Bug fixes

* `getEstimateStats()`: An unrecognized `CI_type` now stops with an informative error
  listing the supported options (`"Wald"`, `"likelihood"`, `"bootstrap"`, `"all"`)
  instead of failing with `object 'statOut' not found`.
* `getEstimateStats()`: Requesting bootstrap CIs (`CI_type = "bootstrap"` or `"all"`)
  without supplying a positive integer `rep` now stops with an informative error
  instead of passing `NA` to `mxBootstrap()`.
* `getEstimateStats()`: `conf.level` is now honoured for likelihood-based and bootstrap
  confidence intervals. Both were fixed at 95% regardless of the argument (`mxCI()` was
  called without `interval`, and the bootstrap quantiles were hard-coded at 2.5% and 97.5%).
* `getFigure()` now works for longitudinal mediation models (`sub_Model = "MED"`), with
  both baseline and longitudinal predictors, in single-group, multiple-group, and mixture
  fits. The plotting helper looked up growth-factor means under the user's variable prefix
  with a `_mean0` suffix, whereas the mediation builders name them by role (`Y_mean`,
  `M_mean`, `X_mean`) and store bilinear spline growth factors as (slope 1, value at knot,
  slope 2) with the knot in `<role>_mug`. Every release up to 0.4 therefore failed with
  `object '<prefix>_mean0' not found`; the 0.4 release note on `getFigure()` for mediation
  models described only the trajectory-loop fix.
* `getIndFS()`: `FS_type = "ML"` and `FS_type = "WeightedML"` now work. `mxFactorScores()`
  failed whenever the data stored in the fitted model contained columns the model does not
  use (covariates, identifiers), and the descriptive model names used by this package are
  not valid identifiers for the scoring model. The data are now reduced to the model's
  manifest and definition variables and the model is renamed before scoring; `"Regression"`
  scores are unchanged.
* `getMIX()` with `paramOut = TRUE`, `cluster_TIC`, and three or more latent classes: the
  class-membership coefficients in the estimates table were mislabelled. Labels were listed
  class by class (`beta20`, `beta21`, ..., `beta30`, ...) while the values were read from the
  coefficient matrix column by column, so for example the row named `beta21` held the value of
  `beta30`. Introduced in 0.4 with the generalised cluster-TIC labelling; two-class models were
  unaffected. Values are now read row by row to match the labels.
* Multiple-group and mixture models with a time-varying covariate and a random-knot
  bilinear spline outcome (`getMGroup()` / `getMIX()` with `sub_Model = "TVC"`,
  `curveFun = "BLS"`, `intrinsic = TRUE`): the class-specific builders set the knot entry of
  the delta-method Jacobian to 0 instead of `mueta2s - mueta1s`, so the outcome's growth-factor
  variances and covariances in the original parameter space (and their standard errors and the
  plotted confidence bands) omitted the knot-variance contribution. The single-group TVC builder
  and all LGCM and MGM builders already had the term.
* `getMGroup()` accepts any coding of the grouping variable. Downstream helpers index the groups
  `1..K`, so groups coded `0/1`, as characters, or as factors failed with `subscript out of
  bounds`. The distinct values are now sorted and recoded to `1..K` (a message reports the
  mapping) and rows with a missing group value are dropped with a warning.
* Initial values tolerate persons with too few observations. The per-person linear and
  quadratic fits in `getUNI.GF()` crashed (`0 (non-NA) cases`) for a person with fewer than two
  (three) complete occasions; multivariate and mediation starts produced `NA` means and
  covariances from such persons; and the k-means step that forms provisional classes for
  `getMIX()` could not run with incomplete growth factors. Such persons now contribute `NA`
  growth factors that are skipped in the means and covariances, and they are assigned to the
  nearest provisional class.
* `getSummary(..., HetModels = TRUE)` returns a plain `data.frame`, as documented and as
  `HetModels = FALSE` already did, instead of a tibble.
* Standard-error lookups through `mxSE()` are now silent. Previously every `getFigure()` call
  and every fit with `paramOut = TRUE` printed OpenMx's "Treating first argument ..." message
  once per looked-up quantity (dozens of lines for multiple-group and mixture models).

## Package

* Added `inst/CITATION` for Liu (2025), *Behavior Research Methods*, 57(3), 87,
  <doi:10.3758/s13428-025-02596-4>.
* The fitted example objects in `inst/extdata` are stored with xz compression
  (3.7 MB to 1.1 MB), which brings the installed package size below 5 MB.

## Documentation

* Corrected the 0.4 release note describing the improved BLS initial values: the
  shipped implementation uses an OLS-based grid search over candidate knot locations;
  the NLS refinement step described previously was replaced by closed-form OLS fits
  before release.
* Removed stale `nls`/`nls.control` imports left over from the pre-grid-search BLS
  starts implementation.

## Testing

* Added regression tests for post-processing functions (`getEstimateStats()`,
  `getSummary()`, `getPosterior()`, `getLatentKappa()`) using synthetic inputs and
  the fitted example objects shipped in `inst/extdata`.
* Added regression tests driven by the shipped example fits for `getFigure()` across all
  seven sub-models (including both mediation designs), the output helpers and automatic
  parameter names for every shipped configuration, `getIndFS()` for all three score types,
  `getLRT()`, `getSummary()`, `getPosterior()` with cluster-forming covariates, and the
  `conf.level` handling in `getEstimateStats()`.
* Added tests that build (without fitting) class-specific TVC models and check the
  delta-method Jacobian, tests for the group recoding, and tests for starts computed from data
  with persons who have too few observations.

---

# nlpsem 0.4

This release note reflects a direct comparison between the current package and
`package_v0.3/nlpsem`.

## Code quality and documentation improvements

### Bug fixes

* Fixed division-by-zero in `getPosterior()` entropy calculation when `nClass = 1`.
* Fixed typo "Accuracys" → "Accuracy" in posterior output display.
* Fixed missing closing parenthesis in `getTVCmodel()` documentation for Jenss-Bayley curve function.
* Fixed typo `memebership` → `membership` in `getSummary()`.
* Added missing NULL check for `class` column in `getPosterior()` when `label = TRUE`.
* Added missing `model` NULL validation in `getEstimateStats()` for likelihood, bootstrap, and combined CIs.
* Fixed nested `any()` logic error in `validate_res_scale()`.
* Fixed `[[]]` list indexing on numeric vector in `getMULTI.initial()`.
* Fixed typo "dose not" → "does not" in `TVC.initial_helper.R`.
* Fixed `getFigure()`: Mediation sub-model loop iterated over `y_var` (always length 1)
  instead of `traj_var`, and `outcome` was not indexed — only the first trajectory was plotted.
* Fixed `TVC.initial_helper.R`: Hardcoded `t_var = "T"` replaced with the actual `t_var`
  function parameter when computing TVC initial values.
* Fixed `MED.output_helper.R`: Swapped y_var/m_var standard error order in 3 output blocks
  (BLS baseline, LIN longitudinal, BLS longitudinal) so SEs match their corresponding estimates.
* Fixed `MED.mxModel_helper.R`: Removed quoted strings from `mxAlgebra()` calls in 4 places.
  `mxAlgebra()` uses non-standard evaluation (bare names); quoted strings are for
  `mxAlgebraFromString()`.
* Fixed array indexing in `TVC.output_helper.R`, `MIX.output_helper.R` (2 places), and
  `MGroup.output_helper.R` (2 places): Loop assumed consecutive records starting at 1;
  replaced with `seq_along()` for correct indexing with arbitrary record indices.
* Fixed `getPosterior()`: Protected entropy calculation from `log(0)` using
  `pmax(prob, .Machine$double.xmin)`.
* Fixed `getLatentKappa()`: Added defensive guard for kappa edge cases (kappa >= 1 or <= 0).
* Fixed `getSummary()`: Empty latent classes in mixture models now handled correctly using
  `factor(membership, levels = 1:nClass)` to ensure proportions vector matches `nClass`.
* Fixed `validate_res_scale()`: List branch now uses `<= 0` / `>= 1` to match non-list branch
  (both enforce exclusive bounds as documented).
* Fixed `MIX.initial_helper.R`: TVC residual (`TVC_starts$residuals`) now also inflated by 2×
  when `decompose != 0`, matching the Y residual inflation for mixture model convergence.
* Fixed `MULTI.initial_helper.R`: Residual covariance matrix off-diagonal indexing used
  `res_cor[traj_i + traj_j - 1]` which produced wrong indices for 4+ longitudinal outcomes.
  Replaced with sequential counter to correctly map pair (i,j) to `res_cor` vector position.
* Fixed `TVC.initial_helper.R`: Added `res_cor = NULL` default to `getTVC.initial()` signature.
  Callers in `MIX.initial_helper.R` and `MGroup.initial_helper.R` omit `res_cor` when
  `decompose = 0` (where it is unused), but the missing default could cause errors if the
  function is refactored.
* Fixed `getMGM()` / `getMULTI.initial()`: Multivariate models now preserve the full `y_var`
  prefixes provided by the user instead of collapsing them to single-letter uppercase prefixes.
* Fixed `MGM.output_helper.R`, `MGroup.output_helper.R`, and `MIX.output_helper.R`:
  between-process covariance output for 3+ outcomes/classes no longer appends cumulative
  estimate blocks or corrupts `paramOut` tables.
* Fixed `MGM.mxModel_helper.R`, `MGroup.SUBMODEL.MGM_helper.R`, and
  `MIX.SUBMODEL.MGM_helper.R`: Fixed-parameter multivariate covariance extraction now works
  for 3+ outcomes instead of assuming exactly 2 longitudinal processes.
* Fixed `MED.Initial_helper.R`: In longitudinal-predictor mediation models, the initial
  `beta_x` and `beta_m` blocks are no longer swapped in the outcome regression starts.
* Fixed `TVC.initial_helper.R`: `decompose = 3` now computes `kappa` starts from the
  change-from-baseline temporal state implied by the latent growth-factor starts, and the
  nonparametric path now returns finite starts instead of falling through to `NA`.
* Fixed intrinsic bilinear spline loading expressions across LGCM, TVC, MGM, mixture, and
  multiple-group models by replacing the nondifferentiable knot expression with a smooth
  approximation, improving optimizer behavior near the knot.
* Fixed `MIX.output_helper.R`: Cluster-TIC coefficient names now scale to any number of
  predictors instead of hard-coding a 3-column naming pattern.

### Code improvements

* Removed unnecessary `I()` wrappers from logical conditions across 11 files.
* Extracted duplicated validation logic into shared internal helpers (`validation_helpers.R`).
* Vectorized row-by-row probability matrix computation in `getPosterior()` using `sweep()`/`rowSums()`.
* Refactored positional list indexing (`starts[[1]][[2]]`) to named access (`starts$Y_starts$covMatrix`) across all helper files for maintainability.
* Added `as.data.frame(dat)` guard in all main functions for tibble compatibility.
* Added `validate_res_cor()` helper and validation calls in `getTVCmodel()`, `getMGM()`, `getMediation()`, `getMIX()`, and `getMGroup()`.
* Removed duplicate `validate_curveFun_mediation()` function (logic already in `validate_sub_Model()`).
* Eliminated code duplication in `getFigure()` by removing redundant `is.null(nClass)` / `!is.null(nClass)` branching.
* Simplified redundant `else if` conditions in `StandardMethods.R`.
* **Data-driven residual initial values**: `getUNI.GF()` now computes per-person mean squared
  residuals from individual model fits and returns the median as an attribute. All initial value
  helpers (`UNI.initial_helper.R`, `MULTI.initial_helper.R`, `TVC.initial_helper.R`,
  `MED.Initial_helper.R`) now use this data-driven estimate instead of the `var() * res_scale`
  heuristic, with the heuristic retained as a fallback. This provides more accurate starting
  values, especially for nonlinear functional forms (EXP, JB, BLS) where the heuristic could
  be orders of magnitude off.
* Refactored duplicated Wald/likelihood/bootstrap CI code in `getEstimateStats()` into local
  helpers (`.add_wald_CI()`, `.get_likelihood_CI()`, `.get_bootstrap_CI()`).
* Changed scalar logical operators from `&`/`|` to `&&`/`||` in `validation_helpers.R` (7 places)
  and `getFigure.R` (3 places) for proper short-circuit evaluation.
* Removed dead code in `MIX.initial_helper.R` (duplicate `order_indices` assignment).
* Added 13 NSE global variable declarations to `global_vars.R` for bare names used in
  `mxAlgebra()` calls.
* **`res_scale` is now optional**: With data-driven residual estimation, users no longer need
  to specify `res_scale`. When omitted (default `NULL`), the data-driven estimate is used;
  if that fails, a fallback heuristic of `var() * 0.1` is applied. Users can still override
  with an explicit `res_scale` value if desired.
* **`res_cor` is now optional**: Residual correlation between longitudinal outcomes is now
  estimated from data when `res_cor` is not specified. The estimation uses within-person
  demeaned residuals correlated across trajectory pairs. If estimation fails, a heuristic of
  `0.3` is applied as fallback. This applies to `getMGM()`, `getMediation()`,
  `getTVCmodel()` (when `decompose != 0`), `getMIX()`, and `getMGroup()`. Users can still
  override with explicit `res_cor` values if desired.
* **Automatic parameter naming**: When `paramOut = TRUE`, the `names` argument is now optional.
  If omitted (default `NULL`), meaningful parameter names are automatically generated based on
  the model configuration (functional form, intrinsic flag, covariates, decompose level, etc.).
  This eliminates the need for users to manually construct long parameter name vectors for every
  model. Users can still supply custom names to override the auto-generated labels. Auto-naming
  is supported across all 7 model functions: `getLGCM()`, `getLCSM()`, `getTVCmodel()`,
  `getMediation()`, `getMGM()`, `getMIX()`, and `getMGroup()`. A safety net in all output
  helpers auto-pads with generic labels if the name count doesn't match the estimate count.
* **Improved BLS initial values**: Bilinear spline (BLS) individual growth factor estimation
  now uses a grid search over candidate knot locations (interval midpoints, thinned to at
  most 15 quantiles) with closed-form OLS piecewise-linear fits, anchored by a pooled-data
  fit, instead of a single midpoint guess refined by NLS. Replacing NLS with OLS grid search
  eliminates NLS convergence failures entirely and yields more accurate starting values for
  BLS models. (Corrected in 0.4.1: this note previously described an intermediate
  implementation that retained an NLS refinement step.)
* Added pooled-anchor nonlinear starts for negative exponential, Jenss-Bayley, and bilinear
  spline growth functions. These starts preserve person-level pseudo growth factors, so latent
  mean and covariance starts are still estimated from subject-specific values rather than from
  a single pooled trajectory.
* Added intrinsic-specific stabilization for nonlinear starts in univariate, TVC, multivariate,
  multiple-group, and mixture wrappers, including more conservative starts for random
  knot/ratio/acceleration components.
* Added robust cluster-TIC softmax starts for mixture models: predictors are standardized
  before the multinomial fit, coefficients are transformed back to the original scale, and the
  result is shrunk toward intercept-only logits with a fallback for unstable fits or constant
  predictors.
* Centralized repeated positive-definite repair logic in `initial_helpers_common.R`, added
  reusable validation helpers in `validation_helpers.R`, and introduced shared helpers for
  output naming and multivariate covariance reduction in `output_helpers_common.R` and
  `multi_cov_helpers.R`.

### Documentation

* Added `@seealso` cross-references linking related functions across all 14 exported functions.
* Added `@noRd` to all 44 internal helper files to suppress unnecessary `.Rd` generation.
* Added introductory vignette (`nlpsem_overview`) with package overview and quick-start guide.
* Expanded `getMGroup_examples` vignette with a second example (quadratic LGCM with TICs).
* Added explanatory comment to `global_vars.R`.
* Added `NEWS.md` changelog.
* Added `\examples{}` sections to `ModelSummary` and `printTable` S4 generics.
* Fixed LaTeX math mode (`$95\%$`) to plain text in `KappaOutput` class documentation.
* Shortened DESCRIPTION title to comply with CRAN < 65 character guideline.
* Fixed typo "loadigs" → "loadings" in comments across 10 helper files.
* Fixed typo "his function" → "This function" in `LGCM.mxModel_helper.R`.
* Fixed typo "model.It" → "model. It" in `getEstimateStats.R` roxygen.
* Fixed typo "Calcualte" → "Calculate" in `getLatentKappa.R` example.
* Fixed grammar "if want to" → "if you want to" in `validation_helpers.R`.
* Fixed misleading comments in `getMIX.R` and `getMGroup.R` (said "latent change score model"
  instead of "mixture model" / "multiple group model").
* Synchronized the vignette code with the regenerated examples in `inst/extdata`.
* Updated package- and function-level references from preprint / "advance online publication"
  wording to the final journal citations where available.
* Added vignette build metadata to `DESCRIPTION` (`knitr`, `rmarkdown`, and `VignetteBuilder`)
  and added the introductory `nlpsem_overview` vignette.

### Testing and infrastructure

* Added a `tests/testthat` suite covering validation helpers, automatic parameter naming,
  nonlinear growth-factor starts, mediation starts, and multivariate regression cases.
* Refreshed the saved example objects in `inst/extdata/*.RData` and aligned the vignettes with
  those regenerated examples.

---

# nlpsem 0.3

## Initial CRAN release

### Model types

* Univariate latent growth curve models (`getLGCM()`) and latent change score
  models (`getLCSM()`), with optional time-invariant covariates (TICs).
* Models with time-varying covariates (`getTVCmodel()`).
* Multivariate longitudinal models (`getMGM()`) and longitudinal mediation
  models (`getMediation()`).
* Multiple-group models (`getMGroup()`) for comparing trajectories across
  manifested groups.
* Longitudinal mixture models (`getMIX()`) for identifying latent trajectory
  classes.

### Supported growth curve functions

* Linear, quadratic, negative exponential, Jenss-Bayley, and bilinear spline
  functional forms.
* Intrinsically nonlinear specification for negative exponential, Jenss-Bayley,
  and bilinear spline models.

### Post-processing and inference

* `getSummary()`: Model fit statistics (number of parameters, -2LL, AIC, BIC).
* `getEstimateStats()`: P-values and confidence intervals (Wald,
  likelihood-based, and bootstrap).
* `getFigure()`: Visualization of estimated growth trajectories with confidence
  intervals.
* `getLRT()`: Likelihood ratio test for nested model comparison.
* `getIndFS()`: Individual factor scores.
* `getPosterior()`: Posterior class probabilities for mixture models.
* `getLatentKappa()`: Cohen's kappa for latent class agreement.

### Data

* Included `RMS_dat`, a processed sample from the ECLS-K:2011 study for
  demonstration purposes.
