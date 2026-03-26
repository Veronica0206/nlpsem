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
  now uses a grid search over 10 candidate knot locations (with OLS piecewise-linear fits)
  instead of a single midpoint guess. The top 3 candidates are tried with NLS, using both
  default Gauss-Newton and port algorithms with bounded knot location. This reduces NLS
  failure rates from ~40-46% to ~25% and lowers median residual variance by ~40-60%,
  yielding more accurate starting values for BLS models.
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
