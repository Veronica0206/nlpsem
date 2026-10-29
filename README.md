# nlpsem

[![CRAN](https://www.r-pkg.org/badges/version/nlpsem)](https://cran.r-project.org/package=nlpsem)
[![GitHub repo](https://img.shields.io/badge/GitHub-nlpsem-black)](https://github.com/Veronica0206/nlpsem)
[![R](https://img.shields.io/badge/R-OpenMx%20package-blue)](https://www.r-project.org/)
[![DOI](https://img.shields.io/badge/DOI-10.3758%2Fs13428--025--02596--4-blue)](https://doi.org/10.3758/s13428-025-02596-4)
[![License: GPL >= 3](https://img.shields.io/badge/License-GPL--3%2B-green.svg)](https://github.com/Veronica0206/nlpsem/blob/main/DESCRIPTION)

`nlpsem` provides tools for fitting linear and nonlinear longitudinal models
within the structural equation modeling (SEM) framework. The package is built
on top of `OpenMx` and is designed for settings where the trajectory shape is
nonlinear and may be intrinsically nonlinear rather than a simple
reparameterization of a linear growth model.

## Development 1.0.0.9006

Quadratic, negative exponential and Jenss-Bayley models are fitted to standardized
data and converted back exactly: the estimates, standard errors and likelihood are
reported in the original units, and `getSpecification()` records the scaling.
`options(nlpsem.standardize = FALSE)` fits on the raw scale. Other functional forms
and mediation models are fitted as before.

The 1.0.0.9005 update made EXP, Jenss-Bayley and bilinear spline starting values use per-person
`stats::nls()` with the 0.3 curve formulas and default Gauss-Newton algorithm.
The shared initializer supplies this method to all model families. Failed NLS
fits use recorded heuristic fallbacks; they are not counted as NLS convergence.
The final SEM optimizer is unchanged. The matched comparison protocol explicitly
uses CSOLNP and identical saved unfitted models on each platform.

The preceding development update fixes mediation wave/covariance starts, one-class mixture
post-processing, posterior normalization and class-count validation, multivariate
summary columns, parameter names and confidence-interval labels. It adds opt-in
`getDiagnostics()` checks for stored covariance admissibility, Hessian and standard
error information, and mixture component sizes. It does not refit or constrain models.

```r
print(fit)
coef(fit)
getDiagnostics(fit)
getSpecification(fit)$specification
plot(fit, draw = FALSE)
```

The [NLS validation record](.github/VALIDATION-1.0.0.9005.md) separates completed
checks from outstanding numerical qualification. The fitting-layer reference values
were regenerated with 1.0.0.9006; each fit is compared with its reference only after
it passes the `getDiagnostics()` checks (see `NEWS.md`). The saved examples were
regenerated with 1.0.0.9006 (`inst/extdata/REGENERATED.txt`); `getDiagnostics()` flags four
of them, which the guides note. This development version is not a V1.0.0 release.

New estimate tables retain full precision. Older serialized tables retain their
historical rounding and need recomputation or refitting to recover precision. S4
slots are preserved. Invalid custom names now fail explicitly instead of recycling.
Factor scoring supports single RAM/LISREL models with their own raw data; mixture
and multiple-group parent objects are rejected with a clear explanation.

## Portfolio positioning

`nlpsem` is a statistical software and measurement-methods asset. It supports nonlinear longitudinal modeling, latent growth and latent change models, mixture modeling, multivariate longitudinal processes, multiple-group comparisons, and psychometric/measurement-science workflows in R.

For AI evaluation roles, the package is relevant because repeated AI behavior across prompts, models, seeds, raters, and deployment conditions can be studied with the same measurement mindset used for longitudinal human-development and latent-process modeling: define the measurement occasions, model systematic variation, quantify reliability, and separate signal from instability.

## For hiring reviewers

This repository demonstrates:

- R package development with documented functions, vignettes, tests, and CRAN-facing metadata
- nonlinear longitudinal modeling in the SEM framework
- latent growth, latent change, multivariate, mediation, multiple-group, and mixture modeling workflows
- psychometric and measurement-science foundations relevant to evaluator reliability and repeated-behavior analysis
- reproducible statistical examples backed by `OpenMx`

The methods implemented in this package are described in:

Liu, J. (2025). *Examination of nonlinear longitudinal processes with latent variables, latent processes, latent changes, and latent classes in the structural equation modeling framework: The R package nlpsem*. Behavior
Research Methods. <https://doi.org/10.3758/s13428-025-02596-4>

## Connection to AI evaluation

The same measurement-science principles used for longitudinal human-development modeling also inform evaluation of AI system behavior over prompts, models, seeds, raters, and deployment conditions. In both settings, the central questions are how behavior changes across repeated conditions, which sources of variation are meaningful, and how much confidence a reviewer should place in an observed pattern.

## Main modeling scenarios

The package supports four broad longitudinal modeling settings:

1. Univariate longitudinal models with growth factors, with or without
   time-invariant covariates (TICs) and time-varying covariates (TVCs)
2. Multivariate longitudinal models for assessing correlated or potentially
   directional development across multiple processes
3. Multiple-group longitudinal models for comparing manifested groups
4. Longitudinal mixture models for identifying latent trajectory classes

## Core functions

### Model-fitting functions

- `getLGCM()`: latent growth curve models
- `getLCSM()`: latent change score models
- `getTVCmodel()`: models with time-varying covariates
- `getMGM()`: multivariate growth and change models
- `getMediation()`: longitudinal mediation models
- `getMGroup()`: multiple-group models
- `getMIX()`: longitudinal mixture models

### Post-processing functions

- `getSummary()`: model fit summary across one or more fitted models
- `getEstimateStats()`: p-values and confidence intervals
- `getFigure()`: trajectory plots
- `getLRT()`: likelihood ratio tests for nested models
- `getIndFS()`: individual factor scores
- `getPosterior()`: posterior class probabilities for mixture models
- `getLatentKappa()`: agreement between latent class assignments

## Supported functional forms

For latent growth curve models, `nlpsem` supports:

- linear (`"LIN"`)
- quadratic (`"QUAD"`)
- negative exponential (`"EXP"`)
- Jenss-Bayley (`"JB"`)
- bilinear spline (`"BLS"`)

For latent change score models, `nlpsem` supports:

- nonparametric (`"NonP"`)
- quadratic (`"QUAD"`)
- negative exponential (`"EXP"`)
- Jenss-Bayley (`"JB"`)

Negative exponential, Jenss-Bayley, and bilinear spline models can be fit in
their intrinsically nonlinear forms by setting `intrinsic = TRUE`.

## Installation

The published version is available from CRAN:

```r
install.packages("nlpsem")
```

To install the development source from a local checkout:

```r
# install.packages("remotes")
remotes::install_local(".")
```

## Quick example

```r
library(nlpsem)

mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)

data("RMS_dat")
RMS_dat0 <- RMS_dat

# Re-baseline the measurement occasions
baseT <- RMS_dat0$T1
for (i in 1:9) {
  RMS_dat0[[paste0("T", i)]] <- RMS_dat0[[paste0("T", i)]] - baseT
}

# Fit a bilinear spline LGCM with a fixed knot
Math_BLS <- getLGCM(
  dat = RMS_dat0,
  t_var = "T",
  y_var = "M",
  curveFun = "BLS",
  intrinsic = FALSE,
  records = 1:9,
  paramOut = TRUE
)

# Inspect estimates
printTable(Math_BLS)
```

## Earlier changes in version 0.4

Version 0.4 keeps the modeling scope of version 0.3 while substantially
improving robustness, usability, and documentation.

Highlights include:

- stronger initial-value calculations for nonlinear forms, especially `EXP`,
  `JB`, and `BLS`
- smoothed intrinsic bilinear spline loadings for better optimizer behavior
  near the knot
- data-driven residual variance and residual correlation starts
- automatic parameter naming when `paramOut = TRUE`
- expanded validation and cleaner handling of multivariate, mediation, TVC,
  mixture, and multiple-group edge cases
- synchronized vignettes and saved example objects
- a new `testthat` test suite and cleaner CRAN-facing documentation

See `NEWS.md` for a more comprehensive summary of changes from version 0.3.

## Documentation

The package includes an overview vignette and scenario-specific example
vignettes:

- `vignette("phase_c_workflow")`
- `vignette("nlpsem_overview")`
- `vignette("getLGCM_examples")`
- `vignette("getLCSM_examples")`
- `vignette("getTVCmodel_examples")`
- `vignette("getMGM_examples")`
- `vignette("getMediation_examples")`
- `vignette("getMGroup_examples")`
- `vignette("getMIX_examples")`

## License

GPL (>= 3)
