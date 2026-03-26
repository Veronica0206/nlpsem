# nlpsem

`nlpsem` provides tools for fitting linear and nonlinear longitudinal models
within the structural equation modeling (SEM) framework. The package is built
on top of `OpenMx` and is designed for settings where the trajectory shape is
nonlinear and may be intrinsically nonlinear rather than a simple
reparameterization of a linear growth model.

The methods implemented in this package are described in:

Liu, J. (2025). *Examination of nonlinear longitudinal processes in the
framework of structural equation modeling: The R package nlpsem*. Behavior
Research Methods. <https://doi.org/10.3758/s13428-025-02596-4>

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

If the package is available on CRAN:

```r
install.packages("nlpsem")
```

To install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("Veronica0206/nlpsem")
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

## What is new in version 0.4

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

See [NEWS.md](NEWS.md) for a fuller summary of changes from version 0.3.

## Documentation

The package includes an overview vignette and scenario-specific example
vignettes:

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
