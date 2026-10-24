# data-raw

Scripts that produce the data and fitted example objects shipped with the package. This folder is
excluded from the built package (`.Rbuildignore`).

| Script | Produces |
|---|---|
| `RMS_dat.R` | `data/RMS_dat.RData` from the public-use ECLS-K:2011 files (not redistributed). |
| `regenerate_extdata.R` | Driver: loads the package and runs the seven `regen_*.R` scripts. |
| `regen_getLGCM.R` ... `regen_getMIX.R` | One `inst/extdata/*_examples.RData` each; used by the vignettes and the tests. |
| `make_golden_fits.R` | `tests/testthat/golden/fit_golden.rds`, the reference values for the fitting-layer tests. |

Regenerate the example fits only together with a release, and compare the new objects with the old ones
(log-likelihoods and parameter estimates) before replacing them: the vignette text quotes numbers from
them. `inst/extdata/REGENERATED.txt` records the R and OpenMx versions used.

The bivariate random-knot model (`RM_PLGCM.f`) and the three-class bivariate mixture (`RM_BLS_PLGCM3`) have
multimodal likelihoods: in October 2026 a sweep over residual-correlation starts (`res_cor` 0.1 to 0.95) for
`RM_PLGCM.f` ended at -2LL values between 62505 and 62609 (all status 0 or 6), and `RM_BLS_PLGCM3` reached the same
61903.6 solution from the 1.0 default starts and from the 0.3-era calls. The regenerated objects are the fits the
documented calls produce; `REGENERATED.txt` lists how they compare with the previous objects.
