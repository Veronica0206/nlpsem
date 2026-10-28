# data-raw

Scripts that produce the data and fitted example objects shipped with the package. This folder is
excluded from the built package (`.Rbuildignore`).

| Script | Produces |
|---|---|
| `RMS_dat.R` | `data/RMS_dat.RData` from the public-use ECLS-K:2011 files (not redistributed). |
| `regenerate_extdata.R` | Driver: loads the package and runs the seven `regen_*.R` scripts. |
| `regen_getLGCM.R` ... `regen_getMIX.R` | One `inst/extdata/*_examples.RData` each; used by the vignettes and the tests. |
| `make_golden_fits.R` | `tests/testthat/golden/fit_golden.rds`, the reference values for the fitting-layer tests. |
| `approximation-bias/run_approximation_bias.R` | The approximation sizes quoted in the help of `getSimData()`, `getLCSM()` and `getTVCmodel()`; results and design in `approximation-bias/`. |

Regenerate the example fits only together with a release, and compare the new objects with the old ones
(log-likelihoods and parameter estimates) before replacing them: the vignette text quotes numbers from
them. `inst/extdata/REGENERATED.txt` records the R and OpenMx versions used.

The bivariate random-knot model (`RM_PLGCM.f`) and the three-class bivariate mixture (`RM_BLS_PLGCM3`) have
multimodal likelihoods: in October 2026 a sweep over residual-correlation starts (`res_cor` 0.1 to 0.95) for
`RM_PLGCM.f` ended at -2LL values between 62505 and 62609 (all status 0 or 6), and `RM_BLS_PLGCM3` reached the same
61903.6 solution from the 1.0 default starts and from the 0.3-era calls. The regenerated objects are the fits the
documented calls produce; `REGENERATED.txt` lists how they compare with the previous objects.

## Verified data reconstruction

`reconstruct_RMS_dat.R` uses `RMS_dat-variable-map.csv` and verifies SHA256 hashes
from `RMS_dat-provenance.json` before reading the three local public-use SAS exports.
It selects named variables, translates the documented missing codes, joins by unique
CHILDID, sorts by ID, and samples 500 from 1,817 complete cases with seed 20191029.
It reproduced all 500 rows and 49 variables of the shipped data under R 4.5.3.
The script requires haven, digest and jsonlite and writes only to a new output folder.
The external NCES release identifier was not retained with the historical exports;
the provenance file records this limit instead of inferring a release from filenames.
No raw SAS export is added to the package.

```sh
Rscript data-raw/reconstruct_RMS_dat.R /path/to/source/exports /path/to/new/output
```

The source archive excludes data-raw. The dataset help explains the complete-case
teaching sample; a population analysis requires an appropriate sampling design.
