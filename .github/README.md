# Package automation

These workflows are configuration for a future GitHub run. Adding them locally does
not establish that Windows, Linux, macOS, or hosted documentation checks have passed.

| Workflow | Trigger | Work performed |
|---|---|---|
| `R-CMD-check.yaml` | Pull requests, pushes to `main`/`master`, or manual dispatch | R-release package checks on Windows, Linux, and macOS, including normal tests, examples, and vignette checks. `NOT_CRAN=false` excludes tests marked `skip_on_cran()`. Warnings fail the job; check logs are retained. |
| `full-fitting-tests.yaml` | Manual dispatch only | The complete test suite on Linux R release with `NOT_CRAN=true`, including the slower fitting references. Console evidence is saved as an artifact, including on failure. |
| `pkgdown.yaml` | Pull requests, pushes to `main`/`master`, or manual dispatch | Builds reference pages and articles into a downloadable documentation artifact. No deployment step is configured. |

All workflows request only `contents: read`. Checkout credentials are not retained.
Action major versions follow the official
[r-lib examples](https://github.com/r-lib/actions/tree/v2/examples) and
[upload-artifact usage](https://github.com/actions/upload-artifact#usage).

## Local checks before publishing these workflows

1. Parse the YAML and validate GitHub Actions syntax with `actionlint` when available.
2. Run `R CMD build` and `R CMD check --as-cran --no-manual` on the source archive.
   Confirm that `.github`, `_pkgdown.yml`, generated `docs`, and the redundant full
   license text in `LICENSE.md` are excluded from the R source archive.
3. Run the complete test suite explicitly:

   ```r
   Sys.setenv(NOT_CRAN = "true")
   testthat::test_local(".", reporter = "summary", stop_on_failure = TRUE)
   ```

4. Generate the documentation after the Rd files have been regenerated:

   ```r
   pkgdown::build_site(preview = FALSE, examples = FALSE, new_process = TRUE)
   ```

   The development site is written under `docs/dev/`. Reference examples remain
   visible, while the package check executes examples and the site build executes
   article chunks. Inspect the home, reference, and fitted-object workflow pages.
5. Once the workflows are present in a GitHub repository, inspect each platform's
   actual run and manually dispatch the long suite. Keep hosted results distinct
   from local results; do not infer them from configuration alone.

The `License: GPL (>= 3)` declaration in `DESCRIPTION` remains authoritative.
`LICENSE.md` makes the existing license readable in the repository and documentation;
R already ships the standard GPL text, so the redundant copy is excluded from the
package archive in accordance with the
[R package licensing guidance](https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Licensing).
