# Phase C structural baseline

This is the historical 688-model snapshot from nlpsem 1.0.0.9002, whose fitting
source is preserved at commit 31f3d61bb256b12b9ec1b345dc364c8680925fe6.
The snapshot was recompressed without changing its deserialized contents.
`provenance.json` records both hashes, the source revision, generator hash, and
baseline environment.

From the repository root:

```sh
Rscript -e 'source("data-raw/snapshot_models.R")' . test-results/model-structures.rds
Rscript data-raw/test_structure_comparator.R
Rscript data-raw/compare_model_structures.R data-raw/structural-baseline/phase-c-1.0.0.9002.rds test-results/model-structures.rds
```

Create `test-results/` before running these commands. The GitHub full fitting
workflow runs them automatically and retains their logs.

The 1.0.0.9003 comparison explicitly permits the reviewed split of a shared
ADD1 matrix into outcome-specific, identical constants for the equal-wave cases
in this historical grid. Their dimensions, values, free masks, labels, and
per-outcome algebra references remain checked exactly.

To accommodate numerical roundoff between macOS arm64 and hosted Linux, numeric
matrix values, lower/upper bounds, and the data checksum may differ elementwise by
at most `1e-10 * max(1, abs(baseline), abs(current))`. This combines a small absolute
allowance near zero with a relative allowance for larger values, including summed
data checksums. Numeric types, dimensions and other attributes, NA/NaN/Inf values,
free masks, labels, formulas, and all fields outside that explicit list remain
exact. ADD1 constants remain exact even within the numerical allowance. All other
differences fail. This allowance does not alter numerical fitting-reference
tolerances.

The comparator requires all 688 model identifiers and reports four separate
counts: exact matches, approved ADD1 renames only, roundoff only, and approved
renames plus roundoff. It also reports the largest accepted absolute and scaled
numeric discrepancy. Its self-tests reject missing or misreferenced constants,
meaningful value changes, type/dimension mismatches, and unrelated mutations.

The snapshot records historical behavior; it is not an independent statistical
reference. The new review regression tests separately cover distinct mediation
times, zero covariance starts, three/four outcomes, unequal/reordered schedules,
and numerical inference. Those tests assert mathematical properties and expected
values, including cases absent from the historical grid. Do not replace the
baseline simply to hide a failed comparison.
