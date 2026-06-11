## Submission (v1.0.0)

This is an update to the package previously accepted as v0.9.0.

### R CMD check results (local, `--as-cran`)

0 errors | 0 warnings | 4 notes

### Notes

1. **Installed package size (9.2 MB):** The package bundles seven benchmark
   datasets (`benchmark_rlifting`, `benchmark_wavethresh`, `benchmark_adlift`,
   `benchmark_nlt`, and three irregular-grid variants) and eight rendered
   vignettes. Each dataset is the result of 1,000 Monte Carlo replications across
   multiple configurations and is necessary for the vignettes to reproduce their
   tables and figures without re-running multi-hour simulations. The compiled
   shared library accounts for ~4.6 MB (C++ core via Rcpp).

2. **Future file timestamps:** The local clock could not be verified against an
   external time source during the check. All source file timestamps are set by
   the development environment and do not reflect a genuine future date.

3. **Non-portable compiler flag `-mno-omit-leaf-frame-pointer`:** This flag is
   set by the system R installation (Ubuntu 24.04 LTS), not by the package.
   It does not affect portability or correctness on CRAN build machines.

4. **HTML validation skipped (`tidy` not found):** The `tidy` command-line tool
   is not installed in the local environment. HTML manual pages have been
   reviewed manually and contain no structural issues.

### Test environments

* Local Linux (Ubuntu 24.04 LTS), R 4.3.3, x86\_64
* win-builder (R-devel) — to be run before submission
* R-hub — to be run before submission

### Summary of changes since v0.9.0

* Native irregular-grid support in all three processing modes.
* Two new boundary extensions: `local_linear` and `one_sided`.
* SureShrink threshold rule (`threshold_method = "sure"`).
* SCAD shrinkage (`shrinkage = "scad"`).
* `tune_alpha_beta()` for automatic recursive-parameter selection.
* `diagnose_wavelet()` verification suite for custom wavelets.
* Seven benchmark datasets and eight vignettes bundled.
* `method` argument deprecated in favour of `threshold_method` + `shrinkage`
  (backward-compatible shim retained).
