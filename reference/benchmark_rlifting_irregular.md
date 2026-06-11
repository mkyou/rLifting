# rLifting Irregular-Grid Benchmark Results

MSE and timing statistics for rLifting denoising on seven irregular-grid
test signals (see
[`benchmark_adlift_irregular`](https://mkyou.github.io/rLifting/reference/benchmark_adlift_irregular.md)
for the signal set), covering all three modes (offline, causal, stream),
six built-in wavelets, five boundary modes, and twelve
threshold/shrinkage method combinations (universal/sure x
hard/soft/semisoft/scad, with and without `tune_alpha_beta` where
applicable). 7-quartile summaries over 1000 simulations per
configuration. Companion to
[`benchmark_adlift_irregular`](https://mkyou.github.io/rLifting/reference/benchmark_adlift_irregular.md)
and
[`benchmark_nlt_irregular`](https://mkyou.github.io/rLifting/reference/benchmark_nlt_irregular.md);
used by the irregular-grid benchmark vignette.

## Usage

``` r
data(benchmark_rlifting_irregular)
```

## Format

A data frame with 7560 rows and 47 columns:

- Signal:

  One of the seven irregular-grid test signals.

- Pkg:

  Always `"rLifting"`.

- Mode:

  `"offline"`, `"causal"`, or `"stream"`.

- Wavelet:

  Lifting scheme: `"haar"`, `"db2"`, `"cdf53"`, `"cdf97"`, `"dd4"`, or
  `"lazy"`.

- Boundary:

  Boundary extension mode: `"symmetric"`, `"periodic"`, `"zero"`,
  `"local_linear"`, or `"one_sided"`.

- Method:

  Composite label combining threshold rule, shrinkage, and tuned/untuned
  status (e.g. `"universal_tuned_soft"`).

- ThresholdMethod:

  Threshold rule: `"universal"` or `"sure"`.

- Shrinkage:

  Shrinkage rule: `"hard"`, `"soft"`, `"semisoft"`, or `"scad"`.

- AlphaUsed, BetaUsed:

  Threshold-recursion parameters used (post-tuning when applicable).

- NoiseSd:

  Per-signal Gaussian noise standard deviation.

- N:

  Number of simulations per configuration (1000).

- MSEpos_min, MSEpos_q1, MSEpos_median, MSEpos_mean, MSEpos_q3,
  MSEpos_max, MSEpos_se:

  MSE statistics with position-aware processing.

- MSEfix_min, MSEfix_q1, MSEfix_median, MSEfix_mean, MSEfix_q3,
  MSEfix_max, MSEfix_se:

  MSE statistics with position-ignoring processing (offline only; `NA`
  for causal/stream).

- Ratio_min, Ratio_q1, Ratio_median, Ratio_mean, Ratio_q3, Ratio_max,
  Ratio_se:

  Per-simulation ratio `MSEpos / MSEfix` summarised over the 1000
  simulations (offline only).

- Timepos_min, Timepos_q1, Timepos_median, Timepos_mean, Timepos_q3,
  Timepos_max, Timepos_se:

  Wall-time statistics with position-aware processing (seconds).

- Timefix_min, Timefix_q1, Timefix_median, Timefix_mean, Timefix_q3,
  Timefix_max, Timefix_se:

  Wall-time statistics with position-ignoring processing (offline only).

## Source

`data-raw/generate_rlifting_irregular_benchmark.R`

## Details

Each configuration is run twice in offline mode: with position-aware
processing (`t = t_phys` passed in, irregular path active) and
position-ignoring (`t = NULL`, uniform-grid treatment). The `Ratio_*`
columns report `MSEpos / MSEfix` per simulation as a direct measure of
the value of irregular handling for that configuration. Causal and
stream modes report only position-aware results
(`MSEpos_*`/`Timepos_*`); the position-ignoring columns are `NA` for
those rows.
