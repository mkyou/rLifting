# rLifting Offline/Causal Benchmark Results

MSE and timing statistics for rLifting denoising across the four
Donoho-Johnstone signals, three modes (offline / causal-batch / stream),
multiple wavelets, boundary modes, threshold rules, and shrinkage
methods. Reported as 7-quartile summaries (min / q1 / median / mean / q3
/ max / se) over 1000 simulations per configuration. Used by
[`vignette("v02-thresholding-and-tuning")`](https://mkyou.github.io/rLifting/articles/v02-thresholding-and-tuning.md)
and
[`vignette("v03-causal-stream")`](https://mkyou.github.io/rLifting/articles/v03-causal-stream.md).

## Usage

``` r
data(benchmark_rlifting)
```

## Format

A data frame with 3600 rows and 40 columns:

- Signal:

  Donoho-Johnstone test signal: `"blocks"`, `"bumps"`, `"doppler"`, or
  `"heavisine"`.

- Pkg:

  Always `"rLifting"`.

- Mode:

  `"offline"`, `"causal"`, or `"stream"`.

- Wavelet:

  Lifting scheme: `"haar"`, `"cdf53"`, etc.

- Boundary:

  Boundary extension mode.

- Method:

  Composite label combining threshold rule and shrinkage.

- ThresholdMethod:

  Threshold rule: `"universal"` or `"sure"`.

- Shrinkage:

  Shrinkage rule: `"hard"`, `"soft"`, `"semisoft"`, or `"scad"`.

- AlphaUsed, BetaUsed:

  Threshold-recursion parameters used.

- Version:

  Generator version tag (`"v2"` for current).

- N:

  Number of simulations per configuration (1000).

- MSE_min, MSE_q1, MSE_median, MSE_mean, MSE_q3, MSE_max, MSE_se:

  Mean-squared-error statistics against the noise-free signal.

- MSE_settled_min, MSE_settled_q1, MSE_settled_median, MSE_settled_mean,
  MSE_settled_q3, MSE_settled_max, MSE_settled_se:

  MSE statistics computed after dropping the warm-up window
  (causal/stream modes only).

- Time_total_min, Time_total_q1, Time_total_median, Time_total_mean,
  Time_total_q3, Time_total_max, Time_total_se:

  Total wall time per call (seconds).

- Per_sample_us_min, Per_sample_us_q1, Per_sample_us_median,
  Per_sample_us_mean, Per_sample_us_q3, Per_sample_us_max,
  Per_sample_us_se:

  Per-sample time (microseconds), i.e. Time_total / signal length.

## Source

`data-raw/generate_rlifting_benchmark_v2.R`
