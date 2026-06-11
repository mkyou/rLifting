# adlift Irregular-Grid Benchmark Results

MSE and timing statistics for the adlift package on seven irregular-grid
test signals: three physically motivated (`linear_phys`, `trend_events`,
`blocks_gapped`) and four Donoho-Johnstone classics re-sampled on
irregular grids (`blocks_dj_irr`, `bumps_dj_irr`, `doppler_dj_irr`,
`heavisine_dj_irr`). 7-quartile summaries over 1000 simulations per
configuration. Companion to
[`benchmark_nlt_irregular`](https://mkyou.github.io/rLifting/reference/benchmark_nlt_irregular.md)
and
[`benchmark_rlifting_irregular`](https://mkyou.github.io/rLifting/reference/benchmark_rlifting_irregular.md);
used by the irregular-grid benchmark vignette.

## Usage

``` r
data(benchmark_adlift_irregular)
```

## Format

A data frame with 672 rows and 20 columns:

- Signal:

  One of the seven irregular-grid test signals (see Description).

- Pkg:

  Always `"adlift"`.

- Wavelet:

  Predictor family: `"AdaptPred"`, `"CubicPred"`, `"LinearPred"`, or
  `"QuadPred"`.

- Boundary:

  Encoded combination of `adlift::fwtnp` options.

- NoiseSd:

  Per-signal Gaussian noise standard deviation (e.g. 0.15 for
  `linear_phys`/`trend_events`, 0.30 for the DJ-irregular signals, 0.50
  for `blocks_gapped`).

- N:

  Number of simulations per configuration (1000).

- MSE_min, MSE_q1, MSE_median, MSE_mean, MSE_q3, MSE_max, MSE_se:

  Reconstruction MSE statistics.

- Time_min, Time_q1, Time_median, Time_mean, Time_q3, Time_max, Time_se:

  Wall-time statistics (seconds).

## Source

`data-raw/generate_adlift_irregular_benchmark.R`
