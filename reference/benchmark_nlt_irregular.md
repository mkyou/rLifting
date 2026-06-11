# nlt Irregular-Grid Benchmark Results

MSE and timing statistics for the nlt package on seven irregular-grid
test signals (see
[`benchmark_adlift_irregular`](https://mkyou.github.io/rLifting/reference/benchmark_adlift_irregular.md)
for the signal set). 7-quartile summaries over 1000 simulations per
configuration. Used by the irregular-grid benchmark vignette as a
reference baseline.

## Usage

``` r
data(benchmark_nlt_irregular)
```

## Format

A data frame with 672 rows and 20 columns:

- Signal:

  One of the seven irregular-grid test signals.

- Pkg:

  Always `"nlt"`.

- Wavelet:

  Predictor family inherited from adlift: `"AdaptPred"`, `"CubicPred"`,
  `"LinearPred"`, or `"QuadPred"`.

- Boundary:

  Encoded combination of nlt/adlift options.

- NoiseSd:

  Per-signal Gaussian noise standard deviation.

- N:

  Number of simulations per configuration (1000).

- MSE_min, MSE_q1, MSE_median, MSE_mean, MSE_q3, MSE_max, MSE_se:

  Reconstruction MSE statistics.

- Time_min, Time_q1, Time_median, Time_mean, Time_q3, Time_max, Time_se:

  Wall-time statistics (seconds).

## Source

`data-raw/generate_nlt_irregular_benchmark.R`
