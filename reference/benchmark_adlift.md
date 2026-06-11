# adlift Benchmark Results

MSE and timing statistics for the adlift package (adaptive lifting on
irregular grids) on the four Donoho-Johnstone signals. 7-quartile
summaries over 1000 simulations per configuration. Used as a reference
baseline in the irregular-grid benchmark vignette.

## Usage

``` r
data(benchmark_adlift)
```

## Format

A data frame with 384 rows and 19 columns:

- Signal:

  Donoho-Johnstone test signal.

- Pkg:

  Always `"adlift"`.

- Wavelet:

  Predictor family: `"AdaptPred"`, `"CubicPred"`, `"LinearPred"`, or
  `"QuadPred"`.

- Boundary:

  Encoded combination of `adlift::fwtnp` options (neighbours,
  interpolation, closest-point, mean/median predictor).

- N:

  Number of simulations per configuration (1000).

- MSE_min, MSE_q1, MSE_median, MSE_mean, MSE_q3, MSE_max, MSE_se:

  Reconstruction MSE statistics.

- Time_min, Time_q1, Time_median, Time_mean, Time_q3, Time_max, Time_se:

  Wall-time statistics (seconds).

## Source

`data-raw/generate_adlift_benchmark.R`
