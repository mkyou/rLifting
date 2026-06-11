# wavethresh Benchmark Results

MSE and timing statistics for the wavethresh package on the four
Donoho-Johnstone signals, across its native wavelets and
boundary/threshold combinations. 7-quartile summaries over 1000
simulations per configuration. Used as a reference baseline in the
offline benchmark vignette.

## Usage

``` r
data(benchmark_wavethresh)
```

## Format

A data frame with rows per (Signal, Wavelet, Boundary) and 19 columns:

- Signal:

  Donoho-Johnstone test signal.

- Pkg:

  Always `"wavethresh"`.

- Wavelet:

  Daubechies filter family identifier (e.g. `"co1"` for
  Daubechies-extremal-phase order 1).

- Boundary:

  Combination of wavethresh threshold policy and shrinkage rule (e.g.
  `"BayesThresh_soft"`, `"cv_hard"`).

- N:

  Number of simulations per configuration (1000).

- MSE_min, MSE_q1, MSE_median, MSE_mean, MSE_q3, MSE_max, MSE_se:

  Reconstruction MSE statistics.

- Time_min, Time_q1, Time_median, Time_mean, Time_q3, Time_max, Time_se:

  Wall-time statistics (seconds).

## Source

`data-raw/generate_wavethresh_benchmark.R`
