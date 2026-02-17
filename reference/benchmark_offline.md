# Offline Benchmark Results

Comparison of execution time and reconstruction error (MSE) between
`rLifting` and other packages (wavethresh, wavelets) using Haar wavelet.

## Usage

``` r
data(benchmark_offline)
```

## Format

A data frame with the following columns:

- Pkg:

  Package name.

- Time:

  Execution time in seconds.

- MSE:

  Mean Squared Error.
