# Causal Benchmark Results

Comparison of execution time between `rLifting`'s optimized causal mode
and a naive sliding-window implementation using `wavethresh`.

## Usage

``` r
data(benchmark_causal)
```

## Format

A list containing:

- rLifting_Time_Avg:

  Average time (seconds) for rLifting.

- Wavethresh_Naive_Time:

  Time (seconds) for naive sliding window.

- Speedup_Factor:

  Ratio of Naive Time to rLifting Time.
