# Calculate Adaptive Threshold (Recursive)

Estimates the optimal noise threshold based on current window
statistics. Implements the recursive formula from Liu et al. (2014).
Accelerated with 'C++'.

## Usage

``` r
compute_adaptive_threshold(lwt_obj, alpha = 0.3, beta = 1.2)
```

## Arguments

- lwt_obj:

  Object returned by [`lwt()`](lwt.md).

- alpha:

  Recursive adjustment parameter (Eq. 9).

- beta:

  Initial threshold scale factor (Eq. 9).

## Value

Object of class `adaptive_thresholds` (a list of thresholds).
