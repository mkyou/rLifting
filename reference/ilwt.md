# Inverse Lifting Wavelet Transform ('C++' Accelerated)

Reconstructs the original signal from wavelet coefficients. Optimized
with 'C++' backend.

## Usage

``` r
ilwt(lwt_obj, scheme = NULL)
```

## Arguments

- lwt_obj:

  Object of class `lwt` returned by
  [`lwt()`](https://mkyou.github.io/rLifting/reference/lwt.md). The
  fields `extension`, `ll_k`, and `t` carried by the object are reused
  to mirror the forward transform; the inverse cannot be invoked with a
  different boundary mode or grid.

- scheme:

  (Optional) `lifting_scheme` object. If NULL, uses the one from
  `lwt_obj`.

## Value

Numeric vector containing the reconstructed signal.

## Examples

``` r
s = c(1, 2, 3, 4)
sch = lifting_scheme("haar")
fwd = lwt(s, sch)
#> Warning: Residual signal at level 1 has only 2.0 samples.
rec = ilwt(fwd)
print(rec)
#> [1] 1 2 3 4
```
