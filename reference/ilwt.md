# Inverse Lifting Wavelet Transform ('C++' Accelerated)

Reconstructs the original signal from wavelet coefficients. Optimized
with 'C++' backend.

## Usage

``` r
ilwt(lwt_obj, scheme = NULL)
```

## Arguments

- lwt_obj:

  Object of class 'lwt' returned by [`lwt()`](lwt.md).

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
print(rec) # Should match s
#> [1] 1 2 3 4
```
