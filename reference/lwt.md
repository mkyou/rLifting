# Lifting Wavelet Transform (Forward)

Performs the Forward Wavelet Transform using the Lifting Scheme.
Optimized with 'C++' backend.

## Usage

``` r
lwt(signal, scheme, levels = 1, extension = "symmetric", t = NULL, ll_k = 4L)
```

## Arguments

- signal:

  Numeric vector containing the input signal.

- scheme:

  A `lifting_scheme` object.

- levels:

  Integer. Number of decomposition levels.

- extension:

  Boundary extension mode: `"symmetric"` (default), `"periodic"`,
  `"zero"`, `"local_linear"` (linear extrapolation from boundary
  samples), or `"one_sided"` (asymmetric filter renormalisation at the
  boundary).

- t:

  Optional numeric vector of sample positions for irregular grids. Must
  be sorted and have the same length as `signal`. When supplied,
  irregular-grid Lagrange interpolation is applied in the predict steps
  and `lwt_obj$t` is stored for use by
  [`ilwt()`](https://mkyou.github.io/rLifting/reference/ilwt.md).
  Ignored by `extension = "one_sided"` (with a warning).

- ll_k:

  Local-linear neighbourhood size, used only when
  `extension = "local_linear"`. Default 4L; minimum 2; clamped to the
  signal length if larger.

## Value

An object of class `lwt`. It is a list containing `coeffs` (list of
details d1..dn and approximation an), `scheme` (the scheme object used),
`levels`, `original_len`, `extension`, `ll_k`, and `t`.

## Examples

``` r
data = c(1, 2, 3, 4, 5, 6, 7, 8)
sch = lifting_scheme("haar")
res = lwt(data, sch, levels = 2)
#> Warning: Residual signal at level 2 has only 2.0 samples.
print(res)
#> --- LWT Decomposition (C++ Accelerated) ---
#> Levels: 2
#> Wavelet: haar
#> Coefficients:
#>   a2: length 2
#>   d1: length 4
#>   d2: length 2
```
