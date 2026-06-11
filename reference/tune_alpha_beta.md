# Tune Adaptive Threshold Parameters via SURE

Selects the recursive-threshold parameters `alpha` and `beta` that
minimise Stein's Unbiased Risk Estimate (SURE) for soft thresholding
applied to the wavelet detail coefficients of the supplied signal.

## Usage

``` r
tune_alpha_beta(
  signal,
  scheme,
  levels = 3,
  extension = "symmetric",
  ll_k = 4L,
  alpha_range = c(0, 10),
  beta_range = c(0.5, 3)
)
```

## Arguments

- signal:

  Numeric vector.

- scheme:

  A `lifting_scheme` object.

- levels:

  Decomposition depth.

- extension:

  Boundary mode (passed to `lwt`).

- ll_k:

  Local-linear neighborhood size (passed to `lwt`).

- alpha_range:

  Bounds for `alpha`; default `c(0, 10)`.

- beta_range:

  Bounds for `beta`; default `c(0.5, 3.0)`.

## Value

A list with components `alpha`, `beta`, `sure` (the minimised SURE
value), and `converged` (logical).

## Details

SURE is computed under the soft-threshold estimator assumption (Donoho &
Johnstone, 1995). The chosen parameters are typically usable with
`shrinkage = "soft"`, `"semisoft"`, `"hard"`, or `"scad"`, since all
four share the same threshold location.
