# Validate Perfect Reconstruction (Stress Test)

Verifies wavelet invertibility against a battery of signals.

## Usage

``` r
validate_perfect_reconstruction(scheme, tol = 1e-09)
```

## Arguments

- scheme:

  Object of class `lifting_scheme`.

- tol:

  Numerical error tolerance (default 1e-9).

## Value

List with global status and maximum error found.
