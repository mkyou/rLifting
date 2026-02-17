# Validate Orthogonality (Energy Conservation)

Verifies Parseval's Theorem. Only applicable for orthogonal wavelets.

## Usage

``` r
validate_orthogonality(scheme, expected = TRUE, tol = 1e-09)
```

## Arguments

- scheme:

  Object of class `lifting_scheme`.

- expected:

  Boolean. If TRUE, expects orthogonality.

- tol:

  Tolerance (default 1e-9).

## Value

List with status and energy ratio (Out/In).
