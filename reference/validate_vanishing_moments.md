# Validate Vanishing Moments

Verifies if the wavelet cancels polynomials of a specific degree.

## Usage

``` r
validate_vanishing_moments(scheme, degree = 0, tol = 1e-09)
```

## Arguments

- scheme:

  Object of class `lifting_scheme`.

- degree:

  Polynomial degree (0=Constant, 1=Ramp, 2=Parabola...).

- tol:

  Residual energy tolerance (default 1e-9).

## Value

List with status and residual energy.
