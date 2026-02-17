# Semisoft Shrinkage (Hyperbolic)

Implementation based on Liu et al. (2014). Combines the stability of
Soft Thresholding with the amplitude precision of Hard Thresholding.
Function: `sign(x) * sqrt(x^2 - lambda^2)` for values above lambda.

## Usage

``` r
threshold_semisoft(x, lambda)
```

## Arguments

- x:

  Vector of coefficients.

- lambda:

  Positive threshold value.

## Value

Processed vector.
