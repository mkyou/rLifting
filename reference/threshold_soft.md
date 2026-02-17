# Soft Thresholding

Sets coefficients below the threshold to zero and shrinks others towards
zero. Reduces noise variance but introduces amplitude bias.

## Usage

``` r
threshold_soft(x, lambda)
```

## Arguments

- x:

  Vector of coefficients.

- lambda:

  Positive threshold value.

## Value

Processed vector.
