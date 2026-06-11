# SCAD Shrinkage (Antoniadis & Fan, 2001)

Smoothly Clipped Absolute Deviation shrinkage. Three-region rule with
continuity at lambda, 2*lambda and a*lambda. Identity above a\*lambda
eliminates the bias that soft thresholding imposes on large
coefficients, while keeping sparsity in the zero region.

## Usage

``` r
threshold_scad(x, lambda, a = 3.7)
```

## Arguments

- x:

  Vector of coefficients.

- lambda:

  Positive threshold value.

- a:

  Shape parameter, a \> 2. Default 3.7 (Fan-Li canonical).

## Value

Processed vector.
