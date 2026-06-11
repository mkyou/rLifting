# General Thresholding Wrapper

General Thresholding Wrapper

## Usage

``` r
threshold(x, lambda, method = "soft", a = 3.7)
```

## Arguments

- x:

  Input vector.

- lambda:

  Threshold value.

- method:

  One of `"hard"`, `"soft"`, `"semisoft"`, or `"scad"`.

- a:

  SCAD shape parameter, ignored unless `method = "scad"`. Default 3.7
  (Fan-Li canonical).

## Value

Numeric vector of the same length as `x` with thresholded coefficients.
