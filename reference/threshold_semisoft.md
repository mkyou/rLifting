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

## References

Liu, Z., Mi, Y., & Mao, Y. (2014). Improved real-time denoising method
based on lifting wavelet transform. *Measurement Science Review*, 14(3),
152–159.
[doi:10.2478/msr-2014-0020](https://doi.org/10.2478/msr-2014-0020)
