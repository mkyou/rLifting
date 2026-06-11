# Noisy Doppler Signal Example

Synthetic Doppler signal contaminated with Gaussian noise. Used in
[`vignette("v01-introduction")`](https://mkyou.github.io/rLifting/articles/v01-introduction.md)
and the boundary-mode comparison.

## Usage

``` r
data(doppler_example)
```

## Format

A data frame with 2048 rows and 3 columns:

- index:

  Time index (1..2048).

- original:

  The pure Doppler signal.

- noisy:

  The signal with added Gaussian noise (sd = 0.5).
