# Noisy Doppler Signal Example

A synthetic dataset containing a Doppler signal contaminated with
Gaussian noise. Used in the "General Usage" vignette.

## Usage

``` r
data(doppler_example)
```

## Format

A data frame with 2048 rows and 3 columns:

- index:

  Time index.

- original:

  The pure Doppler signal.

- noisy:

  The signal with added Gaussian noise (sd=0.5).
