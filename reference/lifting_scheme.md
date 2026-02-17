# Lifting Scheme Constructor

Creates an S3 object containing the prediction (P) and update (U) steps
required for the Lifting Transform.

## Usage

``` r
lifting_scheme(wavelet = "haar", custom_steps = NULL, custom_norm = NULL)
```

## Arguments

- wavelet:

  Wavelet name (string). Options: "haar", "db2", "cdf53", "cdf97",
  "dd4", "lazy".

- custom_steps:

  List of custom steps (optional). If provided, ignores internal lookup.

- custom_norm:

  Normalization vector (optional).

## Value

An object of class `lifting_scheme`.
