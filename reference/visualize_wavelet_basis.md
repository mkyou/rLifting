# Visualize Basis Functions (Scaling and Wavelet)

Plots the waveform by iterating the reconstruction over several levels.

## Usage

``` r
visualize_wavelet_basis(scheme, plot = TRUE, levels = 8)
```

## Arguments

- scheme:

  Object of class `lifting_scheme`.

- plot:

  Boolean.

- levels:

  Number of cascade levels.

## Value

Invisibly returns `NULL`. Called for side effects (plotting).
