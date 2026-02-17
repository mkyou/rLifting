# Generate standard signals for wavelet testing (Donoho-Johnstone Benchmark)

Creates synthetic classic signals used for wavelet validation. Formulas
are based on Donoho and Johnstone and Liu et al.

## Usage

``` r
.generate_signal(type, n = 512)
```

## Arguments

- type:

  Signal type: "const", "ramp", "poly2", "poly3", "random", "impulse",
  "sine", "doppler", "heavisine", "bumps".

- n:

  Signal length (default 512).

## Value

Numeric vector.
