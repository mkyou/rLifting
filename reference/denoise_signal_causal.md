# Causal Batch Denoising (Turbo Simulation)

Processes a complete signal simulating the sequential arrival of data.
Uses the specialized 'C++' class `WaveletEngine` to perform causal
filtering efficiently on a historical dataset.

## Usage

``` r
denoise_signal_causal(
  signal,
  scheme,
  levels = 1,
  window_size = 256,
  alpha = 0.3,
  beta = 1.2,
  method = "semisoft",
  extension = "symmetric",
  update_freq = 1
)
```

## Arguments

- signal:

  Complete vector of the noisy signal.

- scheme:

  `lifting_scheme` object.

- levels:

  Decomposition levels.

- window_size:

  Window size.

- alpha:

  Threshold decay parameter (Eq 9).

- beta:

  Threshold gain factor (Eq 9).

- method:

  Thresholding method ("soft", "hard", "semisoft").

- extension:

  Boundary treatment ('symmetric', 'periodic').

- update_freq:

  Frequency of threshold updates.

## Value

Filtered vector (same length as input).
