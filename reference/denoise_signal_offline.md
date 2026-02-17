# Offline Denoising (Global Batch)

Performs denoising on the entire signal at once using a non-causal
approach. Uses global statistics for recursive threshold calculation
(Eq. 9). This function is fully optimized in 'C++' (Zero-Allocation).

## Usage

``` r
denoise_signal_offline(
  signal,
  scheme,
  alpha = 0.3,
  beta = 1.2,
  levels = 3,
  method = "semisoft",
  extension = "symmetric"
)
```

## Arguments

- signal:

  Numeric vector containing the complete signal.

- scheme:

  A `lifting_scheme` object.

- alpha:

  Recursive threshold parameter.

- beta:

  Threshold scale factor.

- levels:

  Number of decomposition levels.

- method:

  Thresholding method ("hard", "soft", "semisoft").

- extension:

  Extension mode ("symmetric", "periodic", "zero").

## Value

Filtered numeric vector (same length as input).
