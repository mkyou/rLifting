# Create an Adaptive Wavelet Stream Processor ('C++' Core)

Generates a stateful function backed by a high-performance 'C++' Ring
Buffer engine. It implements Sliding Window + Lifting Decomposition +
Adaptive Thresholding in highly efficient time per sample.

## Usage

``` r
new_wavelet_stream(
  scheme,
  window_size = 256,
  levels = 1,
  alpha = 0.3,
  beta = 1.2,
  method = "semisoft",
  extension = "symmetric",
  update_freq = 1
)
```

## Arguments

- scheme:

  A `lifting_scheme` object.

- window_size:

  Sliding window size (W). Must be \> 8.

- levels:

  Decomposition levels (default 1).

- alpha:

  Threshold decay parameter (Eq 9).

- beta:

  Threshold gain factor (Eq 9).

- method:

  Shrinkage method: "hard", "soft", "semisoft".

- extension:

  Boundary handling ('symmetric', 'periodic', 'zero').

- update_freq:

  How often to recompute threshold statistics (default 1).

## Value

A closure function `processor(new_sample)` that accepts a single numeric
value and returns the filtered value immediately.
