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
  threshold_method = "universal",
  shrinkage = NULL,
  a = 3.7,
  method = NULL,
  extension = "symmetric",
  update_freq = 1,
  irregular = FALSE,
  ll_k = 4L
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

  Threshold decay parameter (universal rule only). Ignored when
  `threshold_method = "sure"`.

- beta:

  Threshold gain factor (universal rule only). Ignored when
  `threshold_method = "sure"`.

- threshold_method:

  Threshold-selection rule. One of `"universal"` or `"sure"` (per-level
  SURE-minimising threshold, capped at the universal value; `alpha` and
  `beta` are unused).

- shrinkage:

  Shrinkage rule: `"hard"`, `"soft"`, `"semisoft"` (default), or
  `"scad"`.

- a:

  SCAD shape parameter (must be \> 2; default 3.7 per Fan & Li 2001).
  Used only when `shrinkage = "scad"`.

- method:

  Deprecated. Use `shrinkage` instead.

- extension:

  Boundary handling: `"symmetric"`, `"periodic"`, `"zero"`,
  `"local_linear"`, or `"one_sided"`.

- update_freq:

  How often to recompute threshold statistics (default 1). Set to `0` to
  freeze thresholds at the warm-up estimate (a warning is emitted).
  Negative values are rejected.

- irregular:

  Logical. If TRUE, the returned closure accepts a second argument
  `t_val` (the sample's time position) and applies position-aware
  interpolation in the predict steps.

- ll_k:

  Local-linear neighbourhood size, used only when
  `extension = "local_linear"`. Default 4L; minimum 2; clamped to
  `window_size` if larger.

## Value

A closure `processor(new_sample, t_val = NULL)` that accepts one sample
(and optionally its time position) and returns the filtered value.
