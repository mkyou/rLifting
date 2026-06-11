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
  threshold_method = "universal",
  shrinkage = NULL,
  a = 3.7,
  method = NULL,
  extension = "symmetric",
  update_freq = 1,
  t = NULL,
  ll_k = 4L
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

  Boundary treatment: `"symmetric"`, `"periodic"`, `"zero"`,
  `"local_linear"`, or `"one_sided"`.

- update_freq:

  Frequency of threshold updates. Set to `0` to freeze thresholds at the
  warm-up estimate (a warning is emitted). Negative values are rejected.

- t:

  Optional numeric vector of sample time positions (irregular grid).
  Must be sorted and the same length as `signal`. Ignored by
  `extension = "one_sided"` (with a warning).

- ll_k:

  Local-linear neighbourhood size, used only when
  `extension = "local_linear"`. Default 4L; minimum 2; clamped to
  `window_size` if larger.

## Value

Filtered vector (same length as input).
