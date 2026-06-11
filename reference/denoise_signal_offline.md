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
  threshold_method = "universal",
  shrinkage = NULL,
  a = 3.7,
  method = NULL,
  extension = "symmetric",
  t = NULL,
  ll_k = 4L
)
```

## Arguments

- signal:

  Numeric vector containing the complete signal.

- scheme:

  A `lifting_scheme` object.

- alpha:

  Recursive threshold parameter (universal rule only). Ignored when
  `threshold_method = "sure"`.

- beta:

  Threshold scale factor (universal rule only). Ignored when
  `threshold_method = "sure"`.

- levels:

  Number of decomposition levels.

- threshold_method:

  Threshold-selection rule. One of `"universal"` (Donoho-Johnstone
  universal threshold with the recursive per-level decay parameterised
  by `alpha` and `beta`) or `"sure"` (SureShrink: per-level
  SURE-minimising threshold, capped at the universal value; `alpha` and
  `beta` are unused).

- shrinkage:

  Shrinkage rule applied above the threshold: `"hard"`, `"soft"`,
  `"semisoft"` (default), or `"scad"`.

- a:

  SCAD shape parameter (must be \> 2; default 3.7 per Fan & Li 2001).
  Used only when `shrinkage = "scad"`.

- method:

  Deprecated. Use `shrinkage` instead. If provided, takes precedence
  over `shrinkage` with a deprecation warning.

- extension:

  Extension mode: `"symmetric"`, `"periodic"`, `"zero"`,
  `"local_linear"`, or `"one_sided"`.

- t:

  Optional numeric vector of sample positions for irregular grids. Must
  be sorted and have the same length as `signal`. When supplied,
  irregular-grid Lagrange interpolation is applied in the predict steps
  (the scheme is validated by `.check_irregular_scheme`). Ignored by
  `extension = "one_sided"` (with a warning).

- ll_k:

  Local-linear neighbourhood size, used only when
  `extension = "local_linear"`. Default 4L; minimum 2; clamped to the
  signal length if larger.

## Value

Filtered numeric vector (same length as input).
