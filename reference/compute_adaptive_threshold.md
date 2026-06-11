# Calculate Adaptive Threshold (Universal / Recursive)

Estimates the per-level noise threshold from the finest-level detail
coefficients and applies the recursive Liu et al. (2014) decay across
levels. This is the step-by-step entry point for the universal threshold
rule.

## Usage

``` r
compute_adaptive_threshold(lwt_obj, alpha = 0.3, beta = 1.2)
```

## Arguments

- lwt_obj:

  Object returned by
  [`lwt()`](https://mkyou.github.io/rLifting/reference/lwt.md).

- alpha:

  Recursive adjustment parameter (Eq. 9 of Liu et al., 2014).

- beta:

  Initial threshold scale factor (Eq. 9 of Liu et al., 2014).

## Value

Object of class `adaptive_thresholds` (a list of thresholds).

## Details

To use SureShrink instead, call
[`denoise_signal_offline()`](https://mkyou.github.io/rLifting/reference/denoise_signal_offline.md)
or the causal/stream functions with `threshold_method = "sure"` — the
SURE branch lives in the C++ engine and is not exposed as a standalone R
routine. For automatic selection of `alpha` and `beta`, see
[`tune_alpha_beta`](https://mkyou.github.io/rLifting/reference/tune_alpha_beta.md).

## References

Donoho, D. L., & Johnstone, I. M. (1994). Ideal spatial adaptation by
wavelet shrinkage. *Biometrika*, 81(3), 425–455.

Liu, Z., Mi, Y., & Mao, Y. (2014). Improved real-time denoising method
based on lifting wavelet transform. *Measurement Science Review*, 14(3),
152–159.
[doi:10.2478/msr-2014-0020](https://doi.org/10.2478/msr-2014-0020)
