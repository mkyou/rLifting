# rLifting 1.0.0

## New features

*   **Native irregular-grid support**: pass a sorted `t` vector to
    `denoise_signal_offline()`, `denoise_signal_causal()`, or `t_val` per sample
    to the stream closure from `new_wavelet_stream()`. The predict step uses
    Lagrange interpolation conditioned on the actual sample positions — no
    separate regularisation required.
*   **Two new boundary extensions**: `local_linear` (OLS extrapolation with
    configurable neighbourhood `ll_k`) and `one_sided` (renormalised filter at the
    edge), bringing the total to five alongside the existing `symmetric`, `periodic`,
    and `zero`.
*   **SureShrink threshold rule** (`threshold_method = "sure"`): per-level Stein
    Unbiased Risk Estimator, complementing the existing universal (VisuShrink) rule.
*   **SCAD shrinkage** (`shrinkage = "scad"`): Smoothly Clipped Absolute Deviation
    (Antoniadis & Fan, 2001), added alongside the existing hard, soft, and semisoft
    shrinkages.
*   **`tune_alpha_beta()`**: automatic selection of the recursive threshold
    parameters α and β by minimising SURE over the signal.
*   **`diagnose_wavelet()`**: verification suite for custom wavelets — perfect
    reconstruction, vanishing moments, orthogonality, compact support, and shift
    sensitivity.
*   Seven benchmark datasets bundled: `benchmark_rlifting`, `benchmark_wavethresh`,
    `benchmark_adlift`, `benchmark_nlt`, `benchmark_rlifting_irregular`,
    `benchmark_adlift_irregular`, `benchmark_nlt_irregular` — each based on 1,000
    Monte Carlo replications per configuration.
*   Eight vignettes: introduction, thresholding and tuning, causal/stream,
    boundary modes, irregular grids, extensions, full cross-package benchmarks,
    and a real-world case study (BabyECG infant cardiac monitoring).

## Deprecated

*   The `method` argument in `denoise_signal_offline()` and
    `denoise_signal_causal()` is deprecated in favour of the separate
    `threshold_method` and `shrinkage` arguments. Old code continues to work
    via a backward-compatibility shim; `method` will be removed in a future
    version.

# rLifting 0.9.0

*   Initial release of the package.
*   Implements high-performance Wavelet Lifting Transforms.
*   Features:
    *   Unified offline (batch) denoising via `denoise_signal_offline()`.
    *   Causal (sliding-window) filtering via `denoise_signal_causal()`.
    *   Sample-by-sample stream processing via `new_wavelet_stream()`.
    *   Adaptive recursive thresholding (universal rule, hard/soft/semisoft shrinkage).
    *   Zero-allocation 'C++' core via 'Rcpp'.
    *   Six built-in wavelets: Haar, DB2, CDF 5/3, CDF 9/7, DD4, Lazy.
    *   Three boundary extensions: `symmetric`, `periodic`, `zero`.
