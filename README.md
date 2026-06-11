# rLifting — High-Performance Wavelet Lifting for R

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

`rLifting` is an R package for signal denoising via the Lifting Scheme (Sweldens, 1996). It provides three processing modes — **offline** (full signal), **causal** (sliding window, no look-ahead), and **stream** (sample-by-sample, persistent state) — under a unified API backed by a zero-allocation C++ engine.

---

## Why rLifting?

The established wavelet packages in R (`wavethresh`, `adlift`, `nlt`) are designed for offline use. `rLifting` adds two things none of them offer:

- **Causal and stream modes** with a ring-buffer architecture: fixed memory, O(T) total work, microsecond per-sample latency.
- **Native irregular-grid support**: pass a sorted `t` vector to any pipeline function and the predict step uses Lagrange interpolation at the actual sample positions — no separate regularisation step.

On regular grids, `rLifting` is ~40× faster per sample than `wavethresh` and ~40,000× faster than `adlift`/`nlt`, while trailing wavethresh by ~18% at the geometric mean MSE level. On irregular grids, it is ~20,000× faster than `adlift` with comparable accuracy (within ~10% on six of seven signals). Numbers come from `data(benchmark_rlifting)` and `data(benchmark_rlifting_irregular)` — 1,000 Monte Carlo replications per cell across the four Donoho–Johnstone test signals.

---

## Features

- Three processing modes with the same API: `denoise_signal_offline`, `denoise_signal_causal`, `new_wavelet_stream`.
- Six built-in wavelets: `haar`, `db2`, `cdf53`, `cdf97`, `dd4`, `lazy`. Extensible via `lift_step()` + `custom_wavelet()`.
- Four shrinkage methods: hard, soft, semisoft (default), SCAD (Antoniadis & Fan, 2001).
- Two threshold-selection rules: universal (VisuShrink, recursive α/β decay) and SURE (SureShrink, per-level Stein risk). `tune_alpha_beta()` minimises SURE to select α and β automatically.
- Five boundary extensions: `symmetric` (default), `periodic`, `zero`, `local_linear` (OLS extrapolation, neighbourhood `ll_k`), `one_sided` (renormalised filter at the edge).
- Irregular-grid support in every mode: pass `t` to batch functions, or `t_val` per sample to the stream closure.
- Diagnostic suite: `diagnose_wavelet()` verifies perfect reconstruction, vanishing moments, orthogonality, compact support, and shift sensitivity.

---

## Performance

Numbers from `data(benchmark_rlifting)` — regular grid, N = 1,000, medians over 1,000 Monte Carlo replications.

### Per-sample latency (offline mode)

| Wavelet | Median (µs) |
|:--------|------------:|
| haar    | 0.08        |
| cdf53   | 0.09        |
| db2     | 0.09        |
| dd4     | 0.10        |
| cdf97   | 0.11        |

Causal mode: ~10 µs/sample. Stream mode: ~17 µs/sample (~70% overhead over causal, attributable to the R closure call per sample).

### Speed vs other packages (offline, regular grid, per sample)

| Package    | Per-sample (µs) | vs rLifting |
|:-----------|----------------:|------------:|
| rLifting   | ~0.09           | —           |
| wavethresh | ~3.6            | ~40×        |
| adlift     | ~3,600          | ~40,000×    |
| nlt        | ~3,600          | ~40,000×    |

### Real-world signal (infant cardiac monitoring, N = 2,048)

From `vignette("v08-real-world")` — BabyECG from `wavethresh`:

| Mode              | Per-call (ms) | Per-sample (µs) |
|:------------------|:-------------:|:---------------:|
| Offline (regular) | 0.05          | 0.025           |
| Offline (irregular, 70% retained) | 0.10 | 0.071 |
| Causal (window 256) | 0.25        | 0.12            |
| Stream (full 2,048 samples) | 2.0 | ~10             |

On the held-out 30% of dropped samples (irregular grid), rLifting achieves **0.97 bpm RMSE** vs 1.03 bpm for interpolate-then-denoise and 1.73 bpm for raw interpolation.

---

## Installation

```r
# install.packages("remotes")
remotes::install_github("mkyou/rLifting")
```

---

## Quick start

### Offline denoising

```r
library(rLifting)

scheme = lifting_scheme("cdf53")

clean = denoise_signal_offline(
  noisy_signal, scheme,
  levels = 4,
  shrinkage   = "semisoft",
  extension   = "symmetric",
  alpha = 0.3, beta = 1.2
)
```

### Causal denoising (sliding window, no look-ahead)

```r
clean_causal = denoise_signal_causal(
  noisy_signal, scheme,
  window_size = 255, levels = 4,
  shrinkage   = "semisoft"
)
```

### Stream processing (sample by sample)

```r
processor = new_wavelet_stream(
  lifting_scheme("haar"),
  window_size = 255, levels = 4,
  shrinkage   = "semisoft", update_freq = 1
)

out = numeric(length(noisy_signal))
for (i in seq_along(noisy_signal)) {
  out[i] = processor(noisy_signal[i])
}
```

### Irregular grid

```r
# Batch (offline or causal): pass t alongside the signal
clean_irr = denoise_signal_offline(
  y_irr, lifting_scheme("cdf53"),
  t = t_phys, levels = 4,
  extension = "local_linear",
  threshold_method = "sure", shrinkage = "scad"
)

# Stream: pass t_val per sample
proc_irr = new_wavelet_stream(
  lifting_scheme("cdf53"), irregular = TRUE,
  window_size = 255, levels = 4
)
for (i in seq_along(y_irr)) {
  out[i] = proc_irr(y_irr[i], t_val = t_phys[i])
}
```

### Automatic parameter tuning

```r
tuned = tune_alpha_beta(signal, lifting_scheme("cdf53"), levels = 4)

clean = denoise_signal_offline(
  noisy_signal, lifting_scheme("cdf53"),
  levels = 4,
  alpha = tuned$alpha, beta = tuned$beta
)
```

---

## Documentation

Eight vignettes form a self-contained learning path, ordered from tour to reference.

| # | Identifier | Topic |
|:--|:-----------|:------|
| 1 | `v01-introduction` | Package overview; three modes; wavelet and threshold choice; irregular grids; quick decision guides |
| 2 | `v02-thresholding-and-tuning` | MAD noise estimate; universal vs SURE; α/β recursion; `tune_alpha_beta()` walkthrough; shrinkage choice |
| 3 | `v03-causal-stream` | Causality penalty; `window_size` and `update_freq` heuristics; per-sample latency by mode and wavelet |
| 4 | `v04-boundary-modes` | The five boundary extensions; empirical impact; wavelet-dependent recommendations |
| 5 | `v05-irregular-grids` | Lagrange interpolation in predict steps; wavelet eligibility; robust defaults from benchmark |
| 6 | `v06-extensions` | Built-in wavelet reference; custom wavelets; diagnostics; low-level pipeline |
| 7 | `v07-benchmarks` | Full empirical comparison vs `wavethresh`, `adlift`, `nlt`; speed and MSE across modes, signals, wavelets, boundaries |
| 8 | `v08-real-world` | End-to-end case study on infant cardiac monitoring (BabyECG); offline, irregular-grid, and stream modes |

```r
vignette("v01-introduction", package = "rLifting")
vignette("v02-thresholding-and-tuning", package = "rLifting")
vignette("v03-causal-stream", package = "rLifting")
vignette("v04-boundary-modes", package = "rLifting")
vignette("v05-irregular-grids", package = "rLifting")
vignette("v06-extensions", package = "rLifting")
vignette("v07-benchmarks", package = "rLifting")
vignette("v08-real-world", package = "rLifting")
```

### Design notes (`inst/notes/`)

Five implementation-reference documents covering the C++ internals, ordered from broad to deep:

| Doc | Scope |
|:----|:------|
| `00-design-overview` | Architectural index: package layers, central abstractions, four modes, pointers to subsystem docs |
| `01-lifting-scheme-and-transform` | `LiftingStep` struct; polyphase decomposition; predict/update math; irregular path with full Lagrange formula |
| `02-adaptive-thresholding` | MAD estimator; universal and SURE rules; α/β recursion; four shrinkage formulas; `tune_alpha_beta()` pipeline |
| `03-zero-allocation-engine` | Ring buffer mechanics; pre-allocated workspaces; XPtr finalizer; per-sample hot path |
| `04-boundary-and-threshold` | Four mandatory boundary code paths; local-linear OLS; `nth_element` MAD; shrinkage implementations |

---

## Roadmap

- **Multivariate denoising**: joint denoising of correlated signals incorporating covariance structure.
- **Coefficient access**: expose intermediate LWT coefficients from all pipeline functions for inspection and custom post-processing.

---

## License

MIT — see `LICENSE`.
