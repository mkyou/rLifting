# rLifting — Design Overview

This document describes the design decisions of the `rLifting` package, from the broadest architectural choices down to specific implementation details. Each section may expand into its own document as depth grows.

---

## 1. What rLifting Is

`rLifting` is an R package for signal filtering and denoising via the Lifting Wavelet Transform (LWT). Its main differentiators from packages such as `wavethresh` and `adlift` are:

- **Speed**: C++ core (Rcpp) with no allocations in the hot path;
- **Causal / stream mode**: sample-by-sample filter with a sliding window, suitable for real-time processing;
- **Irregular grids**: native support for non-equispaced sampling — Lagrange interpolation is embedded inside the predict steps, so no separate grid-regularisation pass is needed;
- **Extensibility**: any lifting scheme can be defined in R and executed by the C++ engine.

---

## 2. Two-Layer Architecture

```
User (R)
    │
    ├── lifting_scheme / lift_step / custom_wavelet   ← defines the wavelet lifting scheme
    │
    ├── lwt / ilwt                                    ← step-by-step decomposition / reconstruction
    ├── compute_adaptative_threshold / threshold      ← step-by-step thresholding
    │
    ├── denoise_signal_offline                        ← batch denoising (non-causal)
    ├── denoise_signal_causal                         ← causal denoising (batch over history)
    └── new_wavelet_stream                            ← sample-by-sample processor (closure)
            │
            ▼
    C++ (Rcpp)
    ├── lwt_cpp / ilwt_cpp                            ← multi-level LWT
    ├── denoise_offline_cpp                           ← LWT + threshold + ILWT in one pass
    ├── WaveletEngine (XPtr)                          ← ring buffer + LWT + threshold + ILWT
    └── utilities: apply_filter_cpp, get_val_safe,
                   interp_predict, get_t_extrap, onesided_conv
```

The R layer validates inputs, converts types, and dispatches to C++. The C++ layer pre-allocates all working buffers at construction time and performs no heap allocations in the hot path.

---

## 3. The Lifting Scheme as Central Abstraction

Every public function takes a `lifting_scheme` object, which encapsulates:

- `steps`: a sequence of predict (P) and update (U) steps, each with `type`, `coeffs`, `start_idx`, and `degree`;
- `normalization`: vector `[norm_approx, norm_detail]`;
- `wavelet`: name of the wavelet.

Built-in wavelets: `haar`, `db2`, `cdf53`, `cdf97`, `dd4`, `lazy`. Users may define custom wavelets via `custom_wavelet()`.

The `degree` field governs irregular-grid behaviour and is inferred automatically in `lifting_scheme()` based on the sum of the predict step coefficients.

→ *See: `04-irregular-grid-design.md` (section 5 — wavelet classification and degree inference)*

---

## 4. Modes of Operation

### 4.1 Step-by-step (manual pipeline)

The user may assemble the pipeline manually: `lwt` → `compute_adaptive_threshold` → `threshold` → `ilwt`. Each function dispatches to its own C++ routine. This path is useful for inspecting intermediate coefficients, applying custom thresholding logic, or experimenting with the transform independently of denoising. It is batch-only — causal processing requires a sliding window that must be managed manually, which is what sections 4.3–4.4 handle internally.

### 4.2 Offline (non-causal)

`denoise_signal_offline` has access to the full signal. Both regular and irregular modes route to `denoise_offline_cpp` — a single C++ pass (LWT + threshold + ILWT) with no R overhead. The irregular path stores per-level `t` position vectors and applies Lagrange interpolation in predict steps with `degree >= 0`; it is approximately 1.7× slower than the regular path at n = 1024 due to the interpolation cost, but 4.4× faster than the naive R-level routing (`lwt → threshold → ilwt`).

### 4.3 Causal Batch

`denoise_signal_causal` simulates sequential sample arrival over a historical signal. Calls `run_causal_batch_cpp`, which instantiates a `WaveletEngine` and processes samples one by one.

Both regular and irregular modes run **entirely in C++** — `WaveletEngine::push_and_process` accepts a `t_val` argument and performs position-aware interpolation natively.

### 4.4 Stream (sample-by-sample)

`new_wavelet_stream` returns an R closure encapsulating a `WaveletEngine` via `Rcpp::XPtr`. Each call to the closure pushes one sample into the ring buffer and returns the filtered value. No R-level allocations occur in the hot path.

Both regular and irregular modes run **entirely in C++**.

→ *See: `02-modes-of-operation.md` (full detail on all modes)*

---

## 5. Boundary Handling

Five extension modes, passed as an integer from R to C++:

| Integer | Name | Behaviour |
|:--------|:-----|:----------|
| 1 | `symmetric` | Mirror reflection (default) |
| 2 | `periodic` | Periodisation |
| 3 | `zero` | Zero-padding |
| 4 | `local_linear` | Local linear extrapolation |
| 5 | `one_sided` | Renormalised one-sided filter |

Adding or modifying a mode requires changes in **four C++ locations**: `utils.h` (`get_val_safe`), `utils.cpp` (`apply_filter_cpp`), `offline.cpp`, and `WaveletEngine.h`. Mode 5 is the only one that requires additional logic beyond `get_val_safe`, since it changes filter normalisation rather than just the value at a single out-of-bounds index.

→ *See: `03-boundary-modes.md` (full detail on all modes, wavelet interaction, and causal implications)*

---

## 6. Irregular Grids

Signals with physically non-equispaced time positions are supported by passing a `t` vector to public functions. The predict step adapts its coefficients via Lagrange interpolation when the lifting scheme is interpolating (`degree >= 0`). Non-interpolating wavelets (`db2`, `cdf97`) ignore `t` and emit a warning.

→ *See: `04-irregular-grid-design.md` (full detail)*

---

## 7. Adaptive Thresholding

The threshold is estimated from the finest-level detail coefficients via a MAD estimator:

$$\hat{\sigma} = \frac{\text{MAD}(d_1)}{0.6745}$$

$$\lambda_1 = \beta \cdot \hat{\sigma} \cdot \sqrt{2 \log n}$$

$$\lambda_k = \lambda_{k-1} \cdot \frac{k - 1}{k + \alpha - 1}$$

Three shrinkage methods are available: `hard`, `soft`, and `semisoft`. In causal mode, the threshold is recomputed every `update_freq` samples from the finest-level details of the current window.

→ *See: `05-adaptive-thresholding.md` (MAD estimator, parameters α and β, shrinkage methods, offline vs causal behaviour)*

---

*Detail documents to be expanded:*
- `implementation/wavelet-engine.md` — ring buffer, zero-allocation design, XPtr
