# rLifting — Design Overview

This document is the architectural index of the `inst/notes/` tree. It describes the design decisions of the `rLifting` package — broad architectural choices, central abstractions, and pointers — and links out to the subsystem and implementation references in the same directory. For user-facing tutorials and worked examples, see the vignette series (`vignette("v01-introduction")`, `vignette("v02-thresholding-and-tuning")`, `vignette("v03-causal-stream")`, ...).

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
    ├── lwt / ilwt                                    ← multi-level decomposition / reconstruction
    ├── compute_adaptive_threshold / threshold        ← step-by-step thresholding
    ├── tune_alpha_beta                               ← SURE-based recursive-threshold tuner
    │
    ├── denoise_signal_offline                        ← batch denoising (non-causal)
    ├── denoise_signal_causal                         ← causal denoising (batch over history)
    ├── new_wavelet_stream                            ← sample-by-sample processor (closure)
    │
    └── diagnose_wavelet / validate_*                 ← scheme validation and diagnostics
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

→ *See: `vignette("v05-irregular-grids")` for wavelet classification and degree inference on non-uniform sampling.*

---

## 4. Modes of Operation

### 4.1 Step-by-step (manual pipeline)

The user may assemble the pipeline manually: `lwt` → `compute_adaptive_threshold` → `threshold` → `ilwt`. Each function dispatches to its own C++ routine. This path is useful for inspecting intermediate coefficients, applying custom thresholding logic, or experimenting with the transform independently of denoising. It is batch-only — causal processing requires a sliding window that must be managed manually, which is what sections 4.3–4.4 handle internally.

### 4.2 Offline (non-causal)

`denoise_signal_offline` has access to the full signal. Both regular and irregular modes route to `denoise_offline_cpp` — a single C++ pass (LWT + threshold + ILWT) with no R overhead. The irregular path stores per-level `t` position vectors and applies Lagrange interpolation in predict steps with `degree >= 0`. Per-sample times are on the order of 0.1 μs at n = 1000 for typical configurations — see `data(benchmark_rlifting)` (`Mode == "offline"`). The irregular path adds a per-level Lagrange-interpolation cost: median ~0.26 μs vs ~0.09 μs for the regular path (~3× overhead), measured in `data(benchmark_rlifting_irregular)` (`Mode == "offline"`).

Unlike the step-by-step pipeline of §4.1, `denoise_offline_cpp` carries its own inline lifting loops and does **not** route through `apply_filter_cpp`. This is why the boundary-mode "four code paths" rule (in `utils.h`, `utils.cpp`, `offline.cpp`, and `WaveletEngine.h`) lists `offline.cpp` separately — see `04-boundary-and-threshold.md` Part A.

### 4.3 Causal Batch

`denoise_signal_causal(signal, t = NULL, ...)` simulates sequential sample arrival over a historical signal. The irregular grid is supplied as a **vector `t`** parallel to `signal`. Calls `run_causal_batch_cpp` (`src/fast_stream.cpp`), which instantiates a `WaveletEngine` and iterates in C++, feeding `t[i]` per sample into `WaveletEngine::push_and_process`. Both regular and irregular modes run **entirely in C++**.

### 4.4 Stream (sample-by-sample)

`new_wavelet_stream` returns an R closure `processor(new_sample, t_val = NULL)` encapsulating a `WaveletEngine` via `Rcpp::XPtr`. Each call pushes one sample plus a **scalar `t_val`** (the sample's time position; defaults to the integer step index for regular grids) and returns the filtered value. No R-level allocations occur in the hot path. Both regular and irregular modes run **entirely in C++**.

The two surfaces converge at the C++ method `WaveletEngine::push_and_process(new_val, t_val, ...)`, which always takes a scalar `t_val` per sample. The R-side argument name reflects whether the caller is supplying a whole vector (`t`, batch) or a per-call scalar (`t_val`, stream).

→ *See: `vignette("v01-introduction")` for the user-facing tour and `vignette("v03-causal-stream")` for the deep dive on the two sliding-window modes. The C++ data flow per mode lives in `03-zero-allocation-engine.md`.*

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

→ *See: `vignette("v04-boundary-modes")` for the full per-mode discussion and causal-mode implications; `04-boundary-and-threshold.md` Part A documents the four mandatory C++ code paths.*

---

## 6. Irregular Grids

Signals with physically non-equispaced time positions are supported by passing a `t` vector to public functions. The predict step adapts its coefficients via Lagrange interpolation when the lifting scheme is interpolating (`degree >= 0`). Non-interpolating wavelets (`db2`, `cdf97`) ignore `t` and emit a warning.

→ *See: `vignette("v05-irregular-grids")` for the full user-facing discussion; `01-lifting-scheme-and-transform.md` §8 documents the Lagrange interpolation in the predict step (`interp_predict`) and `get_t_extrap` for boundary position extrapolation.*

---

## 7. Adaptive Thresholding

The threshold is estimated from the finest-level detail coefficients via a MAD estimator:

$$\hat{\sigma} = \frac{\text{MAD}(d_1)}{0.6745}$$

$$\lambda_1 = \beta \cdot \hat{\sigma} \cdot \sqrt{2 \log n}$$

$$\lambda_k = \lambda_{k-1} \cdot \frac{k - 1}{k + \alpha - 1}$$

Four shrinkage methods are available: `hard`, `soft`, `semisoft`, and `scad` (Antoniadis & Fan, 2001). Two threshold-selection rules drive how $\lambda$ is chosen: `universal` (VisuShrink with the recursive $\alpha/\beta$ decay shown above) and `sure` (per-level SureShrink). `tune_alpha_beta()` minimises SURE to pick $\alpha$ and $\beta$ automatically when the universal rule is used. In causal mode, the threshold is recomputed every `update_freq` samples from the finest-level details of the current window.

→ *See: `02-adaptive-thresholding.md` for the MAD-with-$0.6745$ derivation, the α/β recursion and SURE-risk derivation, the four shrinkage formulas, the `tune_alpha_beta()` joint optimiser, and the offline-vs-causal threshold-update protocol. `vignette("v02-thresholding-and-tuning")` covers the empirical reality check.*

---

## 8. Where to look next

The rest of the `inst/notes/` tree is ordered from least technical to most technical:

| Document | Scope |
|:---------|:------|
| `01-lifting-scheme-and-transform.md` | The `lifting_scheme` S3 object (fields, slot semantics, marshalling into C++), built-in wavelet table with `start_idx` / `degree`, R↔C++ boundary (`Rcpp::compileAttributes` and `RcppExports.cpp`), polyphase decomposition, predict/update math with worked examples (cdf53, db2), the three apply-filter sites in the codebase, irregular path with the full Lagrange formula, `interp_predict` cases, `get_t_extrap`, and the `lwt` S3 output object. |
| `02-adaptive-thresholding.md` | MAD with the $\hat\sigma = \mathrm{med}(\lvert d_1\rvert)/0.6745$ derivation and Fisher consistency, the universal rule with α/β recursion, the SURE rule (full risk derivation and per-level $\hat\sigma_k$), the four shrinkage formulas (hard / soft / semisoft / SCAD piecewise), the `tune_alpha_beta()` joint optimiser (objective, two-phase optimisation, box clip), causal `update_freq` protocol, and edge cases. |
| `03-zero-allocation-engine.md` | `WaveletEngine` class layout and constructor allocations, ring buffer indexing and rationale, the per-sample hot path step-by-step (output index semantics, boundary modes in the hot path, irregular path), the XPtr finalizer pattern, allocation audit per public path, and how the offline path (`denoise_offline_cpp`) deliberately differs from the engine. |
| `04-boundary-and-threshold.md` | **Part A:** formal definitions of the five boundary modes (half-sample symmetric vs whole-sample; periodic; zero; local-linear OLS with the $K^2(K^2-1)/12$ closed form; one-sided renormalised filter), the four mandatory C++ code paths, R-side dispatch sites, irregular position extrapolation via `get_t_extrap`. **Part B:** MAD via `nth_element` selection routed through the shared `compute_mad` helper (canonical averaged median across all four call sites), the three parallel shrinkage implementations across the codebase. |

The user-facing tour lives in the vignette series (`vignette("v01-introduction")` onward).
