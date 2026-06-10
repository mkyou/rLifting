# rLifting — Boundary Extension and Adaptive Thresholding (Implementation Reference)

Implementation depth for the five boundary modes and the threshold internals. This note assumes the reader has gone through `vignette("04-boundary-modes")` for the user-facing semantics and empirical comparison; here we pin down the C++ surfaces, the four mandatory code paths, the OLS extrapolation, the renormalised one-sided filter, the `nth_element` median selection, and the piecewise shrinkage formulas.

Related notes: `00-design-overview.md` (architecture index), `01-lifting-scheme-and-transform.md` (lifting/polyphase math), `02-adaptive-thresholding.md` (α/β recursion, SURE rule, tuner), `03-zero-allocation-engine.md` (engine layout, ring buffer).

---

## Part A — Boundary Extension

### A.1 The five modes

Integer enum used end-to-end (R → C++):

| Code | Name | Definition of $x[i]$ for $i \notin [0, n)$ | Mechanism |
|:----:|:-----|:--------------------------------------------|:----------|
| 1 | `symmetric` | $x_{-1} = x_0,\ x_{-2} = x_1,\ \dots,\ x_n = x_{n-1},\ x_{n+1} = x_{n-2},\dots$ | Point-wise |
| 2 | `periodic`  | $x_{i} = x_{((i \bmod n) + n) \bmod n}$ | Point-wise |
| 3 | `zero`      | $x_i = 0$ | Point-wise |
| 4 | `local_linear` | $\hat a + \hat b\, i$ from OLS line through $K=\min(\texttt{ll\_k}, n)$ nearest boundary samples | Point-wise |
| 5 | `one_sided` | undefined (signal is **not** extended); the filter is renormalised over in-bounds taps only | Filter-level |

**Symmetric is half-sample, not whole-sample.** The reflection in `get_val_safe` in `utils.h` is

```cpp
while (i < 0 || i >= n) {
    if (i < 0) i = -1 - i;
    else       i = 2 * n - 1 - i;
}
```

For $i = -1$: $i \leftarrow -1 - (-1) = 0$, so $x_{-1} = x_0$. For $i = n$: $i \leftarrow 2n - 1 - n = n - 1$, so $x_n = x_{n-1}$. Boundary samples are duplicated (half-sample symmetric, type `WPER`/`HSF` in the wavelet literature). Whole-sample symmetric would map $x_{-1} = x_1$ and $x_n = x_{n-2}$; that variant is not implemented.

**Periodic** uses C++ modulo with sign correction (`get_val_safe` in `utils.h`):

```cpp
int idx = i % n; if (idx < 0) idx += n; return x[idx];
```

This is exact periodisation, not the same as periodic-then-symmetric padding used in some DWT libraries.

### A.2 The four mandatory code paths

A single boundary-mode change must propagate to **four C++ sites**. They cannot share a common loop because each runs in a different control-flow context (single-shot batch vs. per-level inline loops vs. per-sample workspace updates).

| File:function | Purpose | Modes 1–4 | Mode 5 | Used by |
|:--------------|:--------|:----------|:-------|:--------|
| `inst/include/rLifting/utils.h::get_val_safe` | Single-index virtual value | Implemented here | N/A (cannot express filter-level change) | All four sites for point-wise modes |
| `apply_filter_cpp` in `src/utils.cpp` | Generic per-position convolution, called from `lwt_cpp`/`ilwt_cpp` | Loop calls `get_val_safe` per tap | Explicit `if (ext_mode == 5)` branch calls `onesided_conv` (in `apply_filter_cpp` in `utils.cpp`) | `lwt`, `ilwt` |
| `denoise_offline_cpp` in `src/offline.cpp` | Inlined LWT+threshold+ILWT in one pass | Inline loop calls `get_val_safe` per tap (predict and update bodies in `denoise_offline_cpp` in `offline.cpp`) | `use_os = (ext_mode == 5)` switches each predict/update body to `onesided_conv` (forward and inverse passes in `denoise_offline_cpp` in `offline.cpp`) | `denoise_signal_offline` |
| `push_and_process` in `inst/include/rLifting/WaveletEngine.h` | Per-sample forward and inverse passes against the ring-buffer workspace | Inline loops call `get_val` (a wrapper around `get_val_safe`) | Same `use_os` switch in forward and inverse passes (`push_and_process` in `WaveletEngine.h`) | `denoise_signal_causal`, `new_wavelet_stream` |

The forward and inverse halves of `push_and_process` are structurally identical but written out as two distinct loops; both must mirror the boundary logic.

**Rule of thumb.** Adding a *point-wise* mode (a new way to invent $x[i]$ at a single virtual index) touches only `get_val_safe`. Adding a *filter-level* mode (anything that alters the filter as a function of which taps are in-bounds) touches all four files.

Forgetting a site does not raise a compile error — it silently degrades that usage mode to the default `get_val_safe` branch (or to whichever branch the switch falls through to). Always grep the codebase for `ext_mode ==` and `use_os` when changing the boundary system.

### A.3 R-side dispatch

Four R wrappers convert the string mode to the integer enum before crossing into C++:

| File | Function | Switch location |
|:-----|:---------|:----------------|
| `R/lwt.R` | `lwt()` | boundary `switch` in `lwt` |
| `R/ilwt.R` | `ilwt()` | boundary `switch` in `ilwt` (reads `lwt_obj$extension`) |
| `R/denoising_offline.R` | `denoise_signal_offline()` | boundary `switch` in `denoise_signal_offline` |
| `R/realtime_denoising.R` | `denoise_signal_causal()` and `new_wavelet_stream()` | boundary `switch` in `denoise_signal_causal` and in `new_wavelet_stream` |

Each switch falls through to `1L` (symmetric) on an unknown name. The R layer also raises a warning if `extension = "one_sided"` is combined with a non-NULL `t` (irregular grid) — `one_sided` bypasses the Lagrange path in C++ regardless of `degree`, since the explicit branch is checked before the irregular branch.

### A.4 Mode 4 (local_linear) — OLS extrapolation

`get_val_safe` in `utils.h` fits a least-squares line through $K = \min(\texttt{ll\_k}, n)$ samples nearest the boundary and evaluates it at the requested out-of-bounds index.

Left boundary ($i < 0$): fit through positions $0, 1, \dots, K-1$ with values $x_0, \dots, x_{K-1}$. Right boundary ($i \geq n$): fit through positions $n-K, \dots, n-1$.

Closed form (consecutive integer abscissae):

$$\hat b \;=\; \frac{K \sum_j j\, x_j \;-\; (\sum_j j)(\sum_j x_j)}{K \sum_j j^2 \;-\; (\sum_j j)^2}, \qquad \hat a \;=\; \frac{\sum_j x_j - \hat b \sum_j j}{K}$$

Output: $\hat a + \hat b \cdot i$.

**Edge cases.**

- If $n < 2$: return $x_0$ (no line can be fit; the early-return guard in `get_val_safe` in `utils.h`).
- `ll_k` is clamped: `int k = std::min(ll_k, n); if (k < 2) k = 2;` — so $K \geq 2$ is enforced in C++. The R layer additionally warns when the user-supplied `ll_k` exceeds the signal length, but does not error.
- Degenerate denominator: $K \sum j^2 - (\sum j)^2 = K^2(K^2 - 1)/12$ for consecutive integers starting at 0, which is positive for $K \geq 2$. The guard `abs(denom) < 1e-15` (in `get_val_safe` in `utils.h`) catches the constant-signal case (all $x_j$ equal — slope is mathematically zero, but the formula collapses to $0/0$ when accumulated in floating point); in that case the boundary value $x_0$ (left) or $x_{n-1}$ (right) is returned.

The accumulators `st`, `sty` etc. are recomputed on every call to `get_val_safe`. This is $O(K)$ per virtual sample and could in principle be hoisted out of the convolution loop, but it has never shown up as a hotspot — the loop only fires for out-of-bounds indices, which is a constant per level regardless of $n$.

### A.5 Mode 5 (one_sided) — filter renormalisation

`one_sided` cannot live in `get_val_safe`. `get_val_safe`'s contract is "given index $i$, return a virtual value"; `one_sided` does not invent a virtual value, it removes taps from the filter and rescales the surviving taps so their sum still equals the original filter sum. The output at position $p$ depends on the *set* of in-bounds indices $\{p + \texttt{start\_idx}, \dots, p + \texttt{start\_idx} + k - 1\}$, which is a property of the filter window, not a function of a single index.

The implementation in `onesided_conv` in `utils.h`:

```cpp
int first = pos + start_idx, last = first + k - 1;
if (first >= 0 && last < n) {          // fast path: window fully in bounds
    double s = 0.0;
    for (int j = 0; j < k; j++) s += x[first + j] * c[j];
    return s;
}
double vs = 0.0, ws = 0.0;             // boundary path
for (int j = 0; j < k; j++) {
    int idx = first + j;
    if (idx >= 0 && idx < n) { vs += x[idx] * c[j]; ws += c[j]; }
}
return (ws > 1e-15) ? vs / ws : 0.0;
```

- **Fast path** (interior positions): identical to standard convolution, no overhead. Triggered when the entire filter window is in bounds.
- **Boundary path**: accumulate value-sum $vs$ and weight-sum $ws$ over the in-bounds taps only, return $vs / ws$. The renormalisation preserves the original filter gain — e.g. for the CDF 5/3 predict $c = (0.5, 0.5)$ with one tap missing, the surviving coefficient becomes $1.0$, doubling the contribution of the boundary sample. This amplification is the source of the catastrophic MSE seen for `one_sided` + long filters in `vignette("04-boundary-modes")` §5.2.
- **Empty-window guard**: `ws > 1e-15` handles the case where every tap is out of bounds (signal shorter than the filter at extreme positions); returns 0. The threshold is a numerical safety margin rather than a mathematical condition — for any realistic filter at least one tap is in bounds.

The four sites (Table A.2) each branch on `ext_mode == 5` *before* the `get_val_safe` loop and switch to `onesided_conv`. In `apply_filter_cpp` the branch is at function entry (in `apply_filter_cpp` in `utils.cpp`); in `offline.cpp` and `WaveletEngine.h` a local `bool use_os = (ext_mode == 5)` is set once per call and checked inside each predict/update step.

### A.6 Position extrapolation for irregular grids

For irregular-grid predict steps with `degree >= 0`, the Lagrange interpolator needs both the values and the physical positions at all neighbour indices, including out-of-bounds ones. Positions are extrapolated by `get_t_extrap` in `utils.h`:

```cpp
if (idx < 0)  return t[0]   + (double)idx           * (t[1] - t[0]);
if (idx >= n) return t[n-1] + (double)(idx - (n-1)) * (t[n-1] - t[n-2]);
```

Position extrapolation is **unconditionally linear**, regardless of `ext_mode`. The other modes are not applicable here:

- Symmetric reflection would produce non-monotonic positions ($t_{-1} = t_0$, $t_{-2} = t_1$, ...) — repeated and reversed times.
- Periodic would wrap times around to earlier values, also non-monotonic.
- Zero would put neighbours at $t = 0$, collapsing the Lagrange denominators.

Linear extrapolation from the nearest boundary spacing is the unique extension that preserves strict monotonicity of $t$, which is what `interp_predict` in `utils.h` requires for finite Lagrange weights.

The split of responsibilities is clean: `get_val_safe` answers *what value lives at the virtual position*; `get_t_extrap` answers *where in time the virtual position lies*. The boundary mode affects only the former.

---

## Part B — Threshold and Shrinkage Internals

This part is the implementation reference. The higher-level discussion of the α/β recursion, the universal-vs-SURE rule choice, and the `tune_alpha_beta()` pipeline lives in `02-adaptive-thresholding.md`; see that note first if you need *why*. Here we document *how*.

### B.1 MAD via `std::nth_element`

The recursive-universal rule starts from $\hat\sigma = \mathrm{MAD}(d_1) / 0.6745$. Computing the median by full sort is $O(n \log n)$; only the median is needed, so `std::nth_element` is used everywhere it appears:

All four call sites route through a single inline helper, `compute_mad` in [`inst/include/rLifting/utils.h`](../../inst/include/rLifting/utils.h), which implements the canonical averaged median (matching R's `mad()`): one `nth_element` to find the upper middle, a second one (over the truncated `[begin, mid)` range) for the lower middle when $n$ is even, then their mean.

| Site | File | Surface |
|:-----|:-----|:--------|
| `compute_mad_cpp` / `compute_thresholds_cpp` | `src/adaptative.cpp` | R-exposed; `compute_mad_cpp(x)` returns the median directly, `compute_thresholds_cpp(d1, max_level, alpha, beta)` returns the recursive $\lambda$ vector. |
| `compute_thresholds_internal` | `src/offline.cpp` | Offline universal rule, single C++ pass. |
| `WaveletEngine::update_thresholds` | `inst/include/rLifting/WaveletEngine.h` | Causal/stream hot path; called every `update_freq` samples after warm-up. |
| `compute_sure_lambda_level` | `inst/include/rLifting/utils.h` | SURE rule, both offline and engine; receives a fresh `abs_d` copy because the calling code full-sorts it afterwards. |

**Why `nth_element`.** Expected $O(n)$ in-place partial partitioning: after the call, position `mid` holds the value that would appear there in a sorted array, with smaller values to its left and larger to its right (unordered). The even-$n$ second call partitions only the first half, so total work stays $O(n)$. No allocation beyond the caller-owned `abs_*` buffer. For $n_1 = W / 2^L$ in causal mode (typical 8–128), the absolute speedup over a full sort is small; the structural win is that the hot path on the engine side does not allocate a sorted copy.

### B.2 The $\sigma$ guard

All three universal-rule MAD sites (`compute_thresholds_cpp` in `adaptative.cpp`, `compute_thresholds_internal` in `offline.cpp`, and `update_thresholds` in `WaveletEngine.h`) test:

```cpp
if (sigma < 1e-15) { /* return zero lambdas */ }
```

The `1e-15` threshold (rather than strict zero) is a floating-point safety margin: theoretically-zero detail coefficients can leave a residual on the order of machine epsilon after multiple lifting steps. Without the guard, a perfectly smooth signal would yield a near-zero $\sigma$, a near-zero $\lambda$, and would still pass the threshold test — wasted computation with no shrinkage. With $\lambda = 0$ explicitly, the condition `abs_val < lam` is false for every non-zero coefficient and they all survive unchanged. Same result, well-defined path.

### B.3 Shrinkage formulas

Four methods. The threshold check `abs_val < lam → 0` is identical across all four; the surviving branch differs. Let $d$ be a coefficient, $\lambda > 0$ the threshold, and $s = \mathrm{sign}(d)$.

**Hard** (`threshold_hard_cpp` in `thresholding.cpp`, also inlined in the engines as the implicit "do nothing" branch):

$$T_H(d, \lambda) \;=\; \begin{cases} 0 & |d| < \lambda \\ d & |d| \geq \lambda \end{cases}$$

**Soft** (`threshold_soft_cpp` in `thresholding.cpp`, engines `denoise_offline_cpp` in `offline.cpp` and `push_and_process` in `WaveletEngine.h`):

$$T_S(d, \lambda) \;=\; \begin{cases} 0 & |d| < \lambda \\ s\,(|d| - \lambda) & |d| \geq \lambda \end{cases}$$

**Semisoft / hyperbolic** (`threshold_semisoft_cpp` in `thresholding.cpp`, engines `denoise_offline_cpp` in `offline.cpp` and `push_and_process` in `WaveletEngine.h`):

$$T_{SS}(d, \lambda) \;=\; \begin{cases} 0 & |d| < \lambda \\ s\,\sqrt{d^2 - \lambda^2} & |d| \geq \lambda \end{cases}$$

The square root is always real on this branch because $|d| \geq \lambda \Rightarrow d^2 \geq \lambda^2$. `lam_sq = lam * lam` is hoisted out of the inner loop (computed once per level) in the engines.

**SCAD** (Antoniadis & Fan 2001; Fan & Li 2001) (`threshold_scad_cpp` in `thresholding.cpp`, engines `denoise_offline_cpp` in `offline.cpp` and `push_and_process` in `WaveletEngine.h`). Three regions with canonical shape parameter $a = 3.7$. The SCAD shape parameter `a` is honoured by all four shrinkage entry points: the standalone `threshold_scad()` wrapper, `denoise_signal_offline()` (via the `scad_a` parameter of `denoise_offline_cpp`), `denoise_signal_causal()`, and `new_wavelet_stream()` (the latter two via the `scad_a` member of `WaveletEngine`, set at engine construction by `create_engine_cpp` / `run_causal_batch_cpp`). All sites validate $a > 2$:

$$T_{SC}(d, \lambda) \;=\; \begin{cases}
0 & |d| \leq \lambda \\
s\,(|d| - \lambda) & \lambda < |d| \leq 2\lambda \\
\dfrac{(a - 1)\, d \;-\; s\, a\lambda}{a - 2} & 2\lambda < |d| \leq a\lambda \\
d & |d| > a\lambda
\end{cases}$$

Continuous at $|d| = \lambda$, $2\lambda$, and $a\lambda$. The fourth region (identity) removes the soft-threshold bias on large coefficients while preserving the sparsity-inducing zero region. In the engine inlining (the SCAD branch in `denoise_offline_cpp` in `offline.cpp` and in `push_and_process` in `WaveletEngine.h`) the identity region has no else branch — `det[i]` is left at its current value, which is `val`.

### B.4 Where the formulas live

Three implementations of each rule are kept in lockstep:

| Form | File | Use |
|:-----|:-----|:----|
| Standalone Rcpp export (`threshold_*_cpp`) | `src/thresholding.cpp` | Direct calls from R (`threshold_hard`/`soft`/`semisoft` wrappers in `R/thresholding.R`); step-by-step pipeline |
| Inlined in offline engine | the shrinkage block of `denoise_offline_cpp` in `src/offline.cpp` | Single-pass denoise |
| Inlined in causal engine | the shrinkage block of `push_and_process` in `inst/include/rLifting/WaveletEngine.h` | Per-sample stream |

The inlined copies avoid a function call per coefficient and avoid materialising a full `NumericVector` return. The dispatch on `method` is a string comparison once per level (outside the coefficient loop) — negligible cost relative to the convolutions. The SCAD shape parameter `a` is per-call configurable at all three sites: `threshold_scad_cpp` in `thresholding.cpp`, `denoise_offline_cpp` in `offline.cpp` (via the `scad_a` parameter), and `push_and_process` in `WaveletEngine.h` (via the `scad_a` member of the engine).

If you change a shrinkage formula, all three copies must change together. The standalone exports are also covered by direct unit tests in `tests/testthat/test-test-thresholding.R`; the engine copies are exercised indirectly through `denoise_signal_offline`, `denoise_signal_causal`, and `new_wavelet_stream` test suites.

### B.5 The threshold cache (`update_freq`)

In `WaveletEngine`, the per-level $\lambda$ values are stored in `current_lambdas` (a pre-allocated vector of `levels` doubles) and reused between updates (the threshold-update gate in `push_and_process` in `WaveletEngine.h`):

```cpp
bool should_update = (update_freq > 0)
                        ? (step_iter % update_freq == 0)
                        : !lambdas_initialized;
if (should_update) {
  update_thresholds(alpha, beta, threshold_method);
  lambdas_initialized = true;
}
```

`update_freq = 0` (warned in `R/realtime_denoising.R`; negative values are rejected with an error) freezes thresholds after warm-up; `lambdas_initialized` ensures the one-time initial update fires before the freeze.

Between updates the cache holds the thresholds from the most recent recomputation. For stationary noise the cache is always valid; for non-stationary noise the lag is at most `update_freq` samples. `step_iter` is incremented in the R closure (the stream closure in `R/realtime_denoising.R`), not inside the engine, so changing `update_freq` between calls is meaningful only for the stream API, where the parameter is captured at creation time and frozen.

`update_thresholds` in `WaveletEngine.h` branches once on `threshold_method`: SURE recomputes all $L$ per-level $\lambda_j = $ `compute_sure_lambda_level(work_detail[j])`; universal recomputes $\lambda_1$ from $\mathrm{MAD}(d_1)$ and applies the α/β recursion (see `02-adaptive-thresholding.md` §3) to fill the rest. Both paths reuse the same pre-allocated workspaces, no allocations.

---

## Cross-references

- **Lifting and polyphase math, including how the predict step calls into the boundary system**: `01-lifting-scheme-and-transform.md`.
- **α/β recursion, universal vs. SURE rule, `tune_alpha_beta()` design**: `02-adaptive-thresholding.md`. This is where to look for *why* the threshold formula is what it is; this note covers the *implementation*.
- **Ring buffer, pre-allocation, XPtr lifetime, per-sample data flow**: `03-zero-allocation-engine.md`.
- **User-facing tour of the five modes, empirical MSE comparisons across signals and wavelets, decision guide**: `vignette("04-boundary-modes")`.
