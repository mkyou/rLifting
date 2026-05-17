# rLifting — Boundary Extension and Adaptive Thresholding

This document describes the implementation of boundary extension and adaptive thresholding: the architectural reason modes 1–4 and mode 5 are handled differently, the design constraints that forced four separate code paths, and the algorithmic choices behind the MAD estimator and its causal update mechanism.

---

## Part A: Boundary Extension

### 1. The Single-Index Abstraction and Its Limits

The primary boundary handling function is `get_val_safe` in `inst/include/rLifting/utils.h`. Its interface is:

```cpp
inline double get_val_safe(
    const std::vector<double>& x, int i, int n, int mode, int ll_k = 2
)
```

It answers a single question: given signal `x` of length `n` and boundary mode `mode`, what is the virtual value at index `i`? For modes 1–4, this is well-defined as a function of `i` alone:

- **Symmetric** (1): reflect `i` through the nearest boundary until it falls in `[0, n)`.
- **Periodic** (2): `i mod n` (with sign correction for negative indices).
- **Zero** (3): return `0.0` unconditionally.
- **Local linear** (4): extrapolate using an OLS line fit through the `ll_k` nearest boundary samples.

All four are *point-wise* extensions: the virtual value at index `i` depends only on `i` and the signal, not on which other indices the filter is currently accessing. This allows `get_val_safe` to be called independently per tap in the convolution loop:

```cpp
for (int j = 0; j < k; j++) {
    sum += get_val_safe(x, i + start_idx + j, n, ext_mode, ll_k) * c[j];
}
```

### 2. Why One-sided Cannot Use get_val_safe

Mode 5 (`one_sided`) breaks the point-wise abstraction. Its behaviour is not "return a virtual value at index `i`" — it is "drop all out-of-bounds taps from the filter and renormalize by the sum of in-bounds weights". The output at position `i` depends on which specific subset of taps `{i + start_idx, ..., i + start_idx + k - 1}` happen to be in-bounds, which is a property of the filter window as a whole.

`onesided_conv` (also in `utils.h`) implements this:

```cpp
inline double onesided_conv(
    const std::vector<double>& x, int n,
    const double* c, int k, int start_idx, int pos
) {
    int first = pos + start_idx;
    int last  = first + k - 1;
    if (first >= 0 && last < n) {   // fast path: all taps in bounds
        double s = 0.0;
        for (int j = 0; j < k; j++) s += x[first + j] * c[j];
        return s;
    }
    double vs = 0.0, ws = 0.0;
    for (int j = 0; j < k; j++) {
        int idx = first + j;
        if (idx >= 0 && idx < n) { vs += x[idx] * c[j]; ws += c[j]; }
    }
    return (ws > 1e-15) ? vs / ws : 0.0;
}
```

The fast path — when the entire filter window falls within `[0, n)` — avoids the renormalization overhead; it is taken for all interior positions and is equivalent to standard convolution. The boundary path accumulates partial sums `vs` (weighted values) and `ws` (weight sum), then returns `vs/ws`. The guard `ws > 1e-15` handles the degenerate case where all taps fall outside the signal — it returns `0.0` rather than a division-by-zero.

The architectural consequence is that `one_sided` cannot be handled in `get_val_safe`. It requires the entire filter window, not a single virtual sample. This forces it to be treated as a separate execution branch at every site that performs convolution.

### 3. The Four Mandatory Code Paths

There are four independent C++ sites where the lifting loop runs:

| Site | File | Used by |
|:-----|:-----|:--------|
| `apply_filter_cpp` | `src/utils.cpp` | `lwt_cpp`, `ilwt_cpp` (step-by-step path) |
| `denoise_offline_cpp` | `src/offline.cpp` | `denoise_signal_offline` |
| `WaveletEngine::push_and_process` (forward) | `WaveletEngine.h` | causal and stream modes |
| `WaveletEngine::push_and_process` (inverse) | `WaveletEngine.h` | causal and stream modes |

Each site must independently branch on `ext_mode == 5` to call `onesided_conv` instead of the `get_val_safe` loop. Adding a new boundary mode that requires the full filter (like `one_sided`) must be done in all four places; adding a new point-wise mode (like a new type of extrapolation) requires only updating `get_val_safe`.

The forward and inverse passes of `WaveletEngine` are counted separately because they are distinct loops in `push_and_process` — the inverse loop reverses step order and sign, but the boundary logic is structurally identical.

### 4. The Local Linear OLS Implementation

For mode 4, `get_val_safe` fits a least-squares line through the `ll_k` nearest boundary samples and extrapolates. For the left boundary (index `i < 0`), it fits through positions `0, 1, ..., ll_k - 1`; for the right boundary (`i >= n`), through positions `n - ll_k, ..., n - 1`. The normal equations for consecutive integer positions simplify to:

$$\hat{b} = \frac{k \sum_{j} j \cdot x[j] - \sum_j j \cdot \sum_j x[j]}{k \sum_j j^2 - (\sum_j j)^2}, \quad \hat{a} = \frac{\sum_j x[j] - \hat{b} \sum_j j}{k}$$

The denominator $k \sum j^2 - (\sum j)^2$ equals $k^2(k^2-1)/12$ for consecutive integers starting at 0, and is always positive for $k \geq 2$. The guard `abs(denom) < 1e-15` catches the degenerate case of a constant signal (all $x[j]$ equal), where slope is undefined — in that case, the endpoint value is returned.

The parameter `ll_k` is clamped to `min(ll_k, n)` inside `get_val_safe`, so it is safe to pass any positive integer. The R layer emits a warning when the user-specified `ll_k` exceeds the signal length, but the C++ layer clamps silently.

### 5. Position Extrapolation: get_t_extrap

For irregular-grid predict steps, the Lagrange interpolation requires the physical time positions of all neighbours, including out-of-bounds ones. `get_t_extrap` provides these independently of `get_val_safe`:

```cpp
inline double get_t_extrap(const std::vector<double>& t, int idx, int n) {
    if (idx >= 0 && idx < n) return t[idx];
    if (idx < 0) {
        double dt = t[1] - t[0];
        return t[0] + (double)idx * dt;
    }
    double dt = t[n-1] - t[n-2];
    return t[n-1] + (double)(idx - (n-1)) * dt;
}
```

Position extrapolation is unconditionally linear, regardless of the `ext_mode` selected for signal values. The reason is that time coordinates must remain strictly monotonic: any other extension (symmetric reflection, periodic wrapping, zero-padding) would produce repeated or reversed positions, causing zero denominators in the Lagrange formula. Linear extrapolation from the nearest boundary spacing is the unique extension that preserves strict ordering.

This means the boundary mode affects only the signal values used in Lagrange interpolation, not the positions at which they are evaluated. The two functions have orthogonal responsibilities: `get_val_safe` answers "what value does the signal have at this virtual position?" and `get_t_extrap` answers "where in time is this virtual position?".

---

## Part B: Adaptive Thresholding

### 6. Two Structurally Identical Implementations

The threshold computation is implemented twice:

- `compute_thresholds_internal` in `src/offline.cpp` — called once per full-signal denoising call by `denoise_offline_cpp`.
- `WaveletEngine::update_thresholds` in `WaveletEngine.h` — called every `update_freq` samples by `push_and_process`.

The two functions are structurally identical — same MAD estimator, same universal threshold formula, same recursive decay. They are separate because `offline.cpp` is designed to be self-contained (no dependency on `WaveletEngine`), and the data each function receives differs: `compute_thresholds_internal` receives the full finest-level detail vector from the single-pass transform; `update_thresholds` reads from the pre-allocated `work_detail[0]` workspace.

The threshold parameters `alpha`, `beta`, and `max_level` are passed as arguments in both cases; neither function holds state.

### 7. O(n) MAD via nth_element

The MAD estimator requires the median of `|d_1|`. A full sort is $O(n \log n)$; `std::nth_element` finds the $k$-th order statistic in $O(n)$ expected time by partial partitioning:

```cpp
std::vector<double> abs_d1(n1);
for (int i = 0; i < n1; i++) abs_d1[i] = std::abs(d1[i]);
int mid = n1 / 2;
std::nth_element(abs_d1.begin(), abs_d1.begin() + mid, abs_d1.end());
double mad = abs_d1[mid];
```

After the call, `abs_d1[mid]` holds the value that would appear at position `mid` in a sorted array, with smaller values to its left and larger values to its right (in no particular order). For $n_1 = W / 2^L$ (the number of finest-level coefficients), typical values are 8–128 elements in causal mode; the $O(n)$ vs $O(n \log n)$ difference is small in absolute terms but avoids allocating a sorted copy.

Note that `mid = n1 / 2` computes the lower median for even-length arrays (integer division). For the MAD estimator, this is conventional and does not affect the asymptotic properties of $\hat{\sigma}$.

### 8. The sigma Guard

```cpp
double sigma = mad / 0.6745;
if (sigma < 1e-15) {
    std::fill(current_lambdas.begin(), current_lambdas.end(), 0.0);
    return;
}
```

When all finest-level detail coefficients are near zero — a smooth signal where the LWT removes nearly all energy from the finest level — the MAD is zero and the estimated noise level is zero. Setting all thresholds to zero causes the threshold condition `abs_val < lam` to be false for all non-zero coefficients, so every coefficient survives unchanged. This is the correct behaviour: if the signal is so smooth that no noise is detected, no thresholding should be applied.

The threshold `1e-15` (rather than strict zero) accounts for floating-point rounding: a theoretically-zero detail coefficient may have a residual of order machine epsilon after multiple lifting steps. Without this guard, a signal with tiny numerical residuals would yield a very small but non-zero `sigma` and then an extremely small `lambda`, which would threshold nothing — the same outcome as setting `lambda = 0`, but reached via a numerically unstable path.

### 9. The Threshold Cache and update_freq

In `WaveletEngine`, thresholds are stored in `current_lambdas` (a pre-allocated vector of `levels` doubles) and reused between updates:

```cpp
if (step_iter % update_freq == 0) update_thresholds(alpha, beta);

for (int j = 0; j < levels; j++) {
    double lam    = current_lambdas[j];
    ...
}
```

Between updates, `current_lambdas` holds the thresholds from the most recent call to `update_thresholds`. This caching is safe as long as the noise level does not change faster than `update_freq` samples. For stationary noise, the cache is always valid. For non-stationary noise, the lag introduced by `update_freq > 1` means the threshold reflects the noise level `update_freq` samples ago — an approximation that trades accuracy for reduced computation.

The check `step_iter % update_freq == 0` is performed before the shrinkage loop on every call to `push_and_process`. `step_iter` is maintained in the R closure (not inside `WaveletEngine`) and passed as an argument — this design allows `update_freq` and other parameters to be changed between samples without reconstructing the engine.

### 10. Shrinkage Implementation

All three shrinkage methods share a common threshold check:

```cpp
if (abs_val < lam) {
    det[i] = 0.0;
} else {
    if (method == "soft") {
        det[i] = (val > 0) ? (abs_val - lam) : -(abs_val - lam);
    } else if (method == "semisoft") {
        double s = std::sqrt(val * val - lam_sq);
        det[i] = (val > 0) ? s : -s;
    }
    // hard: no else branch — det[i] is already val, unchanged
}
```

Hard thresholding requires no computation for surviving coefficients — they are left in place. Soft thresholding subtracts `lam` from the absolute value. Semisoft computes `sqrt(val^2 - lam^2)`, which is guaranteed non-negative since we only reach this branch when `abs_val >= lam`. No NaN guard is needed for the sqrt.

The `lam_sq = lam * lam` precomputation is outside the coefficient loop (once per level), avoiding a multiply inside the inner loop for semisoft. For hard and soft thresholding, `lam_sq` is computed but unused — a minor inefficiency that avoids branching on method before the inner loop.
