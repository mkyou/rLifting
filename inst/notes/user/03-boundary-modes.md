# rLifting — Boundary Handling

At every decomposition level, the predict and update steps apply a filter whose window may extend beyond the signal boundaries. How those out-of-bounds positions are filled defines the boundary mode. This document describes when the problem arises, the five available modes, and how they interact with wavelets, irregular grids, and each operation mode.

---

## 1. When Boundary Handling Matters

During a predict step at position $i$, the filter reads neighbours at indices $i + \text{start\_idx}, \ldots, i + \text{start\_idx} + k - 1$. For positions near the edges, some of those indices fall outside $[0, n)$. The boundary mode determines what value is returned for those out-of-bounds reads.

Two factors amplify the effect:

- **Filter length**: a wavelet with $k$ predict coefficients contaminates up to $k - 1$ samples near each boundary per level.
- **Decomposition depth**: boundary effects propagate inward with each level. At level $j$, the contaminated zone spans roughly $(k - 1) \cdot 2^j$ samples from each edge of the original signal.

In causal/stream mode, both edges of the sliding window are virtual: neither the samples before the oldest nor those after the newest sample exist. Boundary handling affects every processed window.

---

## 2. The Five Modes

### 2.1 Symmetric (default)

Out-of-bounds positions are filled by mirror reflection:

$$x[-1] = x[0],\quad x[-2] = x[1],\quad x[n] = x[n-1],\quad x[n+1] = x[n-2],\ldots$$

**When to use.** Smooth signals without strong trends at the boundary. The default for most use cases.

**Advantage.** No discontinuity is introduced at the boundary; the virtual extension is always consistent with the local signal structure.

**Disadvantage.** If the signal has a strong slope or curvature at the edge, the reflection creates an artificial inflection point.

---

### 2.2 Periodic

Out-of-bounds positions wrap around: $x[-1] = x[n-1]$, $x[n] = x[0]$.

**When to use.** Signals that are genuinely periodic (cardiac cycles, rotating machinery, seasonal data) where the last sample connects naturally to the first.

**Advantage.** Exact for periodic signals; no boundary artifact if the period matches the signal length.

**Disadvantage.** Creates a hard discontinuity if the signal does not match at the endpoints. In causal/stream mode, "wrapping" means using the oldest sample in the window to fill the right boundary — rarely meaningful physically.

---

### 2.3 Zero

Out-of-bounds positions are set to zero.

**When to use.** Signals that are known to be zero outside the observation window (e.g. finite-duration pulses).

**Advantage.** Simplest implementation; no extrapolation assumptions.

**Disadvantage.** Introduces a step discontinuity at both edges, which can produce large detail coefficients (false edge detection) and reconstruction artefacts near the boundaries.

---

### 2.4 Local Linear

Out-of-bounds positions are extrapolated using a slope estimated by ordinary least squares (OLS) over the `ll_k` boundary samples (default 4). For the left boundary, the fit uses positions $0, 1, \ldots, \text{ll\_k}-1$; for the right boundary, positions $n-\text{ll\_k}, \ldots, n-1$. The extrapolated value is:

$$\hat{x}[i] = \hat{a} + \hat{b} \cdot i$$

where $\hat{a}$ and $\hat{b}$ are the OLS intercept and slope. If `ll_k > n`, it is clamped to `n` with a warning.

**When to use.** Signals with a measurable trend at the boundary (ramps, drifting baselines).

**Advantage.** Preserves the local slope; less discontinuous than zero-padding. With `ll_k > 2`, the slope estimate is robust to noise in individual boundary samples.

**Disadvantage.** Can over-extrapolate when the signal has curvature at the boundary. With `ll_k = 2` (minimum), a single noisy boundary sample fully determines the slope — among modes 1–4, this is the most sensitive to boundary noise at the minimum setting.

---

### 2.5 One-sided (`one_sided`)

Rather than extending the signal, this mode changes the filter itself. Out-of-bounds taps are dropped and the remaining coefficients are renormalized by their sum:

$$\hat{y}_i = \frac{\sum_{j \in \text{valid}} c_j \cdot x_{i+j}}{\sum_{j \in \text{valid}} c_j}$$

This is implemented in `onesided_conv` and requires explicit handling in `apply_filter_cpp`, `offline.cpp`, and `WaveletEngine.h` — it is the only mode that modifies the filter rather than the virtual signal values.

**When to use.** When no assumption about out-of-bounds values is acceptable. Particularly natural for causal/stream mode, where the right boundary of the window has no real future samples.

**Advantage.** Never invents data; only uses samples that exist. Robust to boundary noise — a noisy edge sample affects only one tap, not an extrapolated chain.

**Disadvantage.** The effective filter changes near boundaries, so the transform is no longer shift-invariant at the edges. Reconstruction (`ilwt`) must apply the same renormalized filter, which it does automatically.

**Incompatibility with irregular grids.** `one_sided` takes priority over the irregular path in C++: when `ext_mode == 5`, `onesided_conv` is called instead of `interp_predict`. Time positions in `t` are stored but **not used** for Lagrange interpolation. A warning is raised at the R level when `one_sided` is combined with `t` or `irregular = TRUE`.

---

## 3. Interaction with Wavelets

| Wavelet | Predict coefficients | Boundary contacts per level |
|:--------|:---------------------|:---------------------------|
| `haar` | 1 | 0 (no out-of-bounds read) |
| `lazy` | 0 | 0 |
| `cdf53` | 2 (`[0.5, 0.5]`) | 1 per edge |
| `db2` | 2 (`[√3]`-based) | 1 per edge |
| `dd4` | 4 | 2 per edge |
| `cdf97` | 4–6 | 2–3 per edge |

Haar and lazy have no boundary problem — the predict step reads only one neighbour at an in-bounds index. All other wavelets have at least one out-of-bounds read per edge per level; the choice of boundary mode becomes relevant from the second decomposition level onward.

For non-interpolating wavelets (`db2`, `cdf97`), the irregular path is not activated regardless of boundary mode (see `04-irregular-grid-design.md`, section 5). The boundary mode behaves identically for regular and irregular grids in those cases.

---

## 4. Implications per Operation Mode

### 4.1 Step-by-step and Offline

The signal has two real boundaries. Boundary handling affects only the outermost samples — up to $(k-1) \cdot 2^{\text{levels}}$ samples from each edge. For long signals and shallow decomposition, the effect is negligible in the interior.

`symmetric` and `local_linear` are the most common choices. `one_sided` is a safe default when no assumption about the signal beyond the observation window is warranted.

### 4.2 Causal Batch and Stream

The sliding window of size $W$ has **two virtual boundaries** at every step: the left edge (before the oldest sample) and the right edge (after the newest sample).

- The **left boundary** is usually less critical: the ring buffer fills with real data as the signal progresses, and the oldest sample is a real past observation.
- The **right boundary** is always virtual: no future samples exist. `symmetric` reflects the most recent samples back; `periodic` wraps to the oldest sample in the window (rarely physical); `zero` pads with zeros (introduces a step); `local_linear` extrapolates the slope; `one_sided` uses only the available samples.

For causal/stream applications, `one_sided` is the most defensible choice because it makes no assumption about future values. `symmetric` is a reasonable second choice if the signal is smooth near the right edge of the window. `local_linear` is also defensable b

`periodic` should generally be avoided in causal mode unless the signal period is known to equal $W$.

---

## 5. Quick Reference

| Mode | Virtual value | Best for | Avoid when |
|:-----|:-------------|:---------|:-----------|
| `symmetric` | Mirror reflection | Smooth signals | Strong slope at boundary |
| `periodic` | Wrap-around | Genuinely periodic signals | Non-periodic signals in causal mode |
| `zero` | 0 | Finite-duration pulses | Signals with non-zero boundary values |
| `local_linear` | OLS extrapolation (`ll_k` points, default 4) | Trending signals | Signals with boundary curvature |
| `one_sided` | Filter renormalisation | Causal mode; no assumption about outside | Irregular-grid Lagrange (incompatible) |
