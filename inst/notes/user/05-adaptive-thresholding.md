# rLifting — Adaptive Thresholding

This document covers how the adaptive threshold is estimated, what the parameters α and β control, the three shrinkage methods, and how thresholding behaves differently in offline and causal modes.

---

## 1. Why Adaptive Thresholding

A fixed threshold treats all signals and noise levels equally. The adaptive approach estimates noise variance directly from the signal's own finest-level detail coefficients, so the threshold scales automatically with the noise level present in the data. This makes it robust across signals with different amplitudes without requiring manual tuning.

---

## 2. Noise Estimation via MAD

The noise standard deviation is estimated from the detail coefficients at the finest level ($d_1$) using the Median Absolute Deviation (MAD):

$$\hat{\sigma} = \frac{\text{MAD}(d_1)}{0.6745}$$

The constant $0.6745$ is the theoretical MAD of a standard normal distribution ($\Phi^{-1}(0.75)$). It normalises the MAD to a consistent estimator of $\sigma$ under Gaussian noise. The finest level is used because it captures the highest-frequency content — where noise dominates over signal energy in most practical cases.

**Why MAD instead of standard deviation?** MAD is resistant to outliers (large signal coefficients do not inflate the estimate), making the threshold more stable on signals with discontinuities such as blocks or bumps.

---

## 3. The Universal Threshold

The threshold at the finest level is the VisuShrink universal threshold (Donoho & Johnstone, 1994), scaled by $\beta$:

$$\lambda_1 = \beta \cdot \hat{\sigma} \cdot \sqrt{2 \log n}$$

where $n$ is the number of finest-level detail coefficients. The $\sqrt{2 \log n}$ factor comes from the expected maximum of $n$ i.i.d. standard normal variables — it is the level above which a pure-noise coefficient is unlikely to appear. $\beta$ is a scale factor that allows tightening or loosening the threshold relative to this theoretical bound.

**Parameter β:**
- $\beta = 1$ corresponds exactly to the universal threshold.
- $\beta < 1$ lowers the threshold → more coefficients survive → less smoothing, more noise retained.
- $\beta > 1$ raises the threshold → fewer coefficients survive → more smoothing, risk of over-smoothing. Default: 1.2.

---

## 4. Recursive Per-level Thresholds

Applying the same $\lambda_1$ to all decomposition levels is suboptimal: coarser levels contain more signal energy and fewer noise-dominated coefficients. The threshold decays recursively across levels:

$$\lambda_k = \lambda_{k-1} \cdot \frac{k - 1}{k + \alpha - 1}$$

**Parameter α:**
- Controls how fast the threshold decays across levels.
- $\alpha \to \infty$: $\lambda_k \to \lambda_{k-1}$ for all $k$ — flat threshold across levels.
- $\alpha = 0$: maximum decay — $\lambda_k = \lambda_{k-1} \cdot \frac{k-1}{k}$.
- Default: 0.3. A small positive value gives a moderate, physically motivated decay.

The recursion starts at $k = 2$ (level 2 uses level 1's threshold as the base). Level 1 always uses $\lambda_1$.

---

## 5. Shrinkage Methods

All three methods zero out coefficients below $\lambda$ and differ in how they treat coefficients above it.

### 5.1 Hard thresholding

$$\hat{d} = \begin{cases} d & |d| \geq \lambda \\ 0 & |d| < \lambda \end{cases}$$

Keeps large coefficients intact; zeroes the rest. Produces sharp reconstructions but can introduce Gibbs-like ringing at discontinuities because of the abrupt transition at $\lambda$.

### 5.2 Soft thresholding

$$\hat{d} = \begin{cases} \text{sign}(d)(|d| - \lambda) & |d| \geq \lambda \\ 0 & |d| < \lambda \end{cases}$$

Shrinks all surviving coefficients toward zero by $\lambda$. The continuity at $\lambda$ eliminates Gibbs artefacts but introduces a systematic bias — large coefficients are under-estimated by $\lambda$.

### 5.3 Semisoft thresholding (default)

$$\hat{d} = \begin{cases} \text{sign}(d)\sqrt{d^2 - \lambda^2} & |d| \geq \lambda \\ 0 & |d| < \lambda \end{cases}$$

A smooth, bias-reducing compromise between hard and soft. At $|d| = \lambda$ the output is zero (like soft); as $|d| \to \infty$, $\hat{d} \to d$ (like hard). The transition is continuous and differentiable, avoiding the ringing of hard while reducing the bias of soft. This is the recommended default.

| Method | Bias | Ringing | Continuity at λ |
|:-------|:-----|:--------|:----------------|
| Hard | None | Yes | No |
| Soft | $\lambda$ | No | Yes |
| Semisoft | < $\lambda$ | No | Yes |

---

## 6. Offline vs Causal Behaviour

### 6.1 Offline

`denoise_offline_cpp` computes thresholds **once** from all $n / 2^\text{levels}$ finest-level detail coefficients of the full signal. This is the most accurate estimate: the full detail vector is available simultaneously, so the MAD reflects global noise characteristics.

Thresholds are computed in `compute_thresholds_internal` (defined in `offline.cpp`), applied in-place, then the inverse transform runs.

### 6.2 Causal Batch and Stream

The `WaveletEngine` recomputes thresholds from the **current window's** finest-level details every `update_freq` samples (`update_thresholds` in `WaveletEngine.h`). The estimate is derived from $W / 2$ coefficients — the detail subband at the finest decomposition level — fewer than the offline case, so it adapts to local noise changes but is more variable.

**`update_freq`** controls the update rate:
- `update_freq = 1` (default): thresholds recomputed at every sample. Maximum adaptivity; highest overhead.
- `update_freq = k`: recomputed every $k$ samples. Reduces overhead proportionally; acceptable if noise is quasi-stationary over $k$ samples.

For stationary noise, `update_freq = 10–50` is often sufficient. For rapidly changing noise (e.g. bursts in ECG), keep `update_freq = 1`.

The cached threshold vector (`current_lambdas`) persists between updates. If the window has not yet filled (`count < window_size`), the raw sample is returned as-is — thresholding starts only once the first full window is available.

---

## 7. Parameter Guidance

| Scenario | Recommended settings |
|:---------|:--------------------|
| Stationary Gaussian noise, smooth signal | α = 0.3, β = 1.2, semisoft |
| Non-stationary noise, real-time | α = 0.3, β = 1.0–1.2, semisoft, update_freq = 1 |
| Signals with sharp discontinuities | α = 0.3, β = 1.0–1.1, hard |
| Over-smoothed output (too much removed) | Decrease β toward 0.8–1.0 |
| Under-smoothed output (too much retained) | Increase β toward 1.5–2.0 |
| Flat threshold across levels desired | Increase α toward 5–10 |
