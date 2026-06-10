# rLifting — Adaptive Thresholding (Technical Reference)

Technical specification of the threshold subsystem: estimators, formulas, code paths, edge cases. User-facing walkthrough with worked examples and decision tables lives in `vignette("02-thresholding-and-tuning")`. Architectural index: [`00-design-overview.md`](00-design-overview.md). Related notes: [`01-lifting-scheme-and-transform.md`](01-lifting-scheme-and-transform.md) (where the coefficients come from), [`03-zero-allocation-engine.md`](03-zero-allocation-engine.md) (the per-sample hot path), [`04-boundary-and-threshold.md`](04-boundary-and-threshold.md) (MAD-via-`nth_element` and shrinkage kernels).

Symbols used throughout:

- $y \in \mathbb{R}^n$ — observed signal; assumed $y = f + \varepsilon$ with $\varepsilon_i \sim \mathcal{N}(0,\sigma^2)$ i.i.d. unless otherwise noted.
- $d_k \in \mathbb{R}^{m_k}$ — detail subband at decomposition level $k = 1, \dots, K$, with $m_1 \approx n/2$ and $m_k \approx n/2^k$.
- $\hat\sigma$ — global noise-scale estimate from $d_1$ (universal rule); $\hat\sigma_k$ — per-level estimate (SURE rule).
- $\lambda_k$ — threshold applied to subband $d_k$. Shrinkage function $\eta_\lambda: \mathbb{R} \to \mathbb{R}$ acts coefficient-wise: $\hat d_{k,i} = \eta_{\lambda_k}(d_{k,i})$.
- $\alpha,\beta$ — recursion parameters of the universal rule ($\beta$ scales $\lambda_1$; $\alpha$ controls level-to-level decay).

---

## 1. Noise Estimation: MAD with the $0.6745$ Constant

For a Gaussian-mixture detail subband with most mass from noise and a sparse tail from signal, the **Median Absolute Deviation** is preferred to the empirical SD because order statistics are insensitive to a small fraction of large $|d_i|$ from real signal energy. The estimator used is

$$\hat\sigma = \frac{\operatorname{med}_i |d_{1,i}|}{0.6745}, \qquad 0.6745 \approx \Phi^{-1}(3/4),$$

where $\Phi$ is the standard-normal CDF. Under $d_{1,i} \sim \mathcal{N}(0,\sigma^2)$, $|d_{1,i}|$ has its median at $\sigma \cdot \Phi^{-1}(3/4)$, so dividing by $0.6745$ yields a Fisher-consistent estimator of $\sigma$.

The note assumes the standard form $\operatorname{MAD}(d) = \operatorname{med}(|d|)$ (no recentring by the median of $d$). This is valid because $d_1$ is mean-zero by construction: after lifting, predict and update steps remove the polynomial trend that survives in the approximation subband. Skipping the recentring saves one full pass over $d_1$.

**Implementation.** A single helper `compute_mad` ([`inst/include/rLifting/utils.h`](../../inst/include/rLifting/utils.h)) implements the canonical averaged median via `std::nth_element` (O(n) selection, no sort): one call to find the upper middle, a second on the truncated range to find the lower middle, then their mean for even $n$. All four call sites route through it:

- `compute_thresholds_cpp` / `compute_mad_cpp` ([`src/adaptative.cpp`](../../src/adaptative.cpp), R-exposed).
- `WaveletEngine::update_thresholds` ([`inst/include/rLifting/WaveletEngine.h`](../../inst/include/rLifting/WaveletEngine.h)) — causal/stream hot path.
- `compute_thresholds_internal` ([`src/offline.cpp`](../../src/offline.cpp)) — offline universal rule.
- `compute_sure_lambda_level` ([`inst/include/rLifting/utils.h`](../../inst/include/rLifting/utils.h)) — SURE per-level.

Routing through one helper guarantees `compute_adaptive_threshold()`, `denoise_signal_offline()`, `denoise_signal_causal()`, and `new_wavelet_stream()` produce the same $\hat\sigma$ on the same $d_1$.

See [`04-boundary-and-threshold.md`](04-boundary-and-threshold.md) for the algorithmic discussion of `nth_element` versus a full sort and why it matters in the per-sample causal hot path.

---

## 2. Universal Rule: VisuShrink with $\alpha/\beta$ Recursion

**Per-level threshold.** Selected by `threshold_method = "universal"` (default). The finest level uses the Donoho–Johnstone universal threshold scaled by $\beta$:

$$\lambda_1 = \beta \cdot \hat\sigma \cdot \sqrt{2 \log m_1},$$

where $m_1$ is the length of $d_1$. The $\sqrt{2 \log m_1}$ factor is the asymptotic upper envelope of $m_1$ i.i.d. half-normal variates: for pure noise, $\mathbb{P}\bigl(\max_i |d_{1,i}|/\hat\sigma > \sqrt{2 \log m_1}\bigr) \to 0$, so $\beta = 1$ kills essentially all noise-only coefficients.

**Recursion across levels** (Liu, Mi & Mao, 2014):

$$\lambda_k = \lambda_{k-1} \cdot \frac{k - 1}{k + \alpha - 1}, \qquad k = 2, \dots, K.$$

Intuition for $\alpha$:

- $\alpha = 0$: $\lambda_k = \lambda_{k-1}$, flat across levels. Useful when coarse subbands still contain non-negligible noise.
- $\alpha \to \infty$: $\lambda_k \to 0$ for $k \geq 2$. Coarse subbands pass through nearly unthresholded — appropriate when coarse coefficients are signal-dominated and you want to avoid biasing them.
- Default $\alpha = 0.3$ gives $\lambda_k / \lambda_1 \approx \{1,\,0.77,\,0.67,\,0.61,\,0.57\}$ for $k = 1,\dots,5$ — moderate decay.

Intuition for $\beta$: pure multiplicative gain on $\lambda_1$. $\beta < 1$ retains more coefficients (less smoothing, higher variance); $\beta > 1$ kills more (more smoothing, more bias). Default $\beta = 1.2$ slightly overshoots VisuShrink to compensate for the dependence in lifting-domain coefficients (they are not exactly i.i.d. Gaussian even under white noise).

**Code.** [`compute_thresholds_cpp(d1, max_level, alpha, beta)` in `src/adaptative.cpp`](../../src/adaptative.cpp), returned to R as `adaptive_thresholds` by [`R/adaptative_threshold.R`](../../R/adaptative_threshold.R). Inline duplicates in [`compute_thresholds_internal` in `src/offline.cpp`](../../src/offline.cpp) and [`WaveletEngine::update_thresholds` in `WaveletEngine.h`](../../inst/include/rLifting/WaveletEngine.h) — the universal branch is the same formula in all three locations.

---

## 3. SURE Rule: Per-Level Stein Risk

Selected by `threshold_method = "sure"`. Each subband is treated independently. For subband $k$ with coefficients $d_{k,1}, \dots, d_{k,m_k}$ and a per-level scale $\hat\sigma_k = \operatorname{med}(|d_k|)/0.6745$, Stein's Unbiased Risk Estimate for the **soft-threshold** estimator is

$$\operatorname{SURE}(\lambda;\,d_k,\hat\sigma_k) \;=\; m_k\,\hat\sigma_k^2 \;+\; \sum_{i=1}^{m_k} \min(d_{k,i}^2,\,\lambda^2) \;-\; 2\,\hat\sigma_k^2 \cdot \#\{i : |d_{k,i}| \leq \lambda\}.$$

This is an unbiased estimator of the mean squared $\ell_2$ risk of the soft-threshold output (Stein, 1981; Donoho & Johnstone, 1995):

$$\mathbb{E}\bigl[\operatorname{SURE}(\lambda;d_k,\sigma)\bigr] = \mathbb{E}\bigl\lVert \eta^{\text{soft}}_\lambda(d_k) - f_k \bigr\rVert_2^2,$$

so minimising it over $\lambda$ yields the data-adaptive risk-minimising soft threshold for that subband.

**Candidate set.** The risk is piecewise quadratic in $\lambda$ between consecutive order statistics of $|d_k|$, so the minimum lies at one of those order statistics. The implementation sorts $|d_k|$ ascending and evaluates SURE at each of the $m_k$ candidates:

$$\lambda_k^{\text{SURE}} = \arg\min_{\lambda \in \{|d_{k,(1)}|,\dots,|d_{k,(m_k)}|\}} \operatorname{SURE}(\lambda;d_k,\hat\sigma_k).$$

**Cap.** The result is clamped by the universal threshold computed with the **per-level** $\hat\sigma_k$, *not* by the global $\hat\sigma$ from $d_1$:

$$\lambda_k^{\text{final}} = \min\bigl(\lambda_k^{\text{SURE}},\; \hat\sigma_k\sqrt{2\log m_k}\bigr).$$

The cap (without the $\beta$ multiplier) guards against the degenerate case where SURE prefers $\lambda \approx \max |d_k|$ in a signal-dominated subband, which would shrink nearly all coefficients.

**Code.** [`compute_sure_lambda_level` in `inst/include/rLifting/utils.h`](../../inst/include/rLifting/utils.h). The implementation:

1. Computes $\hat\sigma_k$ via `compute_mad` on $|d_k|$ (the shared MAD helper).
2. Performs a full `std::sort` of `abs_d` — required because SURE is evaluated at *every* order statistic in turn, not just the median.
3. Builds the cumulative sum $\sum_{i \leq k} d_{(i)}^2$ in linear time. Using sorted $|d|$, for $\lambda = |d_{(k+1)}|$ the SURE expression reduces to `cumsum_sq[k+1] + λ²·(n-k-1) + n·σ² - 2σ²·(k+1)`, giving $O(n)$ total cost for the candidate sweep after the $O(n \log n)$ sort.
4. Clamps by `sigma * sqrt(2 log n)`.

The SURE branch is wired through [`compute_thresholds_sure_internal` in `src/offline.cpp`](../../src/offline.cpp) for offline use and the `threshold_method == "sure"` branch of [`update_thresholds` in `WaveletEngine.h`](../../inst/include/rLifting/WaveletEngine.h) for causal/stream use.

**Effect of $\alpha,\beta$ in SURE mode.** None. The recursion is bypassed; each $\lambda_k$ is computed independently. The R wrappers still accept `alpha` and `beta`, but they are not consumed by the SURE branch — a no-op for compatibility.

**When to prefer which rule.**

- **Universal**: stationary white noise; smooth or oscillatory signals where the coefficient distribution at each level is close to Gaussian. Cheap, stable, theory-backed.
- **SURE**: signals with sharp transitions whose coefficient distributions at fine levels are heavy-tailed; non-stationary noise variance across subbands. The per-level $\hat\sigma_k$ tracks subband variance independently.

---

## 4. Shrinkage Functions

All four operate coefficient-wise on $d \in \mathbb{R}$ and zero $d$ whenever $|d| < \lambda$. They differ in the active region $|d| \geq \lambda$.

### 4.1 Hard

$$\eta^{\text{hard}}_\lambda(d) = d \cdot \mathbb{1}\{|d| \geq \lambda\}.$$

Bias-free for $|d| \geq \lambda$; discontinuous at $\pm\lambda$ — yields Gibbs-like ringing near discontinuities. Code: [`threshold_hard_cpp` in `src/thresholding.cpp`](../../src/thresholding.cpp).

### 4.2 Soft

$$\eta^{\text{soft}}_\lambda(d) = \operatorname{sign}(d) \cdot \max(|d| - \lambda,\,0).$$

Continuous; biased downward by exactly $\lambda$ for every active coefficient. SURE is unbiased only for the soft estimator, which is why `tune_alpha_beta()` uses it as the optimisation surrogate even when the user picks a different shrinkage at denoise time. Code: [`threshold_soft_cpp` in `src/thresholding.cpp`](../../src/thresholding.cpp).

### 4.3 Semisoft (hyperbolic; Liu et al., 2014)

$$\eta^{\text{semisoft}}_\lambda(d) = \operatorname{sign}(d) \cdot \sqrt{\max(d^2 - \lambda^2,\,0)}.$$

Continuous at $\lambda$, asymptotically the identity as $|d| \to \infty$ (bias $\to 0$). Intermediate between hard and soft; package default. Code: [`threshold_semisoft_cpp` in `src/thresholding.cpp`](../../src/thresholding.cpp).

### 4.4 SCAD (Antoniadis & Fan, 2001)

Three-region rule with continuity at $\lambda$, $2\lambda$, and $a\lambda$:

$$\eta^{\text{scad}}_{\lambda,a}(d) = \begin{cases}
0 & |d| \leq \lambda, \\[2pt]
\operatorname{sign}(d)\,(|d| - \lambda) & \lambda < |d| \leq 2\lambda, \\[2pt]
\dfrac{(a-1)\,d \;-\; a\lambda\,\operatorname{sign}(d)}{a - 2} & 2\lambda < |d| \leq a\lambda, \\[2pt]
d & |d| > a\lambda.
\end{cases}$$

The shape parameter $a > 2$ controls the width of the transition region. The canonical Fan–Li value $a = 3.7$ is derived from a Bayes-risk analysis under Gaussian noise. As $a \downarrow 2$, SCAD approaches soft thresholding (the linear segment vanishes). As $a \to \infty$, the linear segment dominates and SCAD approaches soft for moderate $|d|$ but never reaches hard's behaviour; SCAD interpolates *non-trivially* between soft and hard rather than collapsing onto either limit. Above $a\lambda$ SCAD is the identity, recovering hard's bias-freeness for the largest coefficients while keeping continuity. Code: [`threshold_scad_cpp` in `src/thresholding.cpp`](../../src/thresholding.cpp); reused inline in [`push_and_process` in `WaveletEngine.h`](../../inst/include/rLifting/WaveletEngine.h) for the causal/stream paths and inlined again in [`denoise_offline_cpp` in `src/offline.cpp`](../../src/offline.cpp) for the offline path. The shape parameter `a` is honoured at every entry point: `threshold_scad()` (standalone), `denoise_signal_offline()` (via the `scad_a` parameter of `denoise_offline_cpp`), and `denoise_signal_causal()` / `new_wavelet_stream()` (via the `scad_a` member of `WaveletEngine`, set at construction by `create_engine_cpp` / `run_causal_batch_cpp`).

### 4.5 Comparison

| Rule | Bias on large $d$ | Continuity at $\lambda$ | Continuity at $\infty$ | Extra param |
|:-----|:------------------|:------------------------|:-----------------------|:------------|
| Hard | $0$ | No | — (identity) | — |
| Soft | $\lambda$ | Yes | $-\lambda$ asymptote | — |
| Semisoft | $O(\lambda^2/|d|)$ | Yes | Identity asymptote | — |
| SCAD | $0$ for $|d| > a\lambda$ | Yes | Identity (exact above $a\lambda$) | $a$ (default 3.7) |

The `threshold()` dispatcher ([`R/thresholding.R`](../../R/thresholding.R)) is a thin `switch` over the four C++ kernels; same names accepted by `denoise_signal_offline()`, `denoise_signal_causal()`, and `new_wavelet_stream()`.

---

## 5. Offline vs Causal Threshold Updates

### 5.1 Offline (`denoise_offline_cpp`)

Full signal in hand; thresholds computed **once** from the full $d_1$ vector (universal) or from each full $d_k$ (SURE). Maximal sample size for the MAD estimator. Branch: `threshold_method == "sure"` inside [`denoise_offline_cpp` in `src/offline.cpp`](../../src/offline.cpp). Result is applied in place; ILWT runs immediately after.

### 5.2 Causal (`denoise_signal_causal`, `new_wavelet_stream`)

Thresholds are state, not output. The [`WaveletEngine` in `WaveletEngine.h`](../../inst/include/rLifting/WaveletEngine.h) caches the current threshold vector in the `current_lambdas` member and refreshes it according to:

```cpp
bool should_update = (update_freq > 0)
                        ? (step_iter % update_freq == 0)
                        : !lambdas_initialized;
if (should_update) {
  update_thresholds(alpha, beta, threshold_method);
  lambdas_initialized = true;
}
```

`update_freq <= 0` (validated and warned in `R/realtime_denoising.R`) means "freeze after warm-up" — the `lambdas_initialized` flag gates the single initial update; subsequent samples reuse the cached threshold vector forever. Between updates the cached $\lambda_k$ are reused unchanged — the only cost per sample is the shrinkage scan over `work_detail[j]`.

The MAD inside `update_thresholds` sees only $m_1 = W/2$ coefficients ($W$ = `window_size`), the finest-level details of the current sliding window. This is substantially fewer than offline mode; the estimate is noisier but tracks local non-stationarity.

**Warm-up.** The engine returns the raw sample as-is until the ring buffer is full (early-return guard at the top of [`push_and_process` in `WaveletEngine.h`](../../inst/include/rLifting/WaveletEngine.h)):

```
if (count < window_size) return new_val;
```

No transform, no thresholding, no allocation. Filtering starts at sample $W$.

**`update_freq` trade-off.**

- `update_freq = 1` (default): one MAD + recursion per sample. Maximum adaptivity; per-sample cost is dominated by the threshold update for small $W$.
- `update_freq = k`: cost amortised over $k$ samples; the cached $\lambda_k$ vector lags the true noise by up to $k$ samples. Acceptable for noise stationary over the update interval.
- Practical range: 1 for non-stationary noise (e.g. ECG bursts); 10–50 for stationary or slowly drifting noise; values much larger than $W$ defeat the point of streaming.

Both rules are supported in causal mode: the SURE branch evaluates per-level SURE on the window's `work_detail[j]` and ignores $\alpha,\beta$; the universal branch runs MAD on `work_detail[0]` and applies the recursion.

See [`03-zero-allocation-engine.md`](03-zero-allocation-engine.md) for the broader hot-path design (ring buffer, XPtr finalizer, allocation accounting); `update_thresholds` is the only piece of the per-sample loop that may allocate, and only when MAD's working vector $|d_1|$ is touched — which is why it is guarded by `update_freq`.

---

## 6. `tune_alpha_beta()` — Joint $(\alpha,\beta)$ SURE Minimisation

Source: [`R/tune_alpha_beta.R`](../../R/tune_alpha_beta.R). Selects $(\alpha,\beta)$ for the **universal** rule by minimising SURE for soft thresholding on the full set of detail subbands. Acts on $(\alpha,\beta)$ **jointly**, not marginally.

**Objective.** Given the LWT details $\{d_k\}_{k=1}^K$ of the input signal, let $\lambda_k(\alpha,\beta)$ be the universal recursion of Section 2. The objective is

$$\Psi(\alpha,\beta) \;=\; \sum_{k=1}^{K}\Bigl[\,m_k\,\hat\sigma_k^2 \;+\; \sum_{i=1}^{m_k}\min(d_{k,i}^2,\lambda_k^2) \;-\; 2\,\hat\sigma_k^2\,\#\{i: |d_{k,i}| \leq \lambda_k\}\,\Bigr].$$

Each subband's risk is evaluated with its own MAD-based $\hat\sigma_k$; the recursion itself is driven by $\hat\sigma = \hat\sigma_1$ from $d_1$. This decoupling is deliberate: the recursion sets the threshold *location* from the global noise scale, while the risk is measured against each subband's actual noise.

**Domain.** A box $[\alpha_{\min},\alpha_{\max}] \times [\beta_{\min},\beta_{\max}]$. Defaults `alpha_range = c(0, 10)`, `beta_range = c(0.5, 3.0)` — chosen so the box brackets the practically useful corner of parameter space (very weak to very strong decay; mild under-thresholding to aggressive over-thresholding).

**Optimiser.** Two phases:

1. **Coarse grid**: $21 \times 21$ uniform grid over the box. $\Psi$ is piecewise-constant in $\beta$ at scales below an order-statistic gap of $|d_k|$ and piecewise-smooth in $\alpha$; the grid handles the resulting non-smooth landscape, which would mislead a pure-gradient optimiser. The grid minimum becomes the warm-start.
2. **Nelder-Mead refinement**: `stats::optim(method = "Nelder-Mead", reltol = 1e-8)` from the grid optimum. Nelder-Mead is unconstrained, so the refined point is clipped to the box. If the clipped refinement does not improve over the grid minimum, the grid minimum is kept (clip-and-fallback block in [`tune_alpha_beta` in `R/tune_alpha_beta.R`](../../R/tune_alpha_beta.R)). The `converged` flag reports `opt$convergence == 0L`.

**Return value.** `list(alpha, beta, sure, converged)`. The `sure` field is the minimised $\Psi$ (not a noise-scale-free risk).

**Portability across shrinkages.** SURE is unbiased *only* for soft thresholding. The optimal $(\alpha,\beta)$ are nonetheless usable with `shrinkage = "semisoft"|"hard"|"scad"` because all four shrinkages share the same threshold *location* $\lambda_k$; the optimised threshold targets the right magnitude even if the surrogate risk is slightly biased for the actual estimator. Empirically the carry-over is good (see vignette Section 5).

**Tuning the tuner.**

- Low-SNR signals: narrow `beta_range` upward (e.g. `c(1.0, 3.0)`) — the optimum is rarely below 1.
- Signals dominated by sharp edges: narrow `alpha_range` toward `c(0, 2)` — heavy decay over-attenuates the coarse-level transients.
- Cost: one LWT + $441$ SURE evaluations + a Nelder-Mead trajectory; for $n \sim 10^3$ this is order of seconds.

---

## 7. Edge Cases

- **Empty $d_1$.** `compute_adaptive_threshold` returns `list(d1 = 0)` with class `adaptive_thresholds` ([`compute_adaptive_threshold` in `R/adaptative_threshold.R`](../../R/adaptative_threshold.R)). The C++ paths return a zero vector when $\hat\sigma < 10^{-15}$ — guards present in `compute_thresholds_cpp` (`src/adaptative.cpp`), `update_thresholds` (`WaveletEngine.h`), `compute_thresholds_internal` (`src/offline.cpp`), and `compute_sure_lambda_level` (`utils.h`). Downstream, all coefficients pass through unchanged.
- **$n_1$ very small** (e.g. `levels` chosen so $m_1 \leq 4$). $\sqrt{2 \log m_1}$ is poorly calibrated and MAD has high variance; consider raising `update_freq` and `window_size` in causal mode or reducing `levels` offline.
- **SCAD with $a \leq 2$.** Rejected at the kernel boundary (`stop(...)` in [`threshold_scad_cpp` in `src/thresholding.cpp`](../../src/thresholding.cpp)). Neither the engine nor `denoise_offline_cpp` currently validate `a` at construction/entry; values $\leq 2$ produce a degenerate `scad_denom = a - 2` (zero or negative) and silently break the SCAD region-3 formula. R wrappers should clamp or validate.
- **`update_freq = 0`** is accepted with a warning and treated as "freeze after warm-up": thresholds are computed once on the first full window and held fixed. Negative values are rejected by [`R/realtime_denoising.R`](../../R/realtime_denoising.R) (`.validate_update_freq`). The C++ engine carries a `lambdas_initialized` flag to gate the single update.

---

## References

Antoniadis, A., & Fan, J. (2001). Regularization of wavelet approximations. *Journal of the American Statistical Association*, **96**(455), 939–967.

Donoho, D. L., & Johnstone, I. M. (1994). Ideal spatial adaptation by wavelet shrinkage. *Biometrika*, **81**(3), 425–455.

Donoho, D. L., & Johnstone, I. M. (1995). Adapting to unknown smoothness via wavelet shrinkage. *Journal of the American Statistical Association*, **90**(432), 1200–1224.

Liu, Z., Mi, Y., & Mao, Y. (2014). Improved real-time denoising method based on lifting wavelet transform. *Measurement Science Review*, **14**(3), 152–159. DOI: 10.2478/msr-2014-0020.

Stein, C. M. (1981). Estimation of the mean of a multivariate normal distribution. *Annals of Statistics*, **9**(6), 1135–1151.
