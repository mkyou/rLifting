# rLifting — Adding New Wavelets

rLifting is built around the `lifting_scheme` object as a first-class abstraction: every public function accepts one, and the C++ engine reads it at runtime without recompilation. This means any wavelet expressible as a sequence of predict and update steps can be used immediately across all operation modes, boundary extensions, and irregular-grid configurations.

---

## 1. Defining a Custom Wavelet

A wavelet is defined by two things: its lifting steps and its normalization.

### 1.1 Defining Steps

Each step is created with `lift_step`:

```r
lift_step(
  type,        # "predict" or "update"
  coeffs,      # numeric vector of filter coefficients
  start_idx,   # integer offset (or NULL to infer from position)
  position,    # "center", "left", "right" — used only when start_idx is NULL
  degree       # integer (or NULL to infer automatically)
)
```

The `start_idx` parameter determines where the filter window begins relative to the current position. For a predict step at odd index $i$, the filter reads even samples at indices $i + \text{start\_idx},\, i + \text{start\_idx} + 1,\, \ldots,\, i + \text{start\_idx} + k - 1$. When `start_idx = NULL`, it is inferred from `position`:
- `"center"`: window is centred on $i$ — `start_idx = -floor((k-1)/2)`.
- `"right"`: window starts at $i$ — `start_idx = 0`.
- `"left"`: window ends at $i$ — `start_idx = -(k-1)`.

### 1.2 Assembling the Scheme

Steps are combined into a `lifting_scheme` object via `custom_wavelet`:

```r
p <- lift_step("predict", coeffs = c(0.5, 0.5), start_idx = 0)
u <- lift_step("update",  coeffs = c(0.25, 0.25), start_idx = -1)
sch <- custom_wavelet("MyCDF53", list(p, u), normalization = c(sqrt(2), 1/sqrt(2)))
```

The resulting object is structurally identical to those produced by `lifting_scheme("cdf53")` and is accepted by all public functions without modification.

### 1.3 Normalization

The `normalization` vector `[norm_approx, norm_detail]` scales the approximation and detail subbands after each decomposition level. For orthonormal wavelets, both factors are chosen so that the transform is energy-preserving; for biorthogonal wavelets, they enforce a chosen amplitude convention. The normalization does not affect the filter coefficients or the irregular-grid behaviour. If no normalization is needed, pass `c(1, 1)`.

---

## 2. Automatic degree Inference

When a predict step has `sum(coeffs) ≈ 1` (within `1e-10`), it is recognized as an interpolating predictor and `degree` is set to `length(coeffs) - 1`. Otherwise `degree = -1`. This happens automatically in both `lift_step` and `lifting_scheme`; there is no need to specify it manually for standard designs.

The practical consequence:
- `coeffs = c(0.5, 0.5)` → `degree = 1` (linear interpolation) → irregular-grid Lagrange activated.
- `coeffs = c(-1/16, 9/16, 9/16, -1/16)` → `degree = 3` (cubic) → Lagrange activated.
- `coeffs = c(sqrt(3))` → `sum ≈ 1.73` → `degree = -1` → fixed coefficients, irregular grid ignored.

If a custom wavelet has a non-standard predict step that is nonetheless interpolating (e.g. a modified Lagrange design), `degree` can be set explicitly via the `degree` argument of `lift_step`. Setting `degree >= 0` on a non-interpolating step will activate Lagrange interpolation using physically incorrect weights — use only if the step genuinely interpolates.

---

## 3. Behaviour Across Operation Modes

A custom wavelet works identically in all four modes with no additional configuration.

**Step-by-step:** `lwt` and `ilwt` parse the steps at call time and apply them in order. Each predict step may call `apply_filter_cpp` (regular path) or `predict_irregular` (irregular path), depending on `degree` and whether `t` is supplied.

**Offline:** `denoise_offline_cpp` re-parses the scheme at the start of each call. The complete LWT → threshold → ILWT pipeline runs in C++ using the custom steps.

**Causal batch and stream:** The `WaveletEngine` constructor parses the scheme once and stores the compiled steps. All subsequent calls to `push_and_process` use the stored steps without re-parsing. There is no performance difference between a built-in wavelet and a custom one after construction.

In all modes, the number of decomposition levels must be compatible with the signal or window size: at each level, the approximation subband is halved, so `levels` must satisfy `n / 2^levels >= 1` (for offline) or `window_size / 2^levels >= 1` (for causal/stream).

---

## 4. Interaction with Boundary Extension

A custom wavelet interacts with boundary modes exactly as built-in wavelets do: the boundary mode is applied to out-of-bounds reads in predict and update steps, and it is independent of the wavelet design.

The only constraint worth noting is filter support: a wavelet with $k$ predict coefficients generates $k - 1$ out-of-bounds reads per edge per decomposition level. For large $k$ (e.g. $k = 8$), the boundary zone grows substantially with depth — up to $(k-1) \cdot 2^\text{levels}$ samples from each edge at the deepest level. The choice of boundary mode matters more for wide-support wavelets.

For causal/stream mode, `one_sided` is the most robust boundary choice for any wavelet (see `03-boundary-modes.md`, section 4.2), since it makes no assumption about the signal beyond the window edge.

---

## 5. Interaction with Irregular Grids

A custom wavelet automatically benefits from irregular-grid Lagrange interpolation if its predict step has `degree >= 0`. No additional code is required: passing `t` to any public function activates the position-aware path for steps where `degree >= 0` and uses fixed coefficients for steps where `degree = -1`.

If the wavelet has multiple predict steps with different degrees (as in complex multi-step designs), each step is handled independently according to its own `degree`. Update steps always use fixed coefficients regardless of `degree`.

**When to design for irregular grids:** if the intended use is signal processing on physically non-uniform samples (e.g. spatial sensors, variable-rate recording), design the predict step to be interpolating (`sum(coeffs) ≈ 1`). This ensures the Lagrange correction is activated and detail coefficients of smooth signals remain near zero on irregular grids. Non-interpolating custom wavelets (orthogonal, biorthogonal) will silently use fixed coefficients when `t` is supplied, with the same accuracy loss as `db2` or `cdf97` on non-uniform grids.

---

## 6. Interaction with Adaptive Thresholding

The adaptive thresholding pipeline is wavelet-agnostic: it reads the finest-level detail coefficients, estimates noise via MAD, and applies the recursive threshold formula. The wavelet design affects thresholding only indirectly, through the energy distribution of the detail coefficients:

- **Vanishing moments**: a wavelet with more vanishing moments concentrates smooth-signal energy in the approximation subband, leaving detail coefficients closer to zero for smooth signals. This makes the MAD estimator more accurate, because the detail vector is more noise-dominated.
- **Filter support**: wider filters produce a larger boundary zone, which may add non-zero detail coefficients near the edges. For short signals or deep decomposition, this can inflate the MAD estimate slightly.

Neither of these properties requires any change to the thresholding code. The parameters `alpha`, `beta`, `method`, and `update_freq` apply identically to custom wavelets.

---

## 7. Validation

Before deploying a custom wavelet, run `diagnose_wavelet`. It accepts a `lifting_scheme` object (or a built-in wavelet name) and a `config` list:

```r
config <- list(
  is_ortho   = FALSE,   # expected orthogonality (TRUE for haar, db2)
  vm_degrees = 0:1,     # polynomial degrees to test for vanishing moments
  max_taps   = 6        # maximum expected non-zero taps in the impulse response
)
diag <- diagnose_wavelet(sch, config, verbose = TRUE, plot = TRUE)
```

The function runs five checks in sequence and returns an S3 object of class `wavelet_diagnosis`. Setting `verbose = TRUE` prints a formatted table; `plot = TRUE` renders the wavelet ($\psi$) and scaling ($\phi$) basis functions via cascade reconstruction.

### 7.1 Perfect Reconstruction

Tests $\| x - \text{ilwt}(\text{lwt}(x)) \|_\infty < 10^{-9}$ on six signals: random noise, ramp, sine, Doppler, HeaviSine, and bumps (all length 512, `periodic` extension). Periodic extension is used to isolate the mathematical properties of the transform from boundary effects.

**Failure diagnosis.** A failed perfect reconstruction almost always means one of:
- A sign error in a predict or update coefficient (e.g. `+` vs `−` in the update step);
- An incorrect `start_idx` that shifts the filter window and breaks the symmetry needed for exact inversion;
- A missing or incorrect normalization factor (e.g. `c(sqrt(2), 1/sqrt(2))` instead of `c(1, 1)`).

### 7.2 Orthogonality (Energy Conservation)

Verifies Parseval's theorem: $E_\text{out} / E_\text{in} \approx 1$ on a random signal (length 512, `periodic` extension), where $E_\text{out} = \sum a_1^2 + \sum d_1^2$ at level 1. The test passes if the computed ratio matches `config$is_ortho`:
- `is_ortho = TRUE`: expects ratio close to 1.
- `is_ortho = FALSE`: the ratio is reported as informational; the test always passes.

Only Haar and DB2 are orthogonal among the built-in wavelets. CDF 5/3, CDF 9/7, and DD4 are biorthogonal — their energy ratio is typically not 1, but they still have perfect reconstruction.

### 7.3 Vanishing Moments

For each degree $d$ in `config$vm_degrees`, tests whether the finest-level detail coefficients of a degree-$d$ polynomial are near zero. The polynomial signals are: constant (degree 0), ramp (degree 1), parabola (degree 2), cubic (degree 3).

The energy check excludes the outer 15% of the detail vector (the boundary zone) to avoid conflating boundary artefacts with a genuine failure:

```r
cut    <- floor(n * 0.15)
d_core <- d1[(cut + 1):(length(d1) - cut)]
energy <- sum(d_core^2)   # must be < 1e-9
```

A wavelet with $p$ vanishing moments exactly annihilates polynomials up to degree $p - 1$. For example, CDF 5/3 has $p = 2$ (annihilates constants and ramps). Pass `vm_degrees = 0:(p-1)` to verify the expected count.

**Failure diagnosis.** Non-zero energy for a polynomial of degree $d$ means the predict step does not reproduce that polynomial — either the coefficients don't sum to the right value, or the `start_idx` creates an asymmetric window that breaks polynomial reproduction.

### 7.4 Compact Support

Measures the number of non-zero taps in the wavelet's impulse response. Internally, it places a unit impulse in the finest-level detail coefficients of a synthetic `lwt` object and reconstructs via `ilwt` with `zero` extension. The number of non-zero values in the output (at tolerance `1e-10`) must be within `[1, config$max_taps + 2]`.

The `+2` margin accounts for minor floating-point spill at filter boundaries. A value of 0 would mean the wavelet produces no response at all (degenerate scheme); a value exceeding `max_taps + 2` means the effective filter support is wider than expected — typically caused by a large `start_idx` offset.

### 7.5 Shift Sensitivity

Measures the variation in finest-level detail energy when the input is shifted by one sample:

$$\Delta = \frac{|E(x) - E(x_{\text{shift}})|}{E(x)} \times 100\%$$

This check always passes (it is informational, not a pass/fail test). The result quantifies the **shift variance** of the wavelet: orthogonal decimated wavelets are not translation-invariant, so $\Delta > 0$ is expected. A value near 0% indicates the wavelet is close to shift-invariant for the test signal; a value near 100% indicates strong sensitivity.

This metric is useful for comparing wavelets when shift invariance matters (e.g. edge detection applications).

### 7.6 Typical Config by Wavelet

| Wavelet | `is_ortho` | `vm_degrees` | `max_taps` |
|:--------|:-----------|:-------------|:-----------|
| `haar` | `TRUE` | `0` | 2 |
| `db2` | `TRUE` | `0:1` | 4 |
| `cdf53` | `FALSE` | `0:1` | 4 |
| `dd4` | `FALSE` | `0:3` | 8 |
| `cdf97` | `FALSE` | `0:3` | 12 |

For a custom wavelet, `max_taps` should equal the sum of predict and update filter lengths (or slightly more if chained steps overlap in support).

---

## 8. Example: A Custom CDF 5/3 Variant

```r
library(rLifting)

# Standard CDF 5/3
p <- lift_step("predict", c(0.5, 0.5), start_idx = 0)
u <- lift_step("update",  c(0.25, 0.25), start_idx = -1)
sch <- custom_wavelet("MyCDF53", list(p, u), c(sqrt(2), 1/sqrt(2)))

# Works in all modes
x      <- rnorm(256)
res    <- denoise_signal_offline(x, sch, levels = 3)
stream <- new_wavelet_stream(sch, window_size = 63, levels = 3)

# Irregular grid — degree = 1 inferred automatically
t <- cumsum(c(0, abs(rnorm(255, mean = 1, sd = 0.3))))
res_irr <- denoise_signal_offline(x, sch, levels = 3, t = t)

# Validate
config <- list(is_ortho = FALSE, vm_degrees = 0:1, max_taps = 4)
diag   <- diagnose_wavelet(sch, config, verbose = TRUE, plot = FALSE)
print(diag)
```
