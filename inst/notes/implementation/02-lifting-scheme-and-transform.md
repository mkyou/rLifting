# rLifting — Lifting Scheme and Transform

This document describes the `lifting_scheme` object — how it is represented in R, how it crosses the R/C++ boundary, and how it drives the forward and inverse transforms. It also explains why the transform is implemented in multiple separate loops rather than a single shared routine, and how the irregular-grid path is embedded inside the predict step.

---

## 1. The `lifting_scheme` Object

The central abstraction in rLifting is the `lifting_scheme` S3 object. Every public function accepts one, and its `steps` field is consumed directly by C++ without further transformation. The object contains three fields:

- `wavelet`: a string identifier (for display and diagnostic purposes only; it does not affect computation).
- `steps`: a list of predict (P) and update (U) steps, each represented as a named R list with fields `type` (`"predict"` or `"update"`), `coeffs` (numeric vector), `start_idx` (integer offset), and `degree` (integer).
- `normalization`: a length-2 numeric vector `[norm_approx, norm_detail]`.

The choice to represent the scheme as a plain R list — rather than, say, a compiled filter object — was deliberate. It keeps the scheme fully inspectable and modifiable from R, allows users to define custom wavelets by constructing lists directly, and maps with zero overhead to the C++ struct that consumes it.

### 1.1 Degree Inference

The `degree` field is the only non-obvious element of the step format. It governs whether a predict step uses fixed coefficients or position-aware Lagrange interpolation on irregular grids. In the R constructor (`lifting_scheme` and `lift_step`), it is inferred automatically with a single test:

```r
degree = if (type == "predict" && abs(sum(coeffs) - 1) < 1e-10)
    as.integer(length(coeffs) - 1L)
else
    -1L
```

The test `abs(sum(coeffs) - 1) < 1e-10` identifies interpolating predict steps: those whose coefficients form an interpolating polynomial (the predicted value at the target position equals the polynomial exactly). For such steps, the degree of the polynomial equals the number of coefficients minus one. For non-interpolating steps (updates, or predict steps of orthogonal wavelets such as DB2), `degree` is set to `-1`, which disables the irregular path entirely.

This design means the user never specifies `degree` explicitly for standard wavelets; it is a derived property of the filter design.

### 1.2 The Built-in Wavelets

The six built-in wavelets illustrate the range of step configurations:

| Wavelet | Steps | Normalization | Notes |
|:--------|:------|:-------------|:------|
| `lazy` | None | `[1, 1]` | Identity split; no P or U |
| `haar` | P: `[1]`, U: `[0.5]` | `[√2, 1/√2]` | Nearest-neighbour predict |
| `cdf53` | P: `[0.5, 0.5]`, U: `[0.25, 0.25]` | `[√2, 1/√2]` | Linear interpolation predict |
| `dd4` | P: `[-1/16, 9/16, 9/16, -1/16]`, U: `p/2` | `[√2, 1/√2]` | Cubic Lagrange predict |
| `db2` | P: `[√3]`, U: `[√3/4, (√3−2)/4]`, P: `[-1]` | `[(√3+1)/√2, (√3−1)/√2]` | Three-step factorisation; orthogonal |
| `cdf97` | P/U/P/U (four steps) | `[ζ, 1/ζ]` | JPEG 2000 wavelet; biorthogonal |

DB2 and CDF 9/7 have `sum(predict_coeffs) ≠ 1`, so `degree = -1` for all their predict steps. Their factorizations follow Daubechies and Sweldens (1998) and were chosen to maximise vanishing moments given a fixed filter length — not to produce an interpolating predictor.

---

## 2. Crossing the R/C++ Boundary

When a public function calls a C++ routine, the `lifting_scheme` object's `steps` field is passed as an `Rcpp::List`. Each C++ function that executes the transform must parse this list into a `std::vector<LiftingStep>`, where `LiftingStep` is defined in `inst/include/rLifting/utils.h`:

```cpp
struct LiftingStep {
    std::string type;
    std::vector<double> coeffs;
    int start_idx;
    int degree = -1;
};
```

The parsing loop is:

```cpp
for (int i = 0; i < n_steps; i++) {
    List s = steps[i];
    LiftingStep step;
    step.type      = as<std::string>(s["type"]);
    step.coeffs    = as<std::vector<double>>(s["coeffs"]);
    step.start_idx = s["start_idx"];
    step.degree    = s.containsElementNamed("degree") ? (int)s["degree"] : -1;
    cpp_steps.push_back(step);
}
```

This parsing runs at the start of `lwt_cpp`, `ilwt_cpp`, `denoise_offline_cpp`, and the `WaveletEngine` constructor — four independent sites. The duplication is a deliberate trade-off: centralising the parsing would require passing the parsed vector across the R/C++ boundary or through a shared pointer, which introduces its own overhead and complexity. Since parsing happens once per transform call (or once at engine construction), the cost is negligible. The maintenance cost is acknowledged: any change to the step format must be applied in all four locations.

---

## 3. Two Implementations of the Transform Loop

The transform is implemented in three separate C++ locations:

1. **`lwt_cpp` / `ilwt_cpp`** (`src/lwt.cpp`, `src/ilwt.cpp`) — used by the step-by-step R API (`lwt`, `ilwt`). These call `apply_filter_cpp`, which operates on `Rcpp::NumericVector`.
2. **`denoise_offline_cpp`** (`src/offline.cpp`) — used by `denoise_signal_offline`. Contains its own inline lifting loops operating on `std::vector<double>`.
3. **`WaveletEngine::push_and_process`** (`inst/include/rLifting/WaveletEngine.h`) — used by causal and stream modes. Also uses inline loops on `std::vector<double>`.

The reason for two separate inline implementations (offline and engine) rather than sharing a single loop is that the data structures differ fundamentally: `denoise_offline_cpp` works on locally allocated vectors of varying size (one per level), while `WaveletEngine` works on pre-allocated workspace vectors of fixed size. Abstracting over both would require either dynamic dispatch or template parameters — neither is justified for a fixed set of two usage patterns.

The reason `lwt_cpp` delegates to `apply_filter_cpp` (instead of inlining) is that `apply_filter_cpp` is also the implementation called by `compute_adaptive_threshold` and `threshold` in the step-by-step path — it is the public C++ primitive for filter application. Using it in `lwt_cpp` keeps the step-by-step path consistent and testable. The cost is one additional function call per step and the use of `Rcpp::NumericVector` instead of `std::vector<double>`, which carries SEXP reference-count updates. This overhead is acceptable for the step-by-step path (which is not the performance-critical path) but would be unacceptable inside `denoise_offline_cpp` or `WaveletEngine`.

---

## 4. The Polyphase Decomposition

At each decomposition level, the current approximation signal of length $n$ is split into two subsequences by index parity:

$$\text{even}[i] = x[2i], \quad i = 0, \ldots, \lceil n/2 \rceil - 1$$
$$\text{odd}[i] = x[2i+1], \quad i = 0, \ldots, \lfloor n/2 \rfloor - 1$$

The lengths are $n_e = \lceil n/2 \rceil$ and $n_o = \lfloor n/2 \rfloor$. When $n$ is odd, $n_e = n_o + 1$; when $n$ is even, $n_e = n_o$. The `window_size` parameter in causal mode is forced to odd precisely to guarantee $n_e > n_o$ at every level — this ensures the finest-level even subband always has at least one element, which is required for the predict step to have a valid reference point near the right boundary of the window.

The split is performed by strided indexing — no memory reordering is needed — and the two vectors are passed through the predict and update steps in sequence.

---

## 5. Predict and Update Steps

The predict step estimates each odd sample from its even neighbours and subtracts the estimate to form the detail coefficient:

$$d[i] = \text{odd}[i] - \sum_{k} c_k \cdot \text{even}[i + \text{start\_idx} + k]$$

The update step corrects the even (approximation) subband using the detail coefficients, so that the updated approximation retains the same mean (and, for wavelets with more vanishing moments, higher-order moments) as the original signal:

$$s[i] = \text{even}[i] + \sum_{k} c_k \cdot d[i + \text{start\_idx} + k]$$

Multiple predict/update pairs may be chained within a single decomposition level, as in DB2 (three steps) and CDF 9/7 (four steps). Each step refines the detail or approximation subband incrementally. The order matters: the steps must be applied in the exact order specified in `scheme$steps`, and inverted in exact reverse order during reconstruction.

### 5.1 Normalization

After all predict and update steps at a level are complete, the approximation and detail subbands are scaled by `norm[0]` and `norm[1]` respectively:

```cpp
for (int i = 0; i < n_even; i++) even[i] *= norm_approx;
for (int i = 0; i < n_odd;  i++) odd[i]  *= norm_detail;
```

Normalization does not affect the filter coefficients and plays no role in the `degree` inference or the irregular-grid path. Its purpose is to enforce energy normalization (orthonormal case) or a chosen amplitude convention across levels.

---

## 6. Perfect Reconstruction

The lifting scheme guarantees perfect reconstruction by construction. Each predict step subtracts a linear combination of even samples from odd samples; the inverse adds back the same linear combination. Each update step adds a linear combination of odd (detail) samples to even samples; the inverse subtracts them back. Since both operations are invertible and their composition is the identity, the full multi-level transform is invertible regardless of the choice of coefficients.

Formally, the predict and update steps can be written as lower/upper triangular matrix operations on the polyphase representation, and the inverse is obtained by transposing the sign of each off-diagonal block. This structure is preserved as long as the same coefficients and the same boundary handling are used in the forward and inverse transforms.

In practice, perfect reconstruction is validated by the `diagnose_wavelet` function, which checks that $\| x - \text{ilwt}(\text{lwt}(x)) \|_\infty < 10^{-9}$ for a set of reference signals.

---

## 7. The Irregular Path: Lagrange Interpolation in the Predict Step

For a regular grid, the predict step uses fixed coefficients — the filter assumes uniform spacing between even samples. For an irregular grid, the optimal prediction of $\text{odd}[i]$ at physical position $t_\text{odd}[i]$ depends on the actual positions of the neighbouring even samples $t_\text{even}[j]$.

The irregular path replaces the fixed linear combination with Lagrange polynomial interpolation evaluated at $t_\text{odd}[i]$:

$$d[i] = \text{odd}[i] - p(t_\text{odd}[i])$$

where $p$ is the unique polynomial of degree $k - 1$ passing through the $k$ points $\{(t_\text{even}[j],\, \text{even}[j])\}$ for $j$ in the filter window. The `interp_predict` function in `utils.h` implements this for $k = 1$ (nearest-neighbour), $k = 2$ (linear), and $k \geq 3$ (Lagrange via the standard product formula).

This design embeds the position-aware prediction inside the predict step itself — no pre-interpolation to a regular grid is required. The consequence is that the irregular path produces detail coefficients with the same vanishing-moment properties as the regular path, provided the signal is smooth relative to the grid spacing.

The irregular path is activated per predict step, conditioned on `step.degree >= 0`. Update steps always use fixed coefficients, because the update step corrects the energy balance of the approximation subband — a correction that does not depend on sample positions in the same way.

### 7.1 Out-of-Bounds Position Extrapolation

When the filter window for a predict step extends beyond the signal boundary, `get_val_safe` handles the signal values (using whichever boundary mode is selected), while `get_t_extrap` handles the time positions independently:

```cpp
double get_t_extrap(const std::vector<double>& t, int idx, int n) {
    if (idx >= 0 && idx < n) return t[idx];
    if (idx < 0) {
        double dt = t[1] - t[0];
        return t[0] + (double)idx * dt;
    }
    double dt = t[n-1] - t[n-2];
    return t[n-1] + (double)(idx - (n-1)) * dt;
}
```

Position extrapolation is always linear, regardless of the boundary mode for values. This separation is necessary because applying symmetric or periodic extension to positions would produce duplicate or reversed time coordinates — zero denominators in the Lagrange formula. Linear extrapolation from the nearest boundary spacing is the only extension that preserves the strict ordering required by the interpolation.

### 7.2 Per-level Position Tracking

At each decomposition level, the time-position vector is split in parallel with the signal:

$$t_\text{even}[i] = t[2i], \quad t_\text{odd}[i] = t[2i+1]$$

After the predict and update steps, $t_\text{even}$ becomes the position vector for the next level. These per-level position vectors must be stored during the forward transform and read back during the inverse — a requirement that motivates `t_levels[0..levels]` in `denoise_offline_cpp` and `work_t_approx`/`work_t_detail` in `WaveletEngine`.

The inverse transform reads the stored $t_\text{even}$ and $t_\text{odd}$ to reconstruct the same Lagrange prediction that was subtracted in the forward direction. The subtracted quantity and the added quantity are identical, so perfect reconstruction holds for the irregular path exactly as for the regular path.
