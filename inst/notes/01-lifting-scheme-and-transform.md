# Lifting Scheme and Transform — Technical Reference

Companion to `vignette("01-introduction")`. The vignette gives the tour; this note covers the math, data structures, and code paths.

Cross-references: `00-design-overview.md` (architectural index), `02-adaptive-thresholding.md` (threshold pipeline), `03-zero-allocation-engine.md` (engine/ring buffer), `04-boundary-and-threshold.md` (boundary modes, `get_val_safe`, MAD).

---

## 1. Notation and signal model

Throughout this note:

- $x = (x_0, \ldots, x_{n-1}) \in \mathbb{R}^n$ — input signal, discrete index $i$.
- For irregular grids, $t = (t_0, \ldots, t_{n-1})$ with $t_i < t_{i+1}$ — sample positions; $t_i$ is read from R as a real-valued vector parallel to $x$.
- A single decomposition level produces detail $d \in \mathbb{R}^{n_o}$ and approximation $s \in \mathbb{R}^{n_e}$ with $n_e = \lceil n/2 \rceil$, $n_o = \lfloor n/2 \rfloor$.
- Multi-level decomposition writes $d^{(j)}$ for the detail at level $j$ and $a^{(J)}$ for the coarsest approximation after $J$ levels.

The lifting transform is built from three primitives applied in sequence at every level:

1. **Lazy split** (polyphase decomposition by parity).
2. **Predict** (P): one or more steps subtracting an estimate of odd from odd.
3. **Update** (U): one or more steps adding a correction from detail back into even.

Then normalization scales each subband.

---

## 2. The `lifting_scheme` S3 object

The central data structure. Constructed by `lifting_scheme()` in `R/lifting_scheme.R`. Every public function accepts one; its `steps` field is consumed verbatim by C++.

### 2.1 Fields

```r
structure(
  list(wavelet = <chr>, steps = <list>, normalization = <num[2]>),
  class = "lifting_scheme"
)
```

| Field | Type | Used by C++? | Purpose |
|:------|:-----|:-------------|:--------|
| `wavelet` | character(1) | no | Identifier; printed by `print.lifting_scheme`. Does not affect computation. |
| `steps` | list of named lists | yes (verbatim) | Ordered P/U operations. See §2.2. |
| `normalization` | numeric(2) | yes | `c(norm_approx, norm_detail)` applied at the end of each level. |

### 2.2 Step entries

Each element of `steps` is a list with four named fields:

| Field | Type | C++ field | Meaning |
|:------|:-----|:----------|:--------|
| `type` | `"predict"` or `"update"` | `std::string` | Selects P or U branch. |
| `coeffs` | numeric vector, length $k$ | `std::vector<double>` | Filter taps. |
| `start_idx` | integer | `int` | Offset of the first tap relative to the target index. The window reads neighbours at offsets `start_idx + 0..(k-1)`. |
| `degree` | integer | `int` | $\geq 0$: interpolating predict, polynomial degree (enables irregular path). $-1$: fixed coefficients (always for updates, and for orthogonal-wavelet predicts). |

The struct mirroring this in C++ is in `inst/include/rLifting/utils.h`:

```cpp
struct LiftingStep {
    std::string type;
    std::vector<double> coeffs;
    int start_idx;
    int degree = -1;
};
```

### 2.3 Degree inference

Set automatically by `lifting_scheme()` and `lift_step()`:

```r
degree = if (type == "predict" && abs(sum(coeffs) - 1) < 1e-10)
    as.integer(length(coeffs) - 1L)
else
    -1L
```

The test $\sum_k c_k = 1$ is the partition-of-unity condition for interpolating filters: a polynomial of degree $\leq k-1$ evaluated by $\sum c_k \cdot p(t_k)$ reproduces $p$ exactly at the target if and only if the coefficients sum to one (constant case) plus the higher-order vanishing-moment constraints encoded in the specific $c_k$. The constructor only checks the constant case as a necessary condition; the user is trusted to provide a genuine interpolating filter when `sum == 1`. The degree itself ($k-1$) is the maximum polynomial degree the filter can reproduce.

The `lift_step()` helper accepts an explicit `degree` argument that overrides the inference, and accepts `position = "center"/"left"/"right"` to derive `start_idx` automatically.

### 2.4 Normalization

`normalization = c(K, 1/K)` for orthogonal/biorthogonal wavelets. The reciprocal relation is enforced by convention but not by code; `lifting_scheme()` accepts any pair. After all P/U steps at a level, the approximation is multiplied by `K` and the detail by `1/K`; the inverse divides instead.

---

## 3. Built-in wavelets

All defined in `.get_wavelet_config()` in `R/lifting_scheme.R`.

| Name | Steps | Normalization | Inferred degrees | Reference |
|:-----|:------|:--------------|:-----------------|:----------|
| `lazy` | none | $(1, 1)$ | — | Identity split |
| `haar` | P: $[1]$, $s_0{=}0$;  U: $[0.5]$, $s_0{=}0$ | $(\sqrt 2, 1/\sqrt 2)$ | P: $0$ (nearest-neighbour) | classical |
| `cdf53` | P: $[\tfrac12, \tfrac12]$, $s_0{=}0$;  U: $[\tfrac14, \tfrac14]$, $s_0{=}{-1}$ | $(\sqrt 2, 1/\sqrt 2)$ | P: $1$ (linear) | Cohen–Daubechies–Feauveau (1992) |
| `dd4` | P: $[-\tfrac1{16}, \tfrac9{16}, \tfrac9{16}, -\tfrac1{16}]$, $s_0{=}{-1}$;  U: $\tfrac12 \cdot p$, $s_0{=}{-1}$ | $(\sqrt 2, 1/\sqrt 2)$ | P: $3$ (cubic Lagrange) | Deslauriers–Dubuc (1989) |
| `db2` | P: $[\sqrt 3]$;  U: $[\tfrac{\sqrt 3}{4}, \tfrac{\sqrt 3 - 2}{4}]$;  P: $[-1]$, $s_0{=}{-1}$ | $(\tfrac{\sqrt 3 + 1}{\sqrt 2}, \tfrac{\sqrt 3 - 1}{\sqrt 2})$ | all $-1$ (orthogonal, non-interpolating) | Daubechies–Sweldens (1998) |
| `cdf97` | P/U/P/U with $\alpha, \beta, \gamma, \delta$ | $(\zeta, 1/\zeta)$, $\zeta \approx 1.1496$ | all $-1$ | Cohen–Daubechies–Feauveau (1992); Daubechies–Sweldens (1998) |

`cdf53` and `dd4` are interpolating: their predict coefficients sum to $1$, so `degree` is set automatically and the irregular Lagrange path activates when `t` is supplied. `db2` and `cdf97` are orthogonal/biorthogonal factorizations whose predict coefficients do not sum to $1$; passing `t` triggers a warning from `.check_irregular_scheme()` and the C++ path falls back to fixed coefficients (still bit-exact reconstructible).

---

## 4. R → C++ boundary

### 4.1 Marshalling

`lwt()` passes four things to `lwt_cpp`:

```r
lwt_cpp(
    as.numeric(signal),
    scheme$steps, # Rcpp::List of List
    as.numeric(scheme$normalization), # NumericVector length 2
    as.integer(levels),
    as.integer(ext_int), # 1..5
    t_cpp, # numeric(0) for regular grid
    as.integer(ll_k) # local-linear neighbourhood
)
```

The `steps` list is consumed by C++ exactly as it appears in R; no transformation. In `lwt_cpp` in `src/lwt.cpp`, each step is read field-by-field:

```cpp
List step = steps[k];
std::string type = as<std::string>(step["type"]);
NumericVector coeffs = step["coeffs"];
int start_idx = step["start_idx"];
int degree = step.containsElementNamed("degree")
               ? (int)step["degree"] : -1;
```

The same parsing pattern appears at four sites — `lwt_cpp`, `ilwt_cpp`, `denoise_offline_cpp` (`src/offline.cpp`), and `WaveletEngine` constructor (`inst/include/rLifting/WaveletEngine.h`). Any change to the step layout must touch all four. This is the cost of avoiding a parsed `std::vector<LiftingStep>` crossing the R/C++ boundary; parsing once per call is cheap.

### 4.2 Boundary mode encoding

The R-side string is mapped to an integer by `switch()` in `R/lwt.R`, `R/ilwt.R`, `R/denoising_offline.R`, and `R/realtime_denoising.R`:

| String | Integer | Code path |
|:-------|--------:|:----------|
| `symmetric` | 1 | `get_val_safe` (reflection) |
| `periodic` | 2 | `get_val_safe` |
| `zero` | 3 | `get_val_safe` |
| `local_linear` | 4 | `get_val_safe` (k-point OLS) |
| `one_sided` | 5 | `apply_filter_cpp` early branch; `onesided_conv` |

Modes 1–4 are pure single-index virtualisations and live entirely in `get_val_safe`. Mode 5 renormalises the entire filter window and therefore needs explicit branches in `apply_filter_cpp` and in every inline lifting loop (offline + engine). See `04-boundary-and-threshold.md` §1 for the full table.

### 4.3 Regeneration

Any change to a `// [[Rcpp::export]]` signature requires `Rcpp::compileAttributes()` to regenerate `src/RcppExports.cpp` and `R/RcppExports.R`. The current export signatures are visible in `src/RcppExports.cpp` (`compileAttributes()` is the source of truth). Forgetting to regenerate leaves the R stub calling an old ABI.

---

## 5. Polyphase decomposition (lazy split)

At each level, the current approximation $x$ of length $n$ is split by index parity:

$$
\text{even}[i] = x[2i], \quad i = 0, \ldots, n_e - 1, \qquad n_e = \lceil n/2 \rceil
$$
$$
\text{odd}[i]  = x[2i+1], \quad i = 0, \ldots, n_o - 1, \qquad n_o = \lfloor n/2 \rfloor
$$

For even $n$, $n_e = n_o$; for odd $n$, $n_e = n_o + 1$. The split is strided indexing — no buffer reordering — done in `lwt_cpp` in `src/lwt.cpp`. Position vectors split in parallel when `t` is supplied (same function).

Why odd $n$ matters for streaming: `window_size` is forced odd so that $n_e > n_o$ at every level, guaranteeing the finest-level even subband has at least one element near the right edge of the window (needed as a reference point for the right-boundary predict).

---

## 6. Predict and update steps

### 6.1 Forward equations

A predict step with coefficients $c_0, \ldots, c_{k-1}$ and `start_idx` $s$ computes, for each odd sample,

$$
d[i] \;\leftarrow\; \text{odd}[i] \;-\; \sum_{j=0}^{k-1} c_j \cdot \widetilde{\text{even}}[i + s + j],
$$

where $\widetilde{\text{even}}$ denotes the boundary-extended view (`get_val_safe` for modes 1–4, renormalised filter for mode 5). This is implemented by `apply_filter_cpp` in `src/utils.cpp` when `degree < 0` and by `predict_irregular` in `src/lwt.cpp` when `degree >= 0` and `t` is supplied.

An update step is symmetric:

$$
s[i] \;\leftarrow\; \text{even}[i] \;+\; \sum_{j=0}^{k-1} c_j \cdot \widetilde{\text{odd}}[i + s + j].
$$

Update steps always go through `apply_filter_cpp`; the irregular path is not used because the update is a correction on the approximation subband (energy balance / vanishing-moment shift on the lowpass), not a position-aware interpolation.

After all steps run, normalization is applied (final block of `lwt_cpp` in `src/lwt.cpp`):

```cpp
for (int m = 0; m < n_even; m++) even[m] *= norm[0];
for (int m = 0; m < n_odd;  m++) odd[m]  *= norm[1];
```

The detail $d^{(j)}$ is stored as `coeffs_out["dj"]`; the new approximation becomes `even` for the next level.

### 6.2 Worked example: CDF 5/3

`cdf53` has two steps:

$$
d[i] = \text{odd}[i] - \tfrac{1}{2}\bigl(\text{even}[i] + \text{even}[i+1]\bigr)
$$
$$
s[i] = \text{even}[i] + \tfrac{1}{4}\bigl(d[i-1] + d[i]\bigr)
$$

The predict has `coeffs = c(0.5, 0.5)`, `start_idx = 0`, `degree = 1` (linear). The update has `coeffs = c(0.25, 0.25)`, `start_idx = -1`, `degree = -1`. On a regular grid the predict reduces to linear interpolation between the two flanking even samples; on an irregular grid the same step takes the Lagrange branch and interpolates using actual positions $(t_{\text{even}}[i], t_{\text{even}}[i+1])$ evaluated at $t_{\text{odd}}[i]$.

### 6.3 Worked example: DB2

`db2` has three steps:

$$
d[i] \leftarrow \text{odd}[i] - \sqrt 3 \cdot \text{even}[i]
$$
$$
s[i] \leftarrow \text{even}[i] + \tfrac{\sqrt 3}{4} d[i] + \tfrac{\sqrt 3 - 2}{4} d[i+1]
$$
$$
d[i] \leftarrow d[i] + s[i-1]
$$

(this is equivalent to the textbook form `d[i] -= s[i-1]` after factoring out the sign; the implementation chose `coeffs = c(-1)` with `start_idx = -1`, so `pred[i] = (-1) · s[i-1]` and `d[i] -= pred[i]` yields `+s[i-1]`.)

The three steps are followed by normalization by $((\sqrt 3 + 1)/\sqrt 2,\, (\sqrt 3 - 1)/\sqrt 2)$. The first predict has $\sum c_k = \sqrt 3 \neq 1$, so `degree = -1` and the irregular path is disabled — orthogonality conflicts with interpolation, which is by design.

### 6.4 Inverse

`ilwt_cpp` in `src/ilwt.cpp` runs the steps in **exact reverse order** with signs flipped: each update becomes a subtraction, each predict an addition. Normalization is undone first (`ilwt_cpp` in `src/ilwt.cpp`):

```cpp
for (int i = 0; i < n_even; i++) even[i] /= norm[0];
for (int i = 0; i < n_odd;  i++) odd[i]  /= norm[1];
```

then the step loop iterates `k = n_steps - 1; k >= 0; k--` (same `ilwt_cpp` function). After the loop, even/odd are interleaved back into the parent approximation:

```cpp
for (int i = 0; i < n_even; i++) merged[2*i]     = even[i];
for (int i = 0; i < n_odd;  i++) merged[2*i + 1] = odd[i];
```

Final length is trimmed to `original_len` because the polyphase recursion can pad up to one sample per level on odd lengths.

### 6.5 Perfect reconstruction

The lifting scheme is exactly invertible by construction: every predict step is $d \mathrel{-}= P(\text{even})$ and the inverse is $d \mathrel{+}= P(\text{even})$ with the same $P$; analogously for updates. Composition of invertible elementary matrices on the polyphase representation is invertible regardless of coefficients, so `db2`, `cdf97`, or any user-supplied scheme reconstructs perfectly **as long as the inverse uses the same boundary mode and the same `t`**. `diagnose_wavelet()` in `R/diagnostics.R` validates this numerically.

---

## 7. Multiple implementations of the transform loop

The same algorithm is implemented at three sites because the data structures differ:

| Site | File | Operates on | Calls `apply_filter_cpp`? |
|:-----|:-----|:------------|:--------------------------|
| `lwt_cpp` / `ilwt_cpp` | `src/lwt.cpp`, `src/ilwt.cpp` | `Rcpp::NumericVector` | yes (also for U; P only when `degree < 0`) |
| `denoise_offline_cpp` | `src/offline.cpp` | `std::vector<double>` per level | no — inline loops |
| `WaveletEngine::push_and_process` | `inst/include/rLifting/WaveletEngine.h` | pre-allocated workspace `std::vector<double>` | no — inline loops |

`lwt_cpp` delegates to `apply_filter_cpp` because it is the public C++ primitive also used by adaptive thresholding and the standalone `apply_filter()` R wrapper; keeping the step-by-step path consistent and testable is worth one extra function call per step. `denoise_offline_cpp` and `WaveletEngine` inline the loops to avoid the SEXP refcount overhead of `Rcpp::NumericVector` on hot paths. See `03-zero-allocation-engine.md` for the engine's allocation discipline.

The trade-off: every change to predict/update semantics must be made in all three places (plus the four marshalling sites of §4.1). The mandatory boundary-mode mirroring is enumerated in `04-boundary-and-threshold.md`.

---

## 8. Irregular grid: Lagrange interpolation in the predict step

### 8.1 Motivation

On a regular grid, the cdf53 predict step
$\widehat{\text{odd}}[i] = \tfrac12(\text{even}[i] + \text{even}[i+1])$
is exact linear interpolation because $\text{odd}[i]$ sits exactly halfway between its even neighbours. On an irregular grid this midpoint assumption fails: the odd sample at position $t_{\text{odd}}[i]$ may be anywhere between $t_{\text{even}}[i]$ and $t_{\text{even}}[i+1]$. Using the regular coefficients would inject systematic bias into the detail $d[i]$, raising its variance and degrading thresholding.

### 8.2 Activation

Per predict step, the irregular path is used iff:

1. `t.size() == signal.size()` at the top of `lwt_cpp` (sets `irregular = true`), AND
2. `step.degree >= 0` for that specific step.

If both hold, control goes to `predict_irregular` in `src/lwt.cpp`; otherwise to `apply_filter_cpp`. Update steps **never** take the irregular path, regardless of `degree`.

### 8.3 The Lagrange evaluation

`predict_irregular` builds two parallel windows of length $k$ (the filter length):

```cpp
for (int j = 0; j < k; j++) {
    int idx = i + start_idx + j;
    x_nbr[j] = get_val_safe(even_v, idx, n_even, ext_mode, ll_k);
    t_nbr[j] = get_t_extrap(t_even, idx, n_even);
}
pred[i] = interp_predict(x_nbr, t_nbr, k, t_odd[i]);
```

`interp_predict` in `inst/include/rLifting/utils.h` computes

$$
p(t_*) = \sum_{i=0}^{k-1} x_i \prod_{\substack{j=0 \\ j \neq i}}^{k-1} \frac{t_* - t_j}{t_i - t_j},
$$

the unique polynomial of degree $\leq k-1$ through the $k$ points $(t_i, x_i)$, evaluated at $t_* = t_\text{odd}[i]$. Special cases:

- $k = 1$: returns $x_0$ (nearest neighbour, Haar).
- $k = 2$: explicit linear formula with degenerate-spacing guard ($|t_1 - t_0| < 10^{-15}$).
- $k \geq 3$: the general product formula, with the same denominator guard per term.

The `interp_predict` header comment in `inst/include/rLifting/utils.h` and the `LiftingStep` struct comment now describe the three cases consistently (k=1 nearest, k=2 linear, k≥3 general Lagrange). No legacy inconsistency remains; the loop body matches the documented cases.

### 8.4 Position extrapolation at boundaries: `get_t_extrap`

When `idx = i + start_idx + j` falls outside $[0, n_e)$, `get_val_safe` extends the **values** according to the active boundary mode, but the **positions** are extended by unconditional linear extrapolation (`get_t_extrap` in `inst/include/rLifting/utils.h`):

$$
t[idx] = \begin{cases}
t[0] + idx \cdot (t[1] - t[0]) & \text{if } idx < 0 \\
t[n-1] + (idx - (n-1)) \cdot (t[n-1] - t[n-2]) & \text{if } idx \geq n
\end{cases}
$$

Why a separate mechanism for positions: applying `symmetric` reflection to time positions would generate duplicates ($t[-1] = t[0]$) and would break the strict ordering $t_0 < t_1 < \cdots$ that the Lagrange denominators require — a duplicate $t_j$ produces a zero denominator and is guarded out, but a non-monotone $t$ sequence produces wrong sign-cancellations in the product and an arbitrary interpolant. Linear extrapolation preserves monotonicity and is the unique choice that does so without injecting a model of the underlying grid.

The value/position split is documented in detail in `04-boundary-and-threshold.md` §5.

### 8.5 Per-level position tracking

The level-$j$ approximation is the even sub-sequence of level $j-1$, so positions descend the same way:

$$
t^{(j)}_i = t^{(j-1)}_{2i}.
$$

`lwt_cpp` in `src/lwt.cpp` advances `current_t = t_even` after each level. `ilwt_cpp` in `src/ilwt.cpp` cannot reconstruct this on the fly — by the time the inverse reaches level $j$, it only has `current_app` (the merged signal), not the split — so it pre-builds `t_levels[0..J-1]` by repeated even-subsampling of the full $t$ at the top. The inverse then reads `t_levels[j-1]` at each iteration, splits it into `t_even`/`t_odd`, and reproduces the same Lagrange predictor that was subtracted in the forward direction.

Perfect reconstruction holds for the irregular path because the subtracted and added quantities are computed by the same `interp_predict` call on the same window. The only requirement is that `lwt_obj$t` survives the round trip — `ilwt()` in `R/ilwt.R` takes it from the `lwt` object's `$t` slot, which is set by `lwt()` in `R/lwt.R`. Denoising functions (`denoise_signal_offline`, `WaveletEngine`) maintain per-level positions in their own workspaces (`t_levels[]` in `offline.cpp`; `work_t_approx`/`work_t_detail` in the engine).

---

## 9. Output structure

`lwt()` returns an S3 `lwt` object:

```r
structure(
  list(
    coeffs = list(d1 = ..., d2 = ..., ..., aJ = ...),
    scheme = <lifting_scheme>,
    levels = J,
    original_len = n,
    extension = <chr>,
    ll_k = <int>,
    t = <num | NULL>
  ),
  class = "lwt"
)
```

`coeffs` follows a fixed naming convention: `d1` is the finest-level detail (highest frequency), `dJ` the coarsest detail, `aJ` the coarsest approximation. The order is set in `lwt_cpp` by the construction order of `coeffs_out`. Length per level follows the recursion $n^{(0)} = n$, $n^{(j)} = \lceil n^{(j-1)}/2 \rceil$, $|d^{(j)}| = \lfloor n^{(j-1)}/2 \rfloor$, $|a^{(J)}| = n^{(J)}$. For $n$ divisible by $2^J$ this simplifies to $|d^{(j)}| = n/2^j$ and $|a^{(J)}| = n/2^J$; for odd $n$ the recursion must be applied level by level (e.g.\ $n=7$, $J=2$: $n^{(1)}=4$, $|d^{(1)}|=3$, $n^{(2)}=2$, $|d^{(2)}|=2$, $|a^{(2)}|=2$).

`ilwt()` requires the full `lwt` object — it cannot be invoked with an isolated coefficient list, because `extension`, `ll_k`, and `t` are needed to reproduce the exact forward step. This is enforced by the `inherits(lwt_obj, "lwt")` check at the top of `ilwt()` in `R/ilwt.R`.

---

## 10. Where to look next

- For the next stage of the pipeline (thresholding the `d_j`): `02-adaptive-thresholding.md`.
- For why the engine inlines the loop and how the ring buffer keeps allocations zero: `03-zero-allocation-engine.md`.
- For the full table of boundary modes and the four mandatory code paths: `04-boundary-and-threshold.md`.
- User-facing usage: `vignette("01-introduction")`, `vignette("05-irregular-grids")`, `vignette("06-extensions")`.
