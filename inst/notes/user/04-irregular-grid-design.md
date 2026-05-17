# rLifting — Irregular Grids

This document covers when to use irregular-grid mode, what it costs, how it interacts with each operation mode and each boundary extension, and how to choose between regular and irregular processing.
---

## 1. What Makes a Grid Irregular

A grid is irregular when the physical time positions of samples are non-uniform — the spacing between consecutive samples varies. Examples: ECG with variable heart rate, accelerometer data with dropped packets, seismic sensors at uneven distances, financial tick data.

In the regular path, the predict step uses fixed coefficients that assume uniform spacing. On an irregular grid, those fixed coefficients introduce a residual proportional to the spacing variation: detail coefficients of a smooth signal are no longer near zero, which inflates the adaptive threshold and degrades denoising quality.

The irregular path corrects this by replacing fixed coefficients with **Lagrange interpolation** computed from the physical positions of the neighbours. The transform is then position-aware: the predict step asks "what value does the signal most likely have at position $t_\text{odd}$, given its values at neighbouring even positions $t_\text{even}$?" rather than assuming evenly-spaced samples.

---

## 2. Regular vs Irregular — Trade-offs

| Aspect | Regular | Irregular |
|:-------|:--------|:---------|
| Input | Signal only | Signal + `t` vector (sorted, same length) |
| Predict step | Fixed coefficients | Lagrange interpolation per step |
| Wavelets supported | All built-in | Interpolating only (`haar`, `cdf53`, `dd4`, `lazy`) |
| Speed (offline, n=1024) | baseline | ~1.7× slower |
| Threshold quality on smooth signals | Good for uniform spacing | Better when spacing varies significantly |
| Perfect reconstruction | Always | Always (same interpolation used in inverse) |

**Use regular when:**
- Spacing variation is < 5% — the residual from fixed coefficients is negligible.
- Using `db2` or `cdf97` — these wavelets ignore `t` regardless (non-interpolating; warning emitted).
- Throughput is the priority and spacing variation is acceptable.

**Use irregular when:**
- Spacing variation exceeds ~10% and the signal is smooth relative to the scale of variation.
- Physical positions carry meaning for the prediction (e.g. spatial interpolation, non-uniform sensor arrays).
- Using `cdf53` or `dd4`, which have well-defined positional interpolation.

---

## 3. Interaction with Operation Modes

### 3.1 Step-by-step (Manual)

Pass `t` to `lwt`. The object stores `t` alongside the coefficients; `ilwt` reads it back automatically. Both forward and inverse transforms apply Lagrange interpolation when `degree >= 0`.

The per-level position vectors are recomputed inside `lwt_cpp` and `ilwt_cpp` — no user action required beyond providing `t` once.

### 3.2 Offline

`denoise_signal_offline` passes `t` to `denoise_offline_cpp`, which maintains `t_levels[0..levels]`: per-level even-position vectors built by splitting `t` at each decomposition level, and read back during reconstruction. The overhead relative to the regular path is ~1.7× at n = 1 024 (Lagrange interpolation in each predict step at every level).

### 3.3 Causal Batch

`denoise_signal_causal` receives the full `t` vector and slices it sample by sample inside `run_causal_batch_cpp`. Each sample's position is pushed into `ring_buffer_t` in parallel with the signal value. On linearisation, the position buffer is split into even/odd arrays (`work_t_approx`, `work_t_detail`) before the forward pass and reused during the inverse.

The `WaveletEngine` must be constructed with `irregular = true` (done automatically when `t` is supplied). The overhead per window is proportional to the number of Lagrange evaluations: `n_odd × levels` per forward pass and the same for the inverse.

### 3.4 Stream

The closure signature becomes `processor(new_sample, t_val)`. The user passes the physical recording time of each sample; the engine stores it in `ring_buffer_t[head]` alongside the sample. Everything downstream is identical to Causal Batch — the difference is lifecycle (persistent engine) rather than mechanics.

If `t_val` is omitted, `step_iter` (a sequential integer) is used as a fallback. This is only correct for uniformly-sampled data where absolute positions are unknown; for a genuine irregular grid, always supply `t_val`.

---

## 4. Interaction with Boundary Modes

When the filter window extends beyond the signal (or window) boundary, two things need extending: the **signal values** and the **time positions** of the neighbours.

- **Signal values**: handled by `get_val_safe` using whichever boundary mode is selected (symmetric, periodic, zero, local_linear, one_sided).
- **Time positions**: always extrapolated linearly by `get_t_extrap`, regardless of the boundary mode for values. The spacing used is that of the nearest real pair at the boundary.

This separation is intentional: `get_val_safe` determines what the signal *value* is at a virtual position; `get_t_extrap` determines where that virtual position *is in time*. Using symmetric or periodic extrapolation for positions would produce duplicate or reversed time values, breaking the Lagrange denominator.

| Boundary mode | Value extension | Position extension | Compatible with irregular? |
|:--------------|:----------------|:-------------------|:--------------------------|
| `symmetric` | Mirror reflection | Linear (`get_t_extrap`) | Yes |
| `periodic` | Wrap-around | Linear (`get_t_extrap`) | Yes |
| `zero` | 0 | Linear (`get_t_extrap`) | Yes |
| `local_linear` | OLS extrapolation | Linear (`get_t_extrap`) | Yes |
| `one_sided` | Filter renormalisation | N/A — Lagrange not called | **No** — `t` is ignored; warning raised |

`one_sided` is the only mode incompatible with irregular grids. When `ext_mode == 5`, the C++ engine calls `onesided_conv` before reaching the irregular path, so Lagrange interpolation is never applied even if `t` is supplied. See `03-boundary-modes.md`, section 2.5.

---

## 5. Choosing a Wavelet for Irregular Grids

Only wavelets whose predict step implements an interpolation (`sum(predict_coeffs) ≈ 1`) benefit from positional Lagrange correction. The `degree` field is inferred automatically.

| Wavelet | Interpolating? | Irregular benefit | Notes |
|:--------|:--------------|:------------------|:------|
| `haar` | Yes (`degree = 0`) | Nearest-neighbour correction | Minimal; useful for step-like signals |
| `cdf53` | Yes (`degree = 1`) | Linear interpolation | Best balance for smooth signals on irregular grids |
| `dd4` | Yes (`degree = 3`) | Cubic Lagrange | Best for smooth signals; 2 vanishing moments |
| `lazy` | Yes (`degree = 0`) | Nearest-neighbour | No update step; not recommended for denoising |
| `db2` | No (`degree = -1`) | None — fixed coefficients | Orthogonal; `t` ignored with warning |
| `cdf97` | No (`degree = -1`) | None — fixed coefficients | Biorthogonal; `t` ignored with warning |

For irregular grids, `cdf53` is the default recommendation: it adapts via linear interpolation, has one vanishing moment, and imposes minimal computational overhead. `dd4` is preferable when the signal is known to be smooth and the grid is highly irregular.

---

## 6. Summary

The irregular path is a drop-in extension: supply `t`, and the engine adapts predict steps wherever the wavelet supports it. The cost is ~1.7× in offline mode; causal/stream modes carry similar per-window overhead. The only hard constraint is wavelet choice (`db2`/`cdf97` are excluded) and boundary mode (`one_sided` disables Lagrange silently — use `symmetric` or `local_linear` instead).

→ *See: `implementation/` — implementation documents (to be expanded)*
