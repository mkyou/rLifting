# rLifting — Zero-Allocation Engine

Implementation reference for the C++ computational core. Companion to `vignette("03-causal-stream")` (user-facing tour of the two sliding-window modes). This note documents the `WaveletEngine` class layout, the ring buffer mechanics, the XPtr finalizer pattern, the per-sample hot path, and the asymmetry between the offline path and the engine-based paths.

Cross-references:
- `01-lifting-scheme-and-transform.md` — `LiftingStep` struct, polyphase decomposition, predict/update math.
- `02-adaptive-thresholding.md` — MAD/SURE rules, the `update_freq` protocol consumed in the hot path.
- `04-boundary-and-threshold.md` — the four mandatory C++ boundary code paths (relevant whenever a new boundary mode is added).

---

## 1. What "Zero-Allocation" Means Here

The phrase refers to the **per-sample hot path** (`process_sample_cpp` → `WaveletEngine::push_and_process`). After warm-up there are:

- No R-level allocations: `process_sample_cpp` accepts scalar `double` arguments and returns a `double`. No `NumericVector`, no `SEXP` wrapping per sample.
- No C++ heap allocations: every workspace consumed by `push_and_process` is sized to its final extent in the constructor; all assignments are stride-1 writes into pre-sized buffers. The `resize` calls inside the forward loop of `push_and_process` in `WaveletEngine.h` are no-ops in steady state — see §4.

What is *not* zero-allocation:
- **Construction.** `WaveletEngine::WaveletEngine` allocates the ring buffer plus $2(L+1)$ workspace vectors (plus $2L+1$ position vectors if `irregular`). One-time cost.
- **Offline path.** `denoise_offline_cpp` is a one-shot routine and intentionally allocates locally (see §6).
- **R closure invocation.** `new_wavelet_stream`'s returned closure pays the R function-call overhead per sample (≈ 50–110% above `denoise_signal_causal`'s per-sample cost, measured in `vignette("03-causal-stream")` §5). The C++ side of the call is allocation-free.

---

## 2. WaveletEngine — Class Layout

Defined at `inst/include/rLifting/WaveletEngine.h`. Header-only; downstream packages can link against it without rebuilding rLifting.

```cpp
class WaveletEngine {
public:
  // Configuration (fixed at construction)
  std::vector<LiftingStep> steps;     // copy of the R scheme's $steps
  double norm_approx, norm_detail;    // from scheme$normalization
  int    levels;                      // L
  int    window_size;                 // W (odd; enforced by the R wrapper)
  int    ext_mode;                    // 1..5 (see 04-boundary-and-threshold.md)
  int    ll_k;                        // local-linear neighbourhood
  bool   irregular;                   // toggles the t-position workspaces
  double scad_a;                      // SCAD shape parameter (default 3.7)
  bool   lambdas_initialized;         // gates the single update when update_freq <= 0

  // Ring buffer state
  std::vector<double> ring_buffer;    // W doubles
  std::vector<double> ring_buffer_t;  // W doubles, irregular only
  int head;                           // next write slot in [0, W)
  int count;                          // samples seen, capped at W (warm-up gate)

  // Forward/inverse workspaces (pre-allocated)
  std::vector<double>               work_signal;     // W
  std::vector<double>               work_t;          // W, irregular only
  std::vector<std::vector<double>>  work_approx;     // L+1 jagged levels
  std::vector<std::vector<double>>  work_detail;     // L+1 jagged levels
  std::vector<std::vector<double>>  work_t_approx;   // L+1, irregular only
  std::vector<std::vector<double>>  work_t_detail;   // L,   irregular only

  // Threshold cache
  std::vector<double> current_lambdas; // size L, refreshed every update_freq pushes
};
```

`LiftingStep` is the shared struct from `inst/include/rLifting/utils.h` (`type`, `coeffs`, `start_idx`, `degree`). The same struct is consumed by `denoise_offline_cpp`, `lwt_cpp`, and `ilwt_cpp` — see `01-lifting-scheme-and-transform.md` §2.

### 2.1 Constructor allocations

The `WaveletEngine` constructor in `WaveletEngine.h` takes the lifting steps, normalization, levels, window size, extension mode, irregular flag, local-linear neighbourhood `ll_k`, and the SCAD shape parameter `scad_a_val` (default 3.7). It allocates:

- `ring_buffer.resize(W, 0.0)` — circular sample storage.
- `work_signal.resize(W)` — linearized (chronological) view of the ring.
- `work_approx[j].resize(current_len)` and `work_detail[j].resize(current_len)` for `j = 0..L`, with `current_len` halving each level via `(current_len + 1) / 2`. Total approximately $4W$ doubles, but see §4 on the over-allocation of `work_detail`.
- `current_lambdas.resize(L, 0.0)` — one threshold per detail level.
- If `irregular`: `ring_buffer_t`, `work_t`, `work_t_approx[0..L]`, and `work_t_detail[0..L-1]` mirror the structure for position tracking. Note `work_t_detail` has length `L`, not `L+1` (positions are only needed at the L detail subbands; see the irregular branch of the constructor in `WaveletEngine.h`).

`scad_a` and `lambdas_initialized` are POD scalars stored inline — no additional heap allocations.

For $W = 255$, $L = 4$ the engine holds roughly 8 KB of doubles — small enough to stay resident in L1 on contemporary cores.

---

## 3. Ring Buffer Mechanics

The ring is a single contiguous `std::vector<double>` of length $W$, indexed by a head pointer that wraps modulo $W$. Push is O(1) and triggers no memory traffic outside the buffer itself:

```cpp
// push_and_process in WaveletEngine.h — ring insertion
ring_buffer[head] = new_val;
if (irregular) ring_buffer_t[head] = t_val;
head = (head + 1) % window_size;
if (count < window_size) count++;
```

`count` is a saturating counter: it tracks the warm-up phase. While `count < W` the engine returns the raw sample (warm-up guard in `push_and_process` in `WaveletEngine.h`). The first filtered output is the $W$-th sample, when the buffer first fills. The R wrappers force `W` odd by adding 1 if needed (`denoise_signal_causal` and `new_wavelet_stream` in `realtime_denoising.R`).

### 3.1 Why a ring buffer and not a sliding deque

A `std::deque<double>` with `push_back`/`pop_front` would be O(1) amortized but pays for it with:

- non-contiguous memory (deque uses chunked storage), which defeats the stride-1 access pattern of the LWT;
- pointer chasing across chunk boundaries on every level's downsample;
- additional dereferences in the inner loops.

A `std::vector` with `push_back`/`erase(begin())` would force $W$-element memmoves every sample, an O(W) cost that defeats the linear-in-T total work claim.

The ring achieves O(1) push, contiguous storage, and predictable cache behaviour. The single O(W) cost per sample is the linearization in §4.1, which is a stride-1 copy that runs at memory bandwidth.

### 3.2 Linearization

The LWT consumes a chronologically ordered window. The ring's chronological view starts at `head` (the oldest sample, which is also the next write slot) and wraps. One pass copies it into `work_signal`:

```cpp
// linearization in push_and_process in WaveletEngine.h
for (int i = 0; i < window_size; i++)
    work_signal[i] = ring_buffer[(head + i) % window_size];
work_approx[0] = work_signal;  // initial level-0 approximation
```

The assignment `work_approx[0] = work_signal` (`push_and_process` in `WaveletEngine.h`) is a `std::vector` copy-assignment; both operands are pre-sized to `W`, so it dispatches to a `memcpy` with no allocation. Same for the irregular branch in the same function.

---

## 4. The Per-Sample Hot Path

`push_and_process` in `WaveletEngine.h` is the inner loop. After warm-up:

1. **Insert** the sample into the ring.
2. **Linearize** the ring into `work_signal`; same for `work_t` if irregular.
3. **Forward LWT** — $L$ levels of split + predict + update + scale.
4. **Threshold update** — refresh `current_lambdas` according to the gate below. See `02-adaptive-thresholding.md` §6.2.

   ```cpp
   bool should_update = (update_freq > 0)
                           ? (step_iter % update_freq == 0)
                           : !lambdas_initialized;
   if (should_update) {
     update_thresholds(alpha, beta, threshold_method);
     lambdas_initialized = true;
   }
   ```

   `update_freq <= 0` means "freeze after warm-up" — `lambdas_initialized` gates the single initial update.
5. **Shrinkage** — apply hard/soft/semisoft/scad in-place on `work_detail[0..L-1]`.
6. **Inverse LWT** — reverse the predict/update steps, merge even/odd back into `work_approx[j]` for `j = L-1..0`.
7. **Return** `work_approx[0][window_size - 1]` — the reconstructed value at the most recently inserted position.

All numeric work in steps 3–6 happens inside pre-allocated buffers. The `even.resize(n_even)` / `odd.resize(n_odd)` calls in the forward level loop of `push_and_process` are the only `resize` calls in the hot path:

- `work_approx[j+1]` was sized in the constructor to `ceil(W / 2^{j+1})`, which equals `n_even` for every iteration. The resize is a no-op (size unchanged).
- `work_detail[j]` was sized in the constructor to `ceil(W / 2^j)` (same as `work_approx[j]`), but the hot path needs `n_j / 2` (integer division) = `n_odd`, where `n_j = ceil(W / 2^j)` is the level-j approximation size. The first call downsizes the vector. `std::vector::resize` with a smaller argument never reallocates — it only updates the `size_` field; capacity is retained. From the second call onward the size matches and the call is a no-op.

The over-allocation of `work_detail` in the `WaveletEngine` constructor is the deliberate price for keeping the hot path branch-free of allocation logic. The irregular mirror `work_t_detail`, also sized inside the constructor, is sized to `current_len / 2` from the start, matching the steady-state size.

### 4.1 Output index

The return slot `work_approx[0][window_size - 1]` is the **rightmost** position in the linearized window — i.e., the sample just inserted. After step 6, `work_approx[0]` has been fully reconstructed from the thresholded coefficients, so this is the filtered estimate of the newest sample. The first `W − 1` samples of the linearized window are also filtered, but they correspond to past inputs that were already emitted in previous calls; they are discarded as part of the sliding-window semantics.

### 4.2 Boundary modes in the hot path

The `get_val` helper in `WaveletEngine.h` wraps `get_val_safe` from `utils.h` and handles modes 1–4 inline at every fetch. Mode 5 (`one_sided`) requires renormalisation that cannot be expressed as a single-index virtual value; it has its own branch via `onesided_conv` (invoked in the predict, update, inverse-predict, and inverse-update blocks of `push_and_process` in `WaveletEngine.h`). This is the same four-paths invariant documented in `04-boundary-and-threshold.md` Part A: `WaveletEngine.h` is one of the four files that must mirror every boundary mode.

### 4.3 Irregular path in the hot path

When `irregular = true`, predict steps with `degree >= 0` dispatch to `interp_predict` (Lagrange) inside the forward and inverse predict branches of `push_and_process` in `WaveletEngine.h`. The neighbour positions come from `work_t_approx[j+1]` and the target positions from `work_t_detail[j]`; both are populated from `work_t_approx[j]` at the start of the level (position-split block of `push_and_process`). `get_t_extrap` handles out-of-bounds positions via linear extrapolation (defined in `utils.h`, see `01-lifting-scheme-and-transform.md` §8.4).

---

## 5. The XPtr Finalizer Pattern

`new_wavelet_stream` in `realtime_denoising.R` creates an `Rcpp::XPtr<WaveletEngine>` by calling `create_engine_cpp` in `fast_stream.cpp`:

```cpp
WaveletEngine *engine = new WaveletEngine(steps, norm, levels, window_size,
                                          ext_mode, irregular, ll_k, scad_a);
Rcpp::XPtr<WaveletEngine> ptr(engine, true);  // true = register default finalizer
return ptr;
```

The second argument to the `XPtr` constructor is `register_finalizer`. When it is `true`, Rcpp installs a finalizer that calls `delete engine` when R garbage-collects the external-pointer SEXP. The closure returned by `new_wavelet_stream` captures `engine_ptr` in its environment (closure body of `new_wavelet_stream` in `realtime_denoising.R`); when the closure goes out of scope and is collected, the `XPtr` becomes unreachable, the finalizer fires, the `WaveletEngine` destructor runs, and all `std::vector` members release their heap storage.

Per sample, `process_sample_cpp` in `fast_stream.cpp` reconstitutes the typed XPtr from the SEXP and dispatches:

```cpp
Rcpp::XPtr<WaveletEngine> engine(engine_ptr);
return engine->push_and_process(new_sample, t_val, alpha, beta, method,
                                update_freq, step_iter, threshold_method);
```

One pointer reconstitution, one (non-virtual) method call. No `SEXP` allocation on the return path: Rcpp's `// [[Rcpp::export]]` boilerplate wraps the returned `double` as a length-1 `REALSXP`, which is the only per-sample R-side allocation. This is unavoidable at the R/C++ boundary.

### 5.1 Why not store engine state in R

`WaveletEngine` holds non-trivially-destructible C++ members (`std::vector<LiftingStep>`, jagged vectors) that cannot be marshalled into an R object without serialization. An XPtr is the standard Rcpp idiom for opaque, type-stable C++ state with R-managed lifetime.

---

## 6. The Offline Path Does Not Use WaveletEngine

`denoise_offline_cpp` in `src/offline.cpp` is a separate implementation of LWT → threshold → ILWT, with its own inline loops. The CLAUDE.md "four boundary code paths" rule lists it as a distinct path that must be updated whenever a boundary mode is added.

Why the duplication:

- **Access pattern.** Offline has the full signal of length $n$. The engine processes one sample at a time over a fixed-size window; running it sample-by-sample over an $n$-sample signal would require $n$ linearizations and a per-window MAD instead of a single global MAD, which is the wrong statistic for the non-causal case (see `02-adaptive-thresholding.md` §6.1).
- **Allocation behaviour.** Offline allocates locally per level: `even`, `odd`, `details[j]`, optionally `t_even`/`t_odd`/`t_levels[j]`. These are `std::vector` on the C++ stack frame, freed on return. The cost is paid once per call and is negligible against the $O(n \log n)$ LWT work — pre-allocation would offer no measurable benefit for a one-shot routine. Benchmark evidence in `vignette("01-introduction")` and `02-adaptive-thresholding.md` shows the offline single-pass design is ≈ 4.4× faster than routing `lwt → threshold → ilwt` through R, attributable to eliminating intermediate SEXP boxing rather than to allocation strategy.
- **Inverse reconstruction.** Offline produces a fully-reconstructed signal of length `original_len` (return block of `denoise_offline_cpp` in `offline.cpp`), not a single sample at the rightmost slot. The merging step (inverse-level loop in `denoise_offline_cpp`) allocates a new `merged` vector at each inverse level; this is structural, not a hot-path concern.

The engine's design optimizes for repeated, stateful processing with a fixed memory footprint. The offline design optimizes for one-shot full-signal throughput. The two coexist; consolidating them would require sacrificing one set of constraints.

---

## 7. Causal Batch vs Stream

Both modes share the engine; they differ in how `push_and_process` is driven.

### 7.1 `denoise_signal_causal` → `run_causal_batch_cpp`

`run_causal_batch_cpp` (defined in `fast_stream.cpp`) instantiates a transient stack-allocated `WaveletEngine` and loops over the input signal in C++:

```cpp
WaveletEngine engine(steps, norm, levels, window_size, ext_mode, irreg, ll_k);
for (int i = 0; i < n; i++) {
  double t_val = irreg ? t[i] : 0.0;
  output[i] = engine.push_and_process(signal[i], t_val, alpha, beta, method,
                                      update_freq, i, threshold_method);
}
```

Single R→C++ crossing; the engine lives on the C++ stack and is destroyed when the function returns. No XPtr is needed because there is no R-side handle to manage.

### 7.2 `new_wavelet_stream` → `process_sample_cpp`

The closure pays one R→C++ crossing per sample. The engine persists across calls via XPtr (§5). MSE outputs are bit-equivalent to the batch path when both receive the same sample sequence under matching parameters (the closure's `step_iter` is incremented in R; `run_causal_batch_cpp` increments it in C++; both pass an identical integer to `push_and_process`).

The vignette §5 quantifies the per-sample cost difference: stream adds 50–110% over causal, all of it attributable to the R closure invocation.

---

## 8. Allocation Audit by Path

| Path | Per-call R allocations | Per-call C++ allocations | Per-sample R allocations | Per-sample C++ allocations |
|:-----|:----------------------|:-------------------------|:-------------------------|:---------------------------|
| `denoise_signal_offline` | input copy in, output copy out | per-level `even`/`odd`/`details[j]` | n/a (batch) | n/a |
| `denoise_signal_causal` | input copy in, output buffer alloc | engine constructor (one-shot) | none | none in steady state |
| `new_wavelet_stream` closure | engine constructor (one-shot) | engine constructor (one-shot) | length-1 `REALSXP` for the return value | none in steady state |

"None in steady state" means: after the first call, all `resize` calls inside `push_and_process` are no-ops (§4). The constructor pays the heap cost once.

---

## 9. Cross-Reference Index

When modifying the engine, the following files must stay in sync:

| Concern | File(s) |
|:--------|:--------|
| Boundary modes (single-index) | `inst/include/rLifting/utils.h` (`get_val_safe`) |
| Boundary mode 5 (renormalised) | `WaveletEngine.h`, `src/offline.cpp`, `src/utils.cpp` (`apply_filter_cpp`) |
| `LiftingStep` struct | `LiftingStep` in `inst/include/rLifting/utils.h` (consumed by all four C++ paths) |
| Threshold formulas | `WaveletEngine::update_thresholds`, `compute_thresholds_internal` in `offline.cpp`, `src/adaptative.cpp` |
| R-side mode dispatch | `R/lwt.R`, `R/ilwt.R`, `R/denoising_offline.R`, `R/realtime_denoising.R` |
| XPtr lifetime | `create_engine_cpp` in `src/fast_stream.cpp` (`true` = enable finalizer), closure capture in `new_wavelet_stream` in `R/realtime_denoising.R` |

Detailed boundary-mode rationale in `04-boundary-and-threshold.md` Part A; threshold internals in Part B and in `02-adaptive-thresholding.md`.
