# rLifting — Modes of Operation

Detailed description of the four operation modes: internal mechanics, regular vs irregular paths, and trade-offs.

---

## 1. Quick Reference

| Mode | Function | Causal | State | C++ path |
|:-----|:---------|:-------|:------|:---------|
| Step-by-step | `lwt` / `compute_adaptive_threshold` / `threshold` / `ilwt` | No | No | One call per step |
| Offline | `denoise_signal_offline` | No | No | Single pass in `denoise_offline_cpp` |
| Causal Batch | `denoise_signal_causal` | Yes | Ephemeral | `run_causal_batch_cpp` → `WaveletEngine` |
| Stream | `new_wavelet_stream` | Yes | Persistent | `XPtr<WaveletEngine>` |

---

## 2. Step-by-step (Manual Pipeline)

The user calls `lwt` → `compute_adaptive_threshold` → `threshold` → `ilwt`. Each is an independent R function that validates inputs and dispatches to a dedicated C++ routine.

**Internal mechanics.** `lwt_cpp` performs a multi-level polyphase split (even/odd) followed by predict (P) and update (U) steps. `compute_thresholds_cpp` estimates noise variance from the finest-level details via a MAD estimator and computes per-level thresholds recursively. `threshold_*_cpp` applies hard, soft, or semisoft shrinkage in-place. `ilwt_cpp` reverses the P/U steps and merges even/odd subbands back.

**Regular vs irregular.** In the regular path, predict step coefficients are fixed. In the irregular path (`t` vector supplied to `lwt`/`ilwt`), time positions are stored per level during the forward pass and Lagrange interpolation replaces the fixed convolution in each predict step.

**What you gain.** Full visibility into intermediate representations — detail coefficients after decomposition, threshold values, thresholded coefficients before reconstruction. Custom shrinkage functions can be inserted between steps. Useful for diagnostics (`diagnose_wavelet`), research, or hybrid pipelines that apply different logic per level.

**What you lose.** Four separate R-level dispatch calls instead of one. Batch-only — no internal state, so causal processing would require managing a sliding window manually and calling `lwt`/`ilwt` on each slice, which replicates what `WaveletEngine` does but without its pre-allocated workspaces.

---

## 3. Offline (Non-causal)

`denoise_signal_offline` executes LWT + thresholding + ILWT in a single C++ function call (`denoise_offline_cpp`), with no R-level overhead between steps.

**Internal mechanics.** `denoise_offline_cpp` parses the lifting steps once, then runs the forward LWT level by level (polyphase split → P/U → normalisation), computes per-level thresholds from the finest-level details via the MAD estimator, applies shrinkage in-place, and reconstructs via the inverse lifting steps in reverse order.

**Regular vs irregular.** In the regular path, predict steps use fixed coefficients with boundary extension (`get_val_safe`). In the irregular path (`t` supplied), the function additionally maintains `t_levels[0..levels]` — per-level vectors of time positions. At each forward level it splits `t_levels[j]` into even and odd position vectors; at each inverse level it reads them back. Predict steps use `interp_predict` (Lagrange) instead of fixed convolution. This costs approximately 1.7× more time than the regular path at n = 1024.

**What you gain.** Fastest batch path: a single C++ round-trip regardless of signal length. The threshold is computed from the full set of finest-level detail coefficients, giving the most accurate noise estimate for a given signal.

**What you lose.** Non-causal — the output at any sample depends on the entire signal. Each call is independent; there is no persistent state.

---

## 4. Causal Batch

`denoise_signal_causal` simulates sample-by-sample processing on a complete historical signal. It calls `run_causal_batch_cpp`, which instantiates a `WaveletEngine` and feeds samples one at a time.

**Internal mechanics.** `WaveletEngine` maintains a circular ring buffer of `window_size` samples (forced to odd so that `n_even = n_odd + 1` holds at every decomposition level). For each new sample:

1. The sample is written into `ring_buffer[head]`; `head` advances modulo `window_size`.
2. Until the buffer is full (`count < window_size`), the raw sample is returned as-is.
3. Once full, the ring buffer is linearized into `work_signal` in chronological order.
4. The forward LWT runs over `work_signal` using pre-allocated workspace vectors (`work_approx`, `work_detail`).
5. If `step_iter % update_freq == 0`, thresholds are recomputed from `work_detail[0]` (finest-level details).
6. Shrinkage is applied in-place.
7. ILWT reconstructs into `work_approx[0]`.
8. The method returns `work_approx[0][window_size - 1]` — the reconstructed value of the most recently added sample.

All workspace vectors are allocated once in the constructor; no heap allocations occur inside `push_and_process`.

**Regular vs irregular.** In irregular mode, a parallel ring buffer `ring_buffer_t` stores each sample's time position. On linearisation, positions are split into even/odd arrays (`work_t_approx`, `work_t_detail`) during the forward pass and reused during the inverse. `run_causal_batch_cpp` receives the full `t` vector and slices it per sample internally.

**What you gain.** Causality: each output depends only on the current and previous `window_size - 1` samples. The threshold adapts over time via `update_freq`. Efficient for retrospective causal analysis — one `WaveletEngine` instantiation plus one `push_and_process` per sample.

**What you lose.** Output lags input by `window_size` samples (buffer must fill before the first filtered value is returned). Threshold estimates are derived from `window_size / 2` detail coefficients (the finest-level subband) — fewer than the offline estimate, so potentially noisier on short signals or at startup. State is ephemeral: the `WaveletEngine` is destroyed when the batch call returns.

---

## 5. Stream (Sample-by-sample)

`new_wavelet_stream` returns a stateful R closure backed by a persistent `WaveletEngine`. The engine is heap-allocated by `create_engine_cpp` and returned to R as an `Rcpp::XPtr<WaveletEngine>`. Each closure call invokes `process_sample_cpp`, which calls `push_and_process` on the pointer.

**Internal mechanics.** Identical to Causal Batch — the same `WaveletEngine::push_and_process` is called. The difference is lifecycle: the engine persists for the lifetime of the closure rather than a single batch call. `alpha`, `beta`, `method`, and `update_freq` are passed at each invocation, so they can be changed between samples without reconstructing the engine. The closure maintains a `step_iter` counter (incremented after each sample) that drives the `update_freq` check.

**Regular vs irregular.** In irregular mode (`irregular = TRUE`), the closure signature becomes `processor(new_sample, t_val)`, where `t_val` is the physical recording time of that sample. The engine stores it in a parallel ring buffer (`ring_buffer_t`) and uses it automatically during the forward and inverse Lagrange predict steps — the user only needs to supply the timestamp; everything else is internal. If `t_val` is omitted, `step_iter` (a sequential integer counter) is used as a fallback, which is only correct for uniformly-sampled data.

**What you gain.** Real-time capable: each call returns one filtered value with no batch overhead. Engine state (ring buffer, cached thresholds) persists across calls, so there is no startup cost after the first `window_size` samples. Parameters can be adjusted on-the-fly without reinitialising the engine.

**What you lose.** Same `window_size`-sample initial delay as Causal Batch. Thread safety is the caller's responsibility — concurrent access to the same closure is not safe. The `XPtr` becomes invalid if the closure is garbage-collected; using a raw copy of the pointer after that will crash.

---

## 6. Choosing a Mode

| Need | Recommended mode |
|:-----|:----------------|
| Inspect coefficients or prototype custom logic | Step-by-step |
| Best quality on a complete signal | Offline |
| Reproduce causal behaviour on historical data | Causal Batch |
| Process a live data stream | Stream |
