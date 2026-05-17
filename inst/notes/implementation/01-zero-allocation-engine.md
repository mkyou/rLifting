# rLifting — Zero-Allocation Engine

This document describes the computational decisions behind rLifting's performance: the single-pass offline design, the pre-allocated causal engine, the ring buffer architecture, and the trade-offs involved in each choice. The goal is to explain not only *what* was implemented but *why* each decision was necessary and what it costs.

---

## 1. The Performance Problem in R-Level Wavelet Denoising

A straightforward R implementation of wavelet denoising involves three sequential operations: a forward transform (`lwt`), thresholding (`threshold`), and an inverse transform (`ilwt`). Each step returns an R object; the next step receives it as input. This design is clean and composable, but it imposes overhead that accumulates across the pipeline:

- Each R function call carries SEXP boxing and unboxing — converting between R's internal representation and the C++ types used in computation.
- Intermediate objects (`lwt` output, threshold output) are heap-allocated in R's memory manager, which involves garbage collection pressure and cache misses.
- The `NumericVector` type used by Rcpp wraps an R SEXP: it is not a plain `std::vector<double>`, and passing it between functions incurs reference-count updates.

For offline denoising on a signal of length $n$, the R-level routing (`lwt → threshold → ilwt`) requires three separate C++ dispatches, three intermediate allocations, and the associated GC overhead. Benchmarks on the rLifting implementation confirm this: routing through R costs approximately 4.4× more than executing the equivalent logic in a single C++ call.

For causal denoising, the problem is more severe. The naive approach — applying `lwt`/`ilwt` on a sliding window of size $W$ at each of $T$ samples — incurs $O(T \cdot W)$ total work, because each window processes $W$ samples independently. For $T = 10\,000$ and $W = 256$, this means executing $10\,000$ full LWT+ILWT cycles, each on 256 samples. At rLifting's offline speed of ~51 µs per transform, the naive causal approach would take ~510 ms; the ring-buffer implementation completes in ~1 ms.

---

## 2. The Single-Pass Offline Design

The function `denoise_offline_cpp` (defined in `src/offline.cpp`) fuses the forward LWT, thresholding, and inverse LWT into a single C++ call. No R objects are created between steps; all intermediate state lives in `std::vector<double>` variables on the local stack or heap, managed by standard C++ destructors.

The function receives the signal as a `NumericVector` (Rcpp type), converts it to `std::vector<double>` once at entry, and returns a `NumericVector` at exit. These two conversions — one copy in, one copy out — are the only R-managed allocations in the entire pipeline. Everything between them is plain C++.

Within `denoise_offline_cpp`, the forward LWT allocates per-level `even` and `odd` vectors (`std::vector<double>`) and a `details` array to accumulate detail coefficients. For the irregular path, it additionally allocates `t_levels[0..levels]` to store per-level position vectors. These are local variables; they are freed when the function returns. This is a deliberate trade-off: offline denoising is a one-shot operation with no persistent state, so pre-allocation would provide no benefit — the allocation cost is paid once and is negligible compared to the $O(n \log n)$ LWT computation.

The thresholding step operates in-place on the `details` array. No temporary storage is needed.

The speedup of this design over R-level routing stems from two factors: eliminating intermediate SEXP allocations and their associated GC pressure, and keeping the working data in CPU cache across all three phases. With R-level routing, each phase may trigger GC between steps, evicting the working data from cache. The single-pass design keeps `even`, `odd`, and `details` resident in L2/L3 cache throughout.

---

## 3. The WaveletEngine: Pre-Allocation at Construction Time

For causal and stream processing, the single-pass design is insufficient: the state between samples — the ring buffer and its associated workspaces — must persist across calls. The `WaveletEngine` class (defined in `inst/include/rLifting/WaveletEngine.h`) was designed around a single principle: **all heap allocations happen in the constructor; none happen in the processing loop**.

The constructor receives the wavelet scheme, normalization, number of decomposition levels $L$, window size $W$, boundary mode, and an irregularity flag. It allocates the following:

**Ring buffer and linearization workspace:**
- `ring_buffer`: $W$ doubles — the circular sample buffer.
- `ring_buffer_t`: $W$ doubles — parallel time-position buffer, allocated only when `irregular = true`.
- `work_signal`: $W$ doubles — the linearized (chronological) view of `ring_buffer`.
- `work_t`: $W$ doubles — linearized view of `ring_buffer_t`, irregular mode only.

**LWT workspaces:** Two jagged arrays `work_approx` and `work_detail`, each with $L + 1$ levels. Level $j$ holds $\lceil W / 2^j \rceil$ doubles for `work_approx` and $\lfloor W / 2^j \rfloor$ for `work_detail`. The total allocation is:

$$\sum_{j=0}^{L} \left\lceil \frac{W}{2^j} \right\rceil + \sum_{j=0}^{L} \left\lfloor \frac{W}{2^j} \right\rfloor \approx 4W \text{ doubles}$$

which is roughly $4 \times 256 \times 8 = 8\,\text{KB}$ for $W = 256$, $L = 4$ — comfortably within L1 cache on modern hardware.

**Position workspaces** (irregular mode only): `work_t_approx` and `work_t_detail` mirror the structure of `work_approx` and `work_detail` for time positions.

**Threshold cache:** `current_lambdas`: $L$ doubles — one threshold per decomposition level, updated every `update_freq` samples.

After construction, `push_and_process` — the per-sample hot path — performs no `new`, no `malloc`, and no `resize` calls. Every write goes into a pre-allocated buffer.

---

## 4. The Ring Buffer

The ring buffer is a circular array of $W$ doubles with a single integer head pointer. Insertion is O(1):

```cpp
ring_buffer[head] = new_val;
head = (head + 1) % window_size;
if (count < window_size) count++;
```

The modular arithmetic keeps `head` within `[0, W)` without branching. The `count` variable tracks how many samples have been inserted; until `count == W`, the buffer is not yet full and the raw sample is returned without processing. This warm-up phase spans the first $W$ samples.

Linearization — converting the circular view to a chronologically ordered array — requires one pass of $W$ reads:

```cpp
for (int i = 0; i < window_size; i++)
    work_signal[i] = ring_buffer[(head + i) % window_size];
```

This is $O(W)$ per sample, unavoidable given that the LWT needs a contiguous ordered array. However, it is a simple stride-1 memory access pattern, which is highly cache-friendly and executes close to memory bandwidth limits on modern hardware. For $W = 256$, this copies 2 KB — well within L1 cache.

The ring buffer achieves $O(T)$ total work for a signal of $T$ samples: each sample requires one $O(W)$ linearization plus one $O(W \log W)$ LWT, but $W$ is a constant fixed at construction time, so the per-sample cost is $O(W \log W)$ and the total cost is $O(T \cdot W \log W)$ with $W$ fixed — linear in $T$.

---

## 5. Why offline.cpp Does Not Use WaveletEngine

A natural question is why `denoise_offline_cpp` does not instantiate a `WaveletEngine` internally. The answer is that the engine's design optimizes for a different access pattern.

`WaveletEngine` pre-allocates for a fixed window size $W$ and processes one sample at a time, maintaining a ring buffer. Offline denoising has access to the full signal of length $n$, which may be larger than any reasonable $W$. Processing it sample-by-sample through a `WaveletEngine` of size $n$ would produce the correct result but would (a) waste the pre-allocated workspaces on a single use, (b) incur $n$ ring buffer insertions and $n$ linearizations instead of one LWT on the full signal, and (c) compute a causal (windowed) threshold at each step rather than the global offline threshold — which is mathematically different and generally less accurate.

The offline design is intentionally simpler: allocate locally, process once, return. The engine design is for repeated, stateful processing with a fixed memory footprint.

---

## 6. The Hot Path: What Happens Per Sample

When `push_and_process` is called after warm-up, the execution sequence is:

1. **Insert**: write `new_val` into `ring_buffer[head]`; advance `head`.
2. **Linearize**: copy $W$ values from `ring_buffer` into `work_signal` in chronological order. If irregular, also linearize `ring_buffer_t` into `work_t`.
3. **Forward LWT**: $L$ levels of predict/update steps over pre-allocated `work_approx`/`work_detail` buffers. No allocations.
4. **Threshold update** (conditional on `update_freq`): recompute `current_lambdas` from `work_detail[0]` using `nth_element` for O(n) MAD.
5. **Shrinkage**: apply hard/soft/semisoft in-place on `work_detail[0..L-1]`.
6. **Inverse LWT**: reverse the predict/update steps; merge even/odd back into `work_approx[0]`.
7. **Return**: `work_approx[0][window_size - 1]` — the reconstructed value of the most recently inserted sample.

No heap allocations occur in steps 1–7 after warm-up. The `resize` calls inside the LWT levels of `WaveletEngine` (lines 160–163 of `WaveletEngine.h`) deserve a note: the approximation vectors (`work_approx[j+1]`) are pre-allocated to exactly the right size in the constructor and their `resize` calls are true no-ops. The detail vectors (`work_detail[j]`) are pre-allocated to `window_size` in the constructor and downsized to `window_size / 2` on the first call — downsizing a `std::vector` never allocates memory (it only updates the size bookkeeping), so no heap work occurs. After the first call the sizes stabilise and all subsequent `resize` calls are no-ops.

---

## 7. The XPtr Pattern for Stream Mode

`new_wavelet_stream` in R needs to return a closure that, on each call, invokes `push_and_process` on a persistent `WaveletEngine`. The engine cannot live in R's memory as an R object (it is not serializable and contains non-R types). Instead, `create_engine_cpp` heap-allocates a `WaveletEngine` via `new` and wraps the raw pointer in an `Rcpp::XPtr<WaveletEngine>`:

```cpp
WaveletEngine *engine = new WaveletEngine(steps, norm, levels, window_size, ext_mode, irreg, ll_k);
Rcpp::XPtr<WaveletEngine> ptr(engine, true);  // true = register finalizer
return ptr;
```

The second argument `true` registers a delete finalizer: when R's garbage collector collects the `XPtr` object (i.e., when the R closure goes out of scope), it calls `delete engine`. This prevents memory leaks without requiring manual cleanup.

The R closure captures `engine_ptr` (the `XPtr`) in its environment. Each call to `process_sample_cpp` dereferences the pointer and calls `push_and_process` directly — one pointer dereference and one virtual method call, with no R-managed memory involved. This is the minimal overhead achievable at the R/C++ boundary.

---

## 8. Benchmark Evidence

The design decisions above translate directly to measured performance. On a Doppler signal with $n = 1\,024$, CDF 5/3 wavelet, 4 decomposition levels, 50 Monte Carlo replications:

| Package | Median time | Speedup vs rLifting |
|:--------|------------:|--------------------:|
| rLifting | 51 µs | — |
| wavethresh | 1 925 µs | 37× slower |
| adlift | 3 047 ms | ~59 000× slower |
| nlt | 3 622 ms | ~70 000× slower |

The gap between rLifting and `wavethresh` (~37×) reflects the elimination of R-level allocations and the single-pass C++ design. Both packages use compiled code; the difference is architectural. `adlift` and `nlt` operate primarily in R with adaptive (signal-dependent) predict steps that cannot be vectorised in the same way, explaining the four-to-five order-of-magnitude difference.

The irregular path adds approximately 1.7× overhead relative to the regular path, attributable entirely to the Lagrange interpolation in each predict step (see `02-lifting-scheme-and-transform.md`). The pre-allocation strategy ensures this overhead is in computation, not memory management.

> **Note:** these figures are preliminary. The full benchmark suite (including irregular grids and all competing packages across all signal types) is in progress. This document will be updated when final results are available.
