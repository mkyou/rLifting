## rLifting: high-performance wavelet lifting for R

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental) [![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

`rLifting` is a high-performance R package for wavelet-based signal
denoising. It implements the Lifting Scheme (Sweldens, 1996) with a
zero-allocation C++ core via Rcpp, supporting both offline (batch) and
causal (real-time) processing modes through a unified API.

### Why rLifting?

Most wavelet packages in R are designed for offline use: they require the
entire signal before processing. To perform causal filtering (where the
output at time *t* depends only on data up to *t*), the user must write
a sliding-window loop that re-computes the transform at every point —
an O(N·W) operation.

`rLifting` solves this with a specialized ring-buffer architecture that
provides efficient O(T) processing for the complete signal (where T is the 
number of observations), while maintaining full compatibility with 
standard offline denoising.

### Key features

- **Offline and causal denoising** in a single package — same wavelet
  scheme, same API, different processing guarantees.
- **Zero-allocation C++ engine** via Rcpp. All transforms, thresholding,
  and reconstruction run at native speed.
- **Ring-buffer stream processor** (`new_wavelet_stream`) for real-time
  applications: fixed memory, high-speed updates, microsecond latency.
- **No look-ahead bias.** Causal mode guarantees zero data leakage from
  future to past (verified via counterfactual leakage tests against
  `wavethresh`).
- **Adaptive thresholding** based on Liu, Mi, & Mao (2014, *Meas. Sci. Rev.*
  14(3), DOI: 10.2478/msr-2014-0020): recursive MAD-based
  noise estimation with four shrinkage methods (hard, soft, semisoft, SCAD)
  and two threshold-selection rules (universal / VisuShrink and SureShrink).
  `tune_alpha_beta()` minimises SURE to pick the recursive parameters
  automatically; its benefit is signal-dependent (helps on signals with
  sharp features, can hurt on smooth ones, see `vignette("02-thresholding-and-tuning")`
  §5 for the empirical reality check).
- **Six built-in wavelets**: Haar, DB2, CDF 5/3, CDF 9/7 (JPEG 2000),
  DD4, and Lazy. Plus a simple API for custom wavelets via
  `lift_step()` + `custom_wavelet()`.
- **Irregular grids and edge-aware boundaries.** Pass a sorted `t` vector
  to any pipeline function to apply position-aware Lagrange interpolation
  in the predict steps. Five boundary modes: `symmetric`, `periodic`,
  `zero`, `local_linear` (linear extrapolation with neighbourhood size
  `ll_k`), and `one_sided` (asymmetric renormalised filters at the edge).
- **Diagnostic suite** (`diagnose_wavelet`): automatic verification of
  perfect reconstruction, vanishing moments, orthogonality, compact
  support, and shift sensitivity.

### Performance

Preliminary benchmark on a Doppler signal (N = 1024, CDF 5/3, 4 levels,
50 Monte Carlo replications):

| Package | Median time | Speedup vs rLifting |
|:---|---:|---:|
| rLifting | 51 µs | — |
| wavethresh | 1 925 µs | 37× |
| adlift | 3 047 ms | ~59 000× |
| nlt | 3 622 ms | ~70 000× |

`rLifting` matches the reconstruction accuracy of the adaptive lifting
packages (`adlift`, `nlt`) while running orders of magnitude faster.
`wavethresh` is fast but uses a fixed global threshold, resulting in
higher MSE.

For causal processing, `rLifting`'s ring-buffer architecture is ~700×
faster than a naive sliding-window loop built on `wavethresh`.

> Figures above are preliminary and consistent with
> `inst/notes/03-zero-allocation-engine.md`. The full
> benchmark suite (regular and irregular grids, all competing packages
> across all signal types) is running and will replace these once
> complete.

### Installation

``` r
# install.packages("remotes")
remotes::install_github("mkyou/rLifting")
```

### Quick start

#### Offline denoising

Process an entire signal at once using global statistics.

``` r
library(rLifting)

# Generate a noisy Doppler signal
x <- rLifting:::.generate_signal("doppler", n = 1024)
x_noisy <- x + rnorm(1024, sd = 0.2)

scheme <- lifting_scheme("cdf97")

x_clean <- denoise_signal_offline(
  x_noisy, scheme,
  levels = floor(log2(length(x_noisy))),
  shrinkage = "semisoft"
)

plot(x_noisy, col = "grey", type = "l", main = "Offline denoising")
lines(x_clean, col = "blue", lwd = 2)
lines(x, col = "black", lty = 2)
```

#### Causal denoising

Process a signal without using future data. Useful for financial
backtesting, real-time control, and streaming applications.

``` r
x_causal <- denoise_signal_causal(
  x_noisy, scheme,
  window_size = 256,
  levels = floor(log2(256)),
  shrinkage = "semisoft"
)
```

#### Real-time stream processing

Feed one sample at a time and get a denoised estimate immediately.

``` r
processor <- new_wavelet_stream(
  scheme, window_size = 256,
  levels = floor(log2(256))
)

stream_output <- numeric(length(x_noisy))
for (i in seq_along(x_noisy)) {
  stream_output[i] <- processor(x_noisy[i])
}
```

#### Custom wavelets

Define wavelets by specifying predict and update steps:

``` r
p <- lift_step("predict", coeffs = c(0.5, 0.5), start_idx = 0)
u <- lift_step("update", coeffs = c(0.25, 0.25), start_idx = -1)
my_wavelet <- custom_wavelet("MyCDF53", list(p, u), c(sqrt(2), 1/sqrt(2)))

# Works in any function
result <- denoise_signal_offline(x_noisy, my_wavelet, levels = 5)
```

### Documentation

The vignette series is organised as a learning path. The introduction (01)
forward-references every other vignette for the deep dives; this table
mirrors those promises.

| # | Vignette | Topic | Status |
|:--|:---------|:------|:------:|
| 01 | `introduction` | Package overview, three modes (offline, causal, stream), wavelet choice, basic parameters, irregular grids | available |
| 02 | `thresholding-and-tuning` | Why threshold parameters matter; MAD noise estimate; universal vs SURE; α/β recursion; `tune_alpha_beta()` walkthrough with empirical reality check; choice of shrinkage (hard / soft / semisoft / SCAD); decision guide separating bench-grounded recommendations from heuristics carried from the literature | available |
| 03 | `causal-stream` | Causality penalty across signals; `window_size` and `update_freq` heuristics; per-sample latency by mode and wavelet; causal-specific wavelet recommendations; leakage-check by construction; decision guide separating bench-grounded from heuristic recommendations | available |
| 04 | `boundary-modes` | Semantics of the five boundary extensions (`symmetric`, `periodic`, `zero`, `local_linear`, `one_sided`); `ll_k` for `local_linear`; empirical impact (small in offline, large in causal/stream); wavelet-dependent recommendations; `one_sided` constraints with irregular grids | available |
| 05 | `irregular-grids` | Lagrange interpolation in predict steps for non-uniform `t`; wavelet selection for irregular data; offline `t` vs stream `t_val` mechanisms | in progress |
| 06 | `extensions` | Built-in wavelet reference; custom wavelets via `lift_step()` + `custom_wavelet()`; diagnostics (`diagnose_wavelet` + standalone `validate_*` + `visualize_wavelet_basis`); low-level pipeline composing `lwt`/`ilwt`/`compute_adaptive_threshold`/`threshold`; brief overview of boundary modes and irregular-grid handling with cross-refs to vignettes 04 and 05 | available |
| 07 | `benchmarks` | Full mode-by-mode empirical comparison (offline / causal / stream × signals × wavelets × boundaries × threshold configs); speed and MSE vs `wavethresh`, `adlift`, `nlt`; ring-buffer speed-up over naive sliding-window | in progress |
| 08 | `case-study` | End-to-end denoising workflow on a real signal: diagnose, choose wavelet, tune α/β, denoise across the three modes, residual analysis | in progress |

``` r
vignette("01-introduction", package = "rLifting")
vignette("02-thresholding-and-tuning", package = "rLifting")
vignette("03-causal-stream", package = "rLifting")
vignette("04-boundary-modes", package = "rLifting")
vignette("06-extensions", package = "rLifting")
```

#### Internal design notes (`inst/notes/`)

Where the vignettes demonstrate behaviour with executable examples and
empirical decision guides, the `inst/notes/` docs specify formulas,
document architectural decisions, and serve as the reference for the C++
implementation. Five documents, ordered least technical to most technical:

| # | Doc | Scope | Status |
|:--|:----|:------|:------:|
| 00 | `design-overview` | Architectural index: package layers (R / C++), `lifting_scheme` abstraction, four modes, boundary handling, irregular grids, adaptive threshold; pointers out to subsystem docs and vignettes. | available |
| 01 | `lifting-scheme-and-transform` | The `lifting_scheme` representation and `LiftingStep` struct; R/C++ boundary; degree inference; built-in wavelets; polyphase decomposition; predict/update math; perfect reconstruction; irregular path inside the predict step. | available |
| 02 | `adaptive-thresholding` | Formal specification of the threshold subsystem: MAD estimator, universal and SURE rules, α/β recursion, four shrinkage formulas (hard / soft / semisoft / SCAD), offline vs causal behaviour, `tune_alpha_beta()` pipeline. | available |
| 03 | `zero-allocation-engine` | Why the C++ engine is designed this way: single-pass offline, pre-allocated `WaveletEngine`, ring buffer mechanics, XPtr finalizer for stream mode, per-sample hot path. | available |
| 04 | `boundary-and-threshold` | The four mandatory boundary-mode code paths (modes 1–4 vs mode 5); local-linear OLS; `get_t_extrap`; MAD via `nth_element`; shrinkage implementations. | available |

The two trees together form a two-axis matrix: vignettes (rows of the
matrix) tour the user-facing API, while the design notes (columns) sit
under each subsystem with the formal spec and implementation rationale.

### Roadmap

Based on the current development of the FATE estimator research, the following features are planned:

- **Multivariate denoising:** joint denoising of correlated signals to incorporate covariance structure.
- **Coefficient access:** expose intermediate LWT coefficients as a structured output (e.g. `result$coefficients`) from all pipeline functions, enabling inspection and custom post-processing without re-running the transform.

The technical documentation under `inst/notes/` covers the architectural decisions and implementation rationale for each subsystem (the lifting scheme and transform, the adaptive threshold, the zero-allocation engine, and the boundary-mode code paths). See the Documentation section above for the full layout.

> **Note:** the performance figures in this README and in `inst/notes/` will be revised once the full benchmark suite (regular and irregular grids, all competing packages) is complete.

### Next steps for CRAN submission

A static audit against the [CRAN Repository Policy](https://cran.r-project.org/web/packages/policies.html)
and the [submission checklist](https://cran.r-project.org/web/packages/submission_checklist.html)
flagged the following items to resolve before the next submission:

**Blockers (must fix before `R CMD build`):**

- ~~**Vignettes duplicated and orphan-referenced.**~~ Resolved: legacy
  `introduction.Rmd`, `extensions.Rmd`, `realtime.Rmd`, `benchmark_offline.Rmd`,
  and `benchmark_causal.Rmd` removed. `vignettes/` now contains only the
  numbered series (`01-introduction.Rmd` through `06-extensions.Rmd`, with
  `05-irregular-grids.Rmd` and `07-benchmarks.Rmd` pending).
- **`data/` polluted with benchmark checkpoints.** Add to `.Rbuildignore`:

  ```
  ^data/tmp_.*$
  ^data/.*\.v1backup$
  ```

  Without this, `R CMD build` would ship hundreds of MB of checkpoint
  `.rds` files, breaking the 5 MB data limit.
- **`DESCRIPTION:Date`** is stale (e.g. `2026-03-05`); bump to the
  submission date.
- **`DESCRIPTION:Version`** must be incremented past the last published
  version on CRAN before resubmission.

**Pre-submission tasks:**

- Refresh `cran-comments.md` to describe the new release's contents
  (SCAD shrinkage, SureShrink, `tune_alpha_beta`, irregular-grid
  support, boundary modes `local_linear` and `one_sided`).
- Add a new entry to `NEWS.md` covering the same.
- Run `R CMD check --as-cran` on a current R-devel and document
  `0 errors | 0 warnings | <N> notes` in `cran-comments.md`, with
  explanations for any remaining notes.
- Cross-check on win-builder and R-hub before submitting.

**Confirmed compliant** (no action needed):

- No `:::` access to other packages, no `.Internal`, no network calls,
  no writing outside `tempdir()`, no parallel use above the 2-core
  CRAN cap inside R or tests, no binaries shipped in the source tarball.
- LICENSE, NEWS.md, and cran-comments.md present.
- 35 man pages, none using `\dontrun{}`; long C++-only tests gated
  by `skip_on_cran()`.

### License

MIT License.

