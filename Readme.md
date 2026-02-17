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
provides amortized O(1) per-sample processing, while maintaining full
compatibility with standard offline denoising.

### Key features

- **Offline and causal denoising** in a single package — same wavelet
  scheme, same API, different processing guarantees.
- **Zero-allocation C++ engine** via Rcpp. All transforms, thresholding,
  and reconstruction run at native speed.
- **Ring-buffer stream processor** (`new_wavelet_stream`) for real-time
  applications: fixed memory, constant-time updates, microsecond latency.
- **No look-ahead bias.** Causal mode guarantees zero data leakage from
  future to past (verified via counterfactual leakage tests against
  `wavethresh`).
- **Adaptive thresholding** based on Liu et al. (2014): recursive MAD-based
  noise estimation with three shrinkage methods (hard, soft, semisoft).
- **Six built-in wavelets**: Haar, DB2, CDF 5/3, CDF 9/7 (JPEG 2000),
  DD4, and Lazy. Plus a simple API for custom wavelets via
  `lift_step()` + `custom_wavelet()`.
- **Diagnostic suite** (`diagnose_wavelet`): automatic verification of
  perfect reconstruction, vanishing moments, orthogonality, compact
  support, and shift sensitivity.

### Performance

Benchmarks against `wavethresh`, `adlift`, and `nlt` on a Doppler signal
(N = 1024, Haar, 50 simulations):

| Metric | rLifting | wavethresh | adlift | nlt |
|:---|---:|---:|---:|---:|
| Median time (ms) | < 1 | < 1 | ~600 | ~900 |
| MSE | 0.009 | 0.025 | 0.009 | 0.009 |

`rLifting` matches the reconstruction accuracy of the adaptive lifting
packages (`adlift`, `nlt`) while running orders of magnitude faster.
`wavethresh` is fast but uses a fixed global threshold, resulting in
higher MSE.

For causal processing, `rLifting`'s ring-buffer architecture is ~700×
faster than a naive sliding-window loop built on `wavethresh`.

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
  method = "semisoft"
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
  method = "semisoft"
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

Detailed vignettes are available after installation:

| Vignette | Topic |
|:---------|:------|
| `01_introduction` | Package overview, basic usage |
| `02_benchmark_offline` | Speed and MSE comparison with wavethresh, adlift, nlt |
| `03_benchmark_causal` | Causal speed benchmark + counterfactual leakage test |
| `04_realtime` | Real-time stream processing with latency analysis |
| `05_extensions` | Custom wavelets, diagnostics, thresholding, low-level API |

``` r
vignette("01_introduction", package = "rLifting")
```

### Roadmap

The following features are planned for future versions of `rLifting`:

- **Wavelet packets:** full support for Wavelet Packet Decomposition (WPD)
  for finer frequency resolution.
- **Irregular grids:** native handling of non-equispaced data (Second
  Generation Wavelets) without imputation.
- **Multivariate denoising:** joint denoising of correlated signals.

### License

MIT License.
