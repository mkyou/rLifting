## rLifting: high-performance wavelet lifting for R

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental) [![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT) **rLifting** is a versatile package for wavelet lifting transforms, bridging the gap between robust **offline analysis** (batch) and high-performance **causal processing** (online/streaming).

Unlike traditional packages that rely heavily on R-based loop structures, `rLifting` implements a hybrid architecture with a specialized C++ core, making it suitable for both functional data analysis (FDA) and production environments requiring low latency.

### 🚀 Key features

-   **Unified framework:** seamlessly switch between non-causal (global) and causal (windowed) denoising using the same wavelet schemes.
-   **High-performance core:** zero-allocation C++ engine (via Rcpp) optimized for speed.
-   **Causal & online:** implements a specialized ring buffer architecture for processing streaming data with fixed memory footprint.
-   **Adaptive thresholding:** recursive thresholding algorithms (based on Liu et al., 2014) for dynamic noise estimation.
-   **Custom wavelets:** simple API to construct bespoke wavelets defining predict/update steps.
-   **Reliability:** comprehensive unit testing (140+ tests) covering edge cases and mathematical reconstruction properties.

### ⚡ Performance

Benchmarks performed on `100,000` events simulation (window=1024, levels=3):

| Metric | Result | Description |
|:-----------------------|:-----------------------|:-----------------------|
| **Latency** | **\~11.9 μs** | Microseconds per sample injection (C++ core) |
| **Batch throughput** | **\~7,100** curves/sec | Full offline denoising (1024 points each) |
| **Stream throughput** | **\~23,000** events/sec | Online processing loop within R |

*Results indicate suitability for industrial sensor monitoring, financial time-series, and real-time signal cleaning.*

### 📦 Installation

You can install the development version from GitHub:

``` r
# install.packages("remotes") 
remotes::install_github("mkyou/rLifting")
```

### 🏁 Quick start

#### Offline denoising (global batch)

Ideal for historical data analysis where future points are known.
Uses global statistics for thresholding.

``` r
library(rLifting)

# 1. Generate synthetic data
x = rLifting:::.generate_signal("doppler", n = 1024)
x_noisy = x + rnorm(1024, sd = 0.2)

# 2. Define scheme (e.g., DB2, Haar, CDF97)
scheme = lifting_scheme("db2")

# 3. Denoise (global)
x_clean = denoise_signal_offline(
  x_noisy, 
  scheme, 
  levels = floor(log2(length(x_noisy))),
  method = "semisoft"
)

# Plot
plot(x_noisy, col = "grey", type = "l", main = "Offline denoising")
lines(x_clean, col = "blue", lwd = 2)
```

#### Causal denoising

Ideal for real-time scenarios.
Data is processed through a sliding window without "looking ahead".

``` r
# Initialize the causal engine (stateful)
# Window size = 256
processor = new_wavelet_stream(
  scheme, 
  window_size = 256, 
  levels = floor(log2(256))
)

# Simulate a stream (point-by-point)
stream_output = numeric(length(x_noisy))

for (i in seq_along(x_noisy)) {
  # Feed one sample, get one denoised sample immediately
  stream_output[i] = processor(x_noisy[i])
}

# Visualization: note the natural phase lag in causal filtering
lines(stream_output, col = "red", lwd = 2)
legend("topright", legend=c("Original", "Offline", "Causal"), 
       col=c("grey", "blue", "red"), lty=1)
```

### 📚 Documentation

Detailed vignettes are available (currently in Portuguese, English translation coming soon):

-   **Introduction:** `vignette("01_introduction", package = "rLifting")`
-   **Benchmarks (offline):** `vignette("02_benchmark_offline", package = "rLifting")`
-   **Benchmarks (causal):** `vignette("03_benchmark_causal", package = "rLifting")`
-   **Real-time processing:** `vignette("04_realtime", package = "rLifting")`
-   **Extensions:** `vignette("05_extensions", package = "rLifting")`

### 🗺️ Roadmap

Future versions of `rLifting` will focus on expanding support for complex data structures:

-   **Wavelet packets:** full support for Wavelet Packet Decomposition (WPD) for finer frequency resolution.
-   **Irregular grids:** native handling of non-equispaced data (Second Generation Wavelets) without the need for imputation.
-   **Multivariate denoising:** joint denoising of correlated signals.

### 📄 License

MIT License.
