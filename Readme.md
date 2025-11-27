## rLifting: High-Performance Wavelet Lifting for R

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental) [![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT) **rLifting** is a versatile package for Wavelet Lifting Transforms, bridging the gap between robust **Offline Analysis** (Batch) and high-performance **Causal Processing** (Online/Streaming)
.

Unlike traditional packages that rely heavily on R-based loop structures, `rLifting` implements a hybrid architecture with a specialized C++ core, making it suitable for both functional data analysis (FDA) and production environments requiring low latency.

### 🚀 Key Features

-   **Unified Framework:** seamlessly switch between non-causal (global) and causal (windowed) denoising using the same wavelet schemes.
-   **High-Performance Core:** zero-allocation C++ engine (via Rcpp) optimized for speed.
-   **Causal & Online:** implements a specialized ring buffer architecture for processing streaming data with fixed memory footprint.
-   **Adaptive Thresholding:** recursive thresholding algorithms (based on Liu et al., 2014) for dynamic noise estimation.
-   **Custom Wavelets:** simple API to construct bespoke wavelets defining Predict/Update steps.
-   **Reliability:** comprehensive unit testing (140+ tests) covering edge cases and mathematical reconstruction properties.

### ⚡ Performance

Benchmarks performed on `100,000` events simulation (Window=1024, Levels=3):

| Metric | Result | Description |
|:-----------------------|:-----------------------|:-----------------------|
| **Latency** | **\~11.9 μs** | Microseconds per sample injection (C++ Core) |
| **Batch Throughput** | **\~7,100** curves/sec | Full offline denoising (1024 points each) |
| **Stream Throughput** | **\~23,000** events/sec | Online processing loop within R |

*Results indicate suitability for industrial sensor monitoring, financial time-series, and real-time signal cleaning.*

### 📦 Installation

You can install the development version from GitHub:

``` r
# install.packages("remotes") 
remotes::install_github("mkyou/rLifting")
```

### 🏁 Quick Start

#### Offline Denoising (Global Batch)

Ideal for historical data analysis where future points are known.
Uses global statistics for thresholding.

``` r
library(rLifting)

# 1. Generate synthetic data
x = rLifting:::.generate_signal("doppler", n = 1024)
x_noisy = x + rnorm(1024, sd = 0.2)

# 2. Define scheme (e.g., DB2, Haar, CDF97)
scheme = lifting_scheme("db2")

# 3. Denoise (Global)
x_clean = denoise_signal_offline(
  x_noisy, 
  scheme, 
  levels = 3,
  method = "semisoft"
)

# Plot
plot(x_noisy, col = "grey", type = "l", main = "Offline Denoising")
lines(x_clean, col = "blue", lwd = 2)
```

#### Causal Denoising

Ideal for real-time scenarios.
Data is processed through a sliding window without "looking ahead".

``` r
# Initialize the causal engine (Stateful)
# Window size = 256, Decomposition Levels = 3
processor = new_wavelet_stream(
  scheme, 
  window_size = 256, 
  levels = 3
)

# Simulate a stream (point-by-point)
stream_output = numeric(length(x_noisy))

for (i in seq_along(x_noisy)) {
  # Feed one sample, get one denoised sample immediately
  stream_output[i] = processor(x_noisy[i])
}

# Visualization: Note the natural phase lag in causal filtering
lines(stream_output, col = "red", lwd = 2)
legend("topright", legend=c("Original", "Offline", "Causal"), 
       col=c("grey", "blue", "red"), lty=1)
```

### 📚 Documentation

Detailed vignettes are available (currently in Portuguese, English translation coming soon):

-   General Usage Guide: `vignette("guide_usage", package = "rLifting")`

-   Statistical Benchmarks: `vignette("statistical_benchmark", package = "rLifting")`

-   Performance & Throughput: `vignette("performance_benchmark", package = "rLifting")`

### 📄 License

MIT License.
