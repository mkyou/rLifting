# 3. Causal decomposition: the rLifting advantage

Most wavelet packages in R are designed for offline analysis: they
assume the entire signal is available at once. To perform causal
filtering (where the output at time $t$ depends only on data at times
$0\ldots t$) with these packages, one must implement a sliding-window
loop, re-calculating the transform at every new point.

This approach has $O(N \cdot W)$ complexity (where $W$ is the window
size). `rLifting` implements a specialized causal mode using a
ring-buffer architecture, which updates the transform state in $O(N)$
(amortized $O(1)$ per point).

As in vignette 2, the results are pre-computed by the script
`data-raw/generate_vignette_data.R`.

## 1. Speed and accuracy: optimized vs. naive

We compare
[`rLifting::denoise_signal_causal`](../reference/denoise_signal_causal.md)
against a naive sliding-window implementation using `wavethresh`. Both
methods process a HeaviSine signal ($n = 500$, $\sigma = 0.3$) using the
Haar wavelet.

``` r
library(rLifting)

if (!requireNamespace("knitr", quietly = TRUE)) {
  knitr::opts_chunk$set(eval = FALSE)
  message("Required package 'knitr' is missing. Vignette code will not run.")
} else {
  library(knitr)
}

data("benchmark_causal", package = "rLifting")

df_causal = data.frame(
  Method = c("rLifting (causal)", "Naive loop (wavethresh)"),
  Time_ms = c(benchmark_causal$rLifting_Time_Avg * 1000,
              benchmark_causal$Wavethresh_Naive_Time * 1000),
  MSE = c(benchmark_causal$rLifting_MSE,
          benchmark_causal$Wavethresh_Naive_MSE)
)

kable(df_causal,
      col.names = c("Method", "Time (ms)", "MSE"),
      digits = 4,
      caption = "Causal denoising: speed and reconstruction accuracy.")
```

| Method                  | Time (ms) |    MSE |
|:------------------------|----------:|-------:|
| rLifting (causal)       |    3.7601 | 0.1094 |
| Naive loop (wavethresh) |  808.0936 | 0.2218 |

Causal denoising: speed and reconstruction accuracy.

`rLifting` is approximately 215× faster than the naive loop. Its MSE is
also competitive, demonstrating that the ring-buffer architecture does
not sacrifice accuracy for speed.

## 2. Leakage test (look-ahead bias)

A causal filter guarantees that the output at time $t$ depends only on
data up to and including time $t$. In offline processing, the transform
operates on the entire signal at once, so a future event (such as a
sudden step) can alter the denoised output at earlier time indices
through the wavelet coefficients that span the event boundary. This is
called data leakage, or look-ahead bias.

### Methodology

To quantify leakage, we use a counterfactual test. We generate two
signals that share identical noise for $t < 128$ but differ after:

- Signal A: noise only (no event)
- Signal B: same noise, but with a step of amplitude 5 added at
  $t \geq 128$

We run each filter on both signals and measure the total squared
difference in their output before the step ($t < 128$):

$$\text{Leakage} = \sum\limits_{t = 1}^{127}\lbrack{\widehat{f}}_{B}(t) - {\widehat{f}}_{A}(t)\rbrack^{2}$$

A truly causal filter produces identical output for both signals in
$\lbrack 1,127\rbrack$, so its leakage is exactly zero. An offline
filter may produce different output because the step’s wavelet
coefficients, which have support spanning the boundary, influence the
reconstruction at earlier indices.

To make this effect visible, we use longer-support wavelets (CDF 9/7 for
`rLifting`, Daubechies-8 for `wavethresh`) with `floor(log2(N))`
decomposition levels. Longer filters have wider support, so their
coefficients near the step boundary reach further into the past.

``` r
data("leakage_results", package = "rLifting")

kable(leakage_results,
      col.names = c("Method", "Leakage (SSE)"),
      caption = "Counterfactual leakage: lower is better. Zero means no look-ahead bias.")
```

| Method                     | Leakage (SSE) |
|:---------------------------|--------------:|
| rLifting causal (CDF 9/7)  |       0.00000 |
| rLifting offline (CDF 9/7) |       1.69639 |
| wavethresh offline (D8)    |      26.27047 |

Counterfactual leakage: lower is better. Zero means no look-ahead bias.

The causal mode of `rLifting` achieves exactly zero leakage: its output
before $t = 128$ is identical regardless of whether a step follows or
not. The filter simply has not seen the future data.

Both offline methods show non-zero leakage. The wavelet coefficients
near the step boundary ($t \approx 128$) have support that extends into
the past region, so the large discontinuity introduced by the step
affects the reconstructed values at earlier indices. This is an inherent
property of offline wavelet processing with non-trivial filter support —
the longer the filter, the further the leakage reaches.

This property is critical for financial backtesting (where look-ahead
bias leads to unrealistically optimistic results) and real-time control
systems (where the controller must not react to events that have not yet
occurred).
