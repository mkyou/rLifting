#' Noisy Doppler Signal Example
#'
#' A synthetic dataset containing a Doppler signal contaminated with Gaussian noise.
#' Used in the "General Usage" vignette.
#'
#' @format A data frame with 2048 rows and 3 columns:
#' \describe{
#'   \item{index}{Time index.}
#'   \item{original}{The pure Doppler signal.}
#'   \item{noisy}{The signal with added Gaussian noise (sd=0.5).}
#' }
#' @usage data(doppler_example)
"doppler_example"

#' Offline Benchmark Results
#'
#' Comparison of execution time and reconstruction error (MSE) between `rLifting` and other packages
#' (wavethresh, wavelets) using Haar wavelet.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{Pkg}{Package name.}
#'   \item{Time}{Execution time in seconds.}
#'   \item{MSE}{Mean Squared Error.}
#' }
#' @usage data(benchmark_offline)
"benchmark_offline"

#' Causal Benchmark Results
#'
#' Comparison of execution time between `rLifting`'s optimized causal mode and 
#' a naive sliding-window implementation using `wavethresh`.
#'
#' @format A list containing:
#' \describe{
#'   \item{rLifting_Time_Avg}{Average time (seconds) for rLifting.}
#'   \item{Wavethresh_Naive_Time}{Time (seconds) for naive sliding window.}
#'   \item{Speedup_Factor}{Ratio of Naive Time to rLifting Time.}
#' }
#' @usage data(benchmark_causal)
"benchmark_causal"

#' Leakage (Impulse Response) Results
#'
#' Measurement of energy leakage into the "past" when processing an impulse signal.
#' Used to demonstrate the zero-lookahead property of the causal mode.
#'
#' @format A data frame with:
#' \describe{
#'   \item{Method}{Method description (e.g. "rLifting causal").}
#'   \item{Leakage}{Sum of squared differences (leakage energy).}
#' }
#' @usage data(leakage_results)
"leakage_results"

