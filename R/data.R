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
#' Comparison of execution time and reconstruction MSE across packages, wavelets,
#' boundary modes, and Donoho-Johnstone signals. Pre-computed by
#' \code{data-raw/generate_vignette_data.R}.
#'
#' @format A data frame with 7 columns:
#' \describe{
#'   \item{Signal}{Test signal: \code{"doppler"}, \code{"heavisine"},
#'     \code{"bumps"}, \code{"blocks"}.}
#'   \item{Pkg}{Package: \code{"rLifting"}, \code{"wavethresh"},
#'     \code{"adlift"}, \code{"nlt"}.}
#'   \item{Wavelet}{Wavelet used (\code{NA} for non-rLifting packages).}
#'   \item{Boundary}{Boundary extension mode (\code{NA} for non-rLifting
#'     packages).}
#'   \item{Sim}{Simulation index.}
#'   \item{Time}{Execution time (seconds).}
#'   \item{MSE}{Mean Squared Error against the noise-free signal.}
#' }
#' @usage data(benchmark_offline)
"benchmark_offline"
