#' Noisy Doppler Signal Example
#'
#' Synthetic Doppler signal contaminated with Gaussian noise. Used in
#' \code{vignette("v01-introduction")} and the boundary-mode comparison.
#'
#' @format A data frame with 2048 rows and 3 columns:
#' \describe{
#'   \item{index}{Time index (1..2048).}
#'   \item{original}{The pure Doppler signal.}
#'   \item{noisy}{The signal with added Gaussian noise (sd = 0.5).}
#' }
#' @usage data(doppler_example)
"doppler_example"

#' rLifting Offline/Causal Benchmark Results
#'
#' MSE and timing statistics for rLifting denoising across the four
#' Donoho-Johnstone signals, three modes (offline / causal-batch / stream),
#' multiple wavelets, boundary modes, threshold rules, and shrinkage methods.
#' Reported as 7-quartile summaries (min / q1 / median / mean / q3 / max / se)
#' over 1000 simulations per configuration. Used by
#' \code{vignette("v02-thresholding-and-tuning")} and
#' \code{vignette("v03-causal-stream")}.
#'
#' @format A data frame with 3600 rows and 40 columns:
#' \describe{
#'   \item{Signal}{Donoho-Johnstone test signal: \code{"blocks"},
#'     \code{"bumps"}, \code{"doppler"}, or \code{"heavisine"}.}
#'   \item{Pkg}{Always \code{"rLifting"}.}
#'   \item{Mode}{\code{"offline"}, \code{"causal"}, or \code{"stream"}.}
#'   \item{Wavelet}{Lifting scheme: \code{"haar"}, \code{"cdf53"}, etc.}
#'   \item{Boundary}{Boundary extension mode.}
#'   \item{Method}{Composite label combining threshold rule and shrinkage.}
#'   \item{ThresholdMethod}{Threshold rule: \code{"universal"} or
#'     \code{"sure"}.}
#'   \item{Shrinkage}{Shrinkage rule: \code{"hard"}, \code{"soft"},
#'     \code{"semisoft"}, or \code{"scad"}.}
#'   \item{AlphaUsed, BetaUsed}{Threshold-recursion parameters used.}
#'   \item{Version}{Generator version tag (\code{"v2"} for current).}
#'   \item{N}{Number of simulations per configuration (1000).}
#'   \item{MSE_min, MSE_q1, MSE_median, MSE_mean, MSE_q3, MSE_max,
#'     MSE_se}{Mean-squared-error statistics against the noise-free signal.}
#'   \item{MSE_settled_min, MSE_settled_q1, MSE_settled_median, MSE_settled_mean,
#'     MSE_settled_q3, MSE_settled_max, MSE_settled_se}{MSE statistics computed
#'     after dropping the warm-up window (causal/stream modes only).}
#'   \item{Time_total_min, Time_total_q1, Time_total_median, Time_total_mean,
#'     Time_total_q3, Time_total_max, Time_total_se}{Total wall time per call
#'     (seconds).}
#'   \item{Per_sample_us_min, Per_sample_us_q1, Per_sample_us_median,
#'     Per_sample_us_mean, Per_sample_us_q3, Per_sample_us_max,
#'     Per_sample_us_se}{Per-sample time (microseconds), i.e.
#'     Time_total / signal length.}
#' }
#' @source \code{data-raw/generate_rlifting_benchmark_v2.R}
#' @usage data(benchmark_rlifting)
"benchmark_rlifting"

#' wavethresh Benchmark Results
#'
#' MSE and timing statistics for the \pkg{wavethresh} package on the four
#' Donoho-Johnstone signals, across its native wavelets and boundary/threshold
#' combinations. 7-quartile summaries over 1000 simulations per configuration.
#' Used as a reference baseline in the offline benchmark vignette.
#'
#' @format A data frame with rows per (Signal, Wavelet, Boundary) and 19
#'   columns:
#' \describe{
#'   \item{Signal}{Donoho-Johnstone test signal.}
#'   \item{Pkg}{Always \code{"wavethresh"}.}
#'   \item{Wavelet}{Daubechies filter family identifier
#'     (e.g. \code{"co1"} for Daubechies-extremal-phase order 1).}
#'   \item{Boundary}{Combination of wavethresh threshold policy and shrinkage
#'     rule (e.g. \code{"BayesThresh_soft"}, \code{"cv_hard"}).}
#'   \item{N}{Number of simulations per configuration (1000).}
#'   \item{MSE_min, MSE_q1, MSE_median, MSE_mean, MSE_q3, MSE_max, MSE_se}{
#'     Reconstruction MSE statistics.}
#'   \item{Time_min, Time_q1, Time_median, Time_mean, Time_q3, Time_max,
#'     Time_se}{Wall-time statistics (seconds).}
#' }
#' @source \code{data-raw/generate_wavethresh_benchmark.R}
#' @usage data(benchmark_wavethresh)
"benchmark_wavethresh"

#' adlift Benchmark Results
#'
#' MSE and timing statistics for the \pkg{adlift} package (adaptive lifting on
#' irregular grids) on the four Donoho-Johnstone signals. 7-quartile summaries
#' over 1000 simulations per configuration. Used as a reference baseline in the
#' irregular-grid benchmark vignette.
#'
#' @format A data frame with 384 rows and 19 columns:
#' \describe{
#'   \item{Signal}{Donoho-Johnstone test signal.}
#'   \item{Pkg}{Always \code{"adlift"}.}
#'   \item{Wavelet}{Predictor family: \code{"AdaptPred"}, \code{"CubicPred"},
#'     \code{"LinearPred"}, or \code{"QuadPred"}.}
#'   \item{Boundary}{Encoded combination of \code{adlift::fwtnp} options
#'     (neighbours, interpolation, closest-point, mean/median predictor).}
#'   \item{N}{Number of simulations per configuration (1000).}
#'   \item{MSE_min, MSE_q1, MSE_median, MSE_mean, MSE_q3, MSE_max,
#'     MSE_se}{Reconstruction MSE statistics.}
#'   \item{Time_min, Time_q1, Time_median, Time_mean, Time_q3, Time_max,
#'     Time_se}{Wall-time statistics (seconds).}
#' }
#' @source \code{data-raw/generate_adlift_benchmark.R}
#' @usage data(benchmark_adlift)
"benchmark_adlift"

#' nlt Benchmark Results
#'
#' MSE and timing statistics for the \pkg{nlt} package (nondecimated lifting
#' transform) on the four Donoho-Johnstone signals. 7-quartile summaries over
#' 1000 simulations per configuration. Used as a reference baseline in the
#' irregular-grid benchmark vignette.
#'
#' @format A data frame with 384 rows and 19 columns:
#' \describe{
#'   \item{Signal}{Donoho-Johnstone test signal.}
#'   \item{Pkg}{Always \code{"nlt"}.}
#'   \item{Wavelet}{Predictor family inherited from \pkg{adlift}.}
#'   \item{Boundary}{Encoded combination of nlt/adlift options.}
#'   \item{N}{Number of simulations per configuration (1000).}
#'   \item{MSE_min, MSE_q1, MSE_median, MSE_mean, MSE_q3, MSE_max,
#'     MSE_se}{Reconstruction MSE statistics.}
#'   \item{Time_min, Time_q1, Time_median, Time_mean, Time_q3, Time_max,
#'     Time_se}{Wall-time statistics (seconds).}
#' }
#' @source \code{data-raw/generate_nlt_benchmark.R}
#' @usage data(benchmark_nlt)
"benchmark_nlt"

#' adlift Irregular-Grid Benchmark Results
#'
#' MSE and timing statistics for the \pkg{adlift} package on seven
#' irregular-grid test signals: three physically motivated (\code{linear_phys},
#' \code{trend_events}, \code{blocks_gapped}) and four Donoho-Johnstone classics
#' re-sampled on irregular grids (\code{blocks_dj_irr}, \code{bumps_dj_irr},
#' \code{doppler_dj_irr}, \code{heavisine_dj_irr}). 7-quartile summaries over
#' 1000 simulations per configuration. Companion to
#' \code{\link{benchmark_nlt_irregular}} and
#' \code{\link{benchmark_rlifting_irregular}}; used by the irregular-grid
#' benchmark vignette.
#'
#' @format A data frame with 672 rows and 20 columns:
#' \describe{
#'   \item{Signal}{One of the seven irregular-grid test signals
#'     (see Description).}
#'   \item{Pkg}{Always \code{"adlift"}.}
#'   \item{Wavelet}{Predictor family: \code{"AdaptPred"}, \code{"CubicPred"},
#'     \code{"LinearPred"}, or \code{"QuadPred"}.}
#'   \item{Boundary}{Encoded combination of \code{adlift::fwtnp} options.}
#'   \item{NoiseSd}{Per-signal Gaussian noise standard deviation
#'     (e.g. 0.15 for \code{linear_phys}/\code{trend_events}, 0.30 for the
#'     DJ-irregular signals, 0.50 for \code{blocks_gapped}).}
#'   \item{N}{Number of simulations per configuration (1000).}
#'   \item{MSE_min, MSE_q1, MSE_median, MSE_mean, MSE_q3, MSE_max,
#'     MSE_se}{Reconstruction MSE statistics.}
#'   \item{Time_min, Time_q1, Time_median, Time_mean, Time_q3, Time_max,
#'     Time_se}{Wall-time statistics (seconds).}
#' }
#' @source \code{data-raw/generate_adlift_irregular_benchmark.R}
#' @usage data(benchmark_adlift_irregular)
"benchmark_adlift_irregular"

#' nlt Irregular-Grid Benchmark Results
#'
#' MSE and timing statistics for the \pkg{nlt} package on seven irregular-grid
#' test signals (see \code{\link{benchmark_adlift_irregular}} for the signal
#' set). 7-quartile summaries over 1000 simulations per configuration.
#' Used by the irregular-grid benchmark vignette as a reference baseline.
#'
#' @format A data frame with 672 rows and 20 columns:
#' \describe{
#'   \item{Signal}{One of the seven irregular-grid test signals.}
#'   \item{Pkg}{Always \code{"nlt"}.}
#'   \item{Wavelet}{Predictor family inherited from \pkg{adlift}:
#'     \code{"AdaptPred"}, \code{"CubicPred"}, \code{"LinearPred"},
#'     or \code{"QuadPred"}.}
#'   \item{Boundary}{Encoded combination of nlt/adlift options.}
#'   \item{NoiseSd}{Per-signal Gaussian noise standard deviation.}
#'   \item{N}{Number of simulations per configuration (1000).}
#'   \item{MSE_min, MSE_q1, MSE_median, MSE_mean, MSE_q3, MSE_max,
#'     MSE_se}{Reconstruction MSE statistics.}
#'   \item{Time_min, Time_q1, Time_median, Time_mean, Time_q3, Time_max,
#'     Time_se}{Wall-time statistics (seconds).}
#' }
#' @source \code{data-raw/generate_nlt_irregular_benchmark.R}
#' @usage data(benchmark_nlt_irregular)
"benchmark_nlt_irregular"

#' rLifting Irregular-Grid Benchmark Results
#'
#' MSE and timing statistics for rLifting denoising on seven irregular-grid
#' test signals (see \code{\link{benchmark_adlift_irregular}} for the signal
#' set), covering all three modes (offline, causal, stream), six built-in
#' wavelets, five boundary modes, and twelve threshold/shrinkage method
#' combinations (universal/sure x hard/soft/semisoft/scad, with and without
#' \code{tune_alpha_beta} where applicable). 7-quartile summaries over 1000
#' simulations per configuration. Companion to
#' \code{\link{benchmark_adlift_irregular}} and
#' \code{\link{benchmark_nlt_irregular}}; used by the irregular-grid benchmark
#' vignette.
#'
#' Each configuration is run twice in offline mode: with position-aware
#' processing (\code{t = t_phys} passed in, irregular path active) and
#' position-ignoring (\code{t = NULL}, uniform-grid treatment). The
#' \code{Ratio_*} columns report \code{MSEpos / MSEfix} per simulation as a
#' direct measure of the value of irregular handling for that configuration.
#' Causal and stream modes report only position-aware results
#' (\code{MSEpos_*}/\code{Timepos_*}); the position-ignoring columns are
#' \code{NA} for those rows.
#'
#' @format A data frame with 7560 rows and 47 columns:
#' \describe{
#'   \item{Signal}{One of the seven irregular-grid test signals.}
#'   \item{Pkg}{Always \code{"rLifting"}.}
#'   \item{Mode}{\code{"offline"}, \code{"causal"}, or \code{"stream"}.}
#'   \item{Wavelet}{Lifting scheme: \code{"haar"}, \code{"db2"},
#'     \code{"cdf53"}, \code{"cdf97"}, \code{"dd4"}, or \code{"lazy"}.}
#'   \item{Boundary}{Boundary extension mode: \code{"symmetric"},
#'     \code{"periodic"}, \code{"zero"}, \code{"local_linear"}, or
#'     \code{"one_sided"}.}
#'   \item{Method}{Composite label combining threshold rule, shrinkage, and
#'     tuned/untuned status (e.g. \code{"universal_tuned_soft"}).}
#'   \item{ThresholdMethod}{Threshold rule: \code{"universal"} or
#'     \code{"sure"}.}
#'   \item{Shrinkage}{Shrinkage rule: \code{"hard"}, \code{"soft"},
#'     \code{"semisoft"}, or \code{"scad"}.}
#'   \item{AlphaUsed, BetaUsed}{Threshold-recursion parameters used (post-tuning
#'     when applicable).}
#'   \item{NoiseSd}{Per-signal Gaussian noise standard deviation.}
#'   \item{N}{Number of simulations per configuration (1000).}
#'   \item{MSEpos_min, MSEpos_q1, MSEpos_median, MSEpos_mean, MSEpos_q3,
#'     MSEpos_max, MSEpos_se}{MSE statistics with position-aware processing.}
#'   \item{MSEfix_min, MSEfix_q1, MSEfix_median, MSEfix_mean, MSEfix_q3,
#'     MSEfix_max, MSEfix_se}{MSE statistics with position-ignoring processing
#'     (offline only; \code{NA} for causal/stream).}
#'   \item{Ratio_min, Ratio_q1, Ratio_median, Ratio_mean, Ratio_q3, Ratio_max,
#'     Ratio_se}{Per-simulation ratio \code{MSEpos / MSEfix} summarised over
#'     the 1000 simulations (offline only).}
#'   \item{Timepos_min, Timepos_q1, Timepos_median, Timepos_mean, Timepos_q3,
#'     Timepos_max, Timepos_se}{Wall-time statistics with position-aware
#'     processing (seconds).}
#'   \item{Timefix_min, Timefix_q1, Timefix_median, Timefix_mean, Timefix_q3,
#'     Timefix_max, Timefix_se}{Wall-time statistics with position-ignoring
#'     processing (offline only).}
#' }
#' @source \code{data-raw/generate_rlifting_irregular_benchmark.R}
#' @usage data(benchmark_rlifting_irregular)
"benchmark_rlifting_irregular"
