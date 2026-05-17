#' Create an Adaptive Wavelet Stream Processor ('C++' Core)
#'
#' Generates a stateful function backed by a high-performance 'C++'
#' Ring Buffer engine.
#' It implements Sliding Window + Lifting Decomposition +
#'  Adaptive Thresholding
#' in highly efficient time per sample.
#'
#' @param scheme A \code{lifting_scheme} object.
#' @param window_size Sliding window size (W). Must be > 8.
#' @param levels Decomposition levels (default 1).
#' @param alpha Threshold decay parameter (Eq 9).
#' @param beta Threshold gain factor (Eq 9).
#' @param method Shrinkage method: "hard", "soft", "semisoft".
#' @param extension Boundary handling ('symmetric', 'periodic', 'zero', 'local_linear').
#' @param update_freq How often to recompute threshold statistics (default 1).
#' @param irregular Logical. If TRUE, the returned closure accepts a second
#'   argument \code{t_val} (the sample's time position) and applies
#'   position-aware interpolation in the predict steps.
#'
#' @return A closure \code{processor(new_sample, t_val = NULL)} that accepts
#' one sample (and optionally its time position) and returns the filtered value.
#' @export
new_wavelet_stream = function(
  scheme,
  window_size = 256,
  levels = 1,
  alpha = 0.3,
  beta = 1.2,
  method = "semisoft",
  extension = "symmetric",
  update_freq = 1,
  irregular = FALSE,
  ll_k = 4L
) {

  if (window_size < 8) stop("window_size must be at least 8.")
  if (window_size %% 2L == 0L) window_size = window_size + 1L
  if (extension == "local_linear" && ll_k > window_size)
    warning(sprintf("ll_k (%d) > window_size (%d): clamped to n.", ll_k, window_size))
  if (irregular) {
    if (extension == "one_sided")
      warning("extension 'one_sided' ignores irregular grid positions: Lagrange interpolation will not be applied. Use 'symmetric' or 'local_linear' for irregular-grid processing.")
    .check_irregular_scheme(scheme)
  }

  ext_int = switch(
    extension,
    "symmetric" = 1L, "periodic" = 2L, "zero" = 3L, "local_linear" = 4L, "one_sided" = 5L, 1L
  )

  engine_ptr = create_engine_cpp(
    scheme$steps,
    as.numeric(scheme$normalization),
    as.integer(levels),
    as.integer(window_size),
    as.integer(ext_int),
    as.logical(irregular),
    as.integer(ll_k)
  )

  step_iter = 0

  processor = function(new_sample, t_val = NULL) {
    if (length(new_sample) != 1) {
      stop("Stream processor accepts only one sample at a time.")
    }

    if (is.na(new_sample) || is.infinite(new_sample)) {
      warning("Invalid sample (NA or Inf) received. Returning as is.")
      return(new_sample)
    }

    t_cpp = if (is.null(t_val)) as.numeric(step_iter) else as.numeric(t_val)

    res = process_sample_cpp(
      engine_ptr,
      as.numeric(new_sample),
      t_cpp,
      as.numeric(alpha),
      as.numeric(beta),
      as.character(method),
      as.integer(update_freq),
      as.integer(step_iter)
    )

    step_iter <<- step_iter + 1
    return(res)
  }

  class(processor) = c("wavelet_stream", "function")
  attr(processor, "config") = list(
    wavelet    = scheme$wavelet,
    window_size = window_size,
    levels     = levels,
    method     = method,
    irregular  = irregular
  )

  return(processor)
}

#' Print method for Wavelet Stream Processor
#'
#' @param x Object of class \code{wavelet_stream}.
#' @param ... Additional arguments.
#' @return Invisibly returns \code{x}.
#' @export
print.wavelet_stream = function(x, ...) {
  cfg = attr(x, "config")
  cat("--- Wavelet Stream Processor ---\n")
  cat(sprintf("Wavelet: %s\n", cfg$wavelet))
  cat(sprintf("Window:  %d samples\n", cfg$window_size))
  cat(sprintf("Levels:  %d\n", cfg$levels))
  cat(sprintf("Method:  %s\n", cfg$method))
  invisible(x)
}

#' Causal Batch Denoising (Turbo Simulation)
#'
#' Processes a complete signal simulating the sequential arrival of data.
#' Uses the specialized 'C++' class \code{WaveletEngine} to perform causal
#' filtering efficiently on a historical dataset.
#'
#' @param signal Complete vector of the noisy signal.
#' @param scheme \code{lifting_scheme} object.
#' @param levels Decomposition levels.
#' @param window_size Window size.
#' @param alpha Threshold decay parameter (Eq 9).
#' @param beta Threshold gain factor (Eq 9).
#' @param method Thresholding method ("soft", "hard", "semisoft").
#' @param extension Boundary treatment ('symmetric', 'periodic', 'zero', 'local_linear').
#' @param update_freq Frequency of threshold updates.
#' @param t Optional numeric vector of sample time positions (irregular grid).
#'   Must be sorted and the same length as \code{signal}.
#'
#' @return Filtered vector (same length as input).
#' @export
denoise_signal_causal = function(
  signal,
  scheme,
  levels = 1,
  window_size = 256,
  alpha = 0.3,
  beta = 1.2,
  method = "semisoft",
  extension = "symmetric",
  update_freq = 1,
  t = NULL,
  ll_k = 4L
) {

  if (window_size %% 2L == 0L) window_size = window_size + 1L
  if (extension == "local_linear" && ll_k > window_size)
    warning(sprintf("ll_k (%d) > window_size (%d): clamped to n.", ll_k, window_size))

  if (!is.null(t)) {
    if (length(t) != length(signal)) stop("'t' must have the same length as 'signal'.")
    if (is.unsorted(t))              stop("'t' must be sorted in increasing order.")
    if (extension == "one_sided")
      warning("extension 'one_sided' ignores irregular grid positions: Lagrange interpolation will not be applied. Use 'symmetric' or 'local_linear' for irregular-grid processing.")
    .check_irregular_scheme(scheme)
  }

  ext_int = switch(
    extension,
    "symmetric" = 1L, "periodic" = 2L, "zero" = 3L, "local_linear" = 4L, "one_sided" = 5L, 1L
  )

  t_cpp = if (is.null(t)) numeric(0) else as.numeric(t)

  output = run_causal_batch_cpp(
    as.numeric(signal),
    scheme$steps,
    as.numeric(scheme$normalization),
    as.integer(levels),
    as.integer(window_size),
    as.numeric(alpha),
    as.numeric(beta),
    as.character(method),
    as.integer(ext_int),
    as.integer(update_freq),
    t_cpp,
    as.integer(ll_k)
  )

  return(output)
}
