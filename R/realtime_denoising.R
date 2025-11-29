#' Create an Adaptive Wavelet Stream Processor (C++ Core)
#'
#' Generates a stateful function backed by a high-performance C++
#' Ring Buffer engine.
#' It implements Sliding Window + Lifting Decomposition +
#'  Adaptive Thresholding
#' in constant amortized time (O(1)) per sample.
#'
#' @param scheme A \code{lifting_scheme} object.
#' @param window_size Sliding window size (W). Must be > 8.
#' @param levels Decomposition levels (default 1).
#' @param alpha Threshold decay parameter (Eq 9).
#' @param beta Threshold gain factor (Eq 9).
#' @param method Shrinkage method: "hard", "soft", "semisoft".
#' @param extension Boundary handling ('symmetric', 'periodic', 'zero').
#' @param update_freq How often to recompute threshold statistics (default 1).
#'
#' @return A closure function \code{processor(new_sample)} that accepts
#' a single numeric value and returns the filtered value immediately.
#' @export
new_wavelet_stream = function(
    scheme,
    window_size = 256,
    levels = 1,
    alpha = 0.3,
    beta = 1.2,
    method = "semisoft",
    extension = "symmetric",
    update_freq = 1
) {

  if (window_size < 8) stop("window_size must be at least 8.")

  ext_int = switch(
    extension,
    "symmetric" = 1L, "periodic" = 2L, "zero" = 3L, 1L
  )

  engine_ptr = create_engine_cpp(
    scheme$steps,
    as.numeric(scheme$normalization),
    as.integer(levels),
    as.integer(window_size),
    as.integer(ext_int)
  )

  step_iter = 0

  processor = function(new_sample) {
    if (length(new_sample) != 1) {
      stop("Stream processor accepts only one sample at a time.")
    }

    if (is.na(new_sample) || is.infinite(new_sample)) {
      warning("Invalid sample (NA or Inf) received. Returning as is.")
      return(new_sample)
    }

    res = process_sample_cpp(
      engine_ptr,
      as.numeric(new_sample),
      as.numeric(alpha),
      as.numeric(beta),
      as.character(method),
      as.integer(update_freq),
      as.integer(step_iter)
    )

    step_iter <<- step_iter + 1
    return(res)
  }

  return(processor)
}

#' Causal Batch Denoising (Turbo Simulation)
#'
#' Processes a complete signal simulating the sequential arrival of data.
#' Uses the specialized C++ class \code{WaveletEngine} to perform causal
#' filtering efficiently on a historical dataset.
#'
#' @param signal Complete vector of the noisy signal.
#' @param scheme \code{lifting_scheme} object.
#' @param levels Decomposition levels.
#' @param window_size Window size.
#' @param alpha Threshold decay parameter (Eq 9).
#' @param beta Threshold gain factor (Eq 9).
#' @param method Thresholding method ("soft", "hard", "semisoft").
#' @param extension Boundary treatment ('symmetric', 'periodic').
#' @param update_freq Frequency of threshold updates.
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
    update_freq = 1
) {

  ext_int = switch(
    extension,
    "symmetric" = 1L, "periodic" = 2L, "zero" = 3L, 1L
  )

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
    as.integer(update_freq)
  )

  return(output)
}
