#' Create an Adaptive Wavelet Stream Processor
#'
#' Generates a stateful function that implements the complete method from Liu et al.:
#' Moving Window + Lifting Decomposition + Recursive Adaptive Thresholding.
#'
#' @param scheme A `lifting_scheme` object.
#' @param window_size Window size (W). Must be a power of 2 (e.g., 128, 256)
#'        according to the paper, though the code accepts any integer > 8.
#' @param alpha Threshold decay parameter (Eq 9).
#' @param beta Threshold gain factor (Eq 9).
#' @param method Shrinkage method: "hard", "soft", "semisoft".
#' @param extension Boundary handling ('symmetric', 'periodic').
#'
#' @return A closure function `processor(new_sample)` that accepts a single
#' numeric value and returns the corresponding filtered value. The state (buffer)
#' is maintained between calls.
#' @export
new_wavelet_stream = function(
    scheme,
    window_size = 256,
    alpha = 0.3,
    beta = 1.2,
    method = "semisoft",
    extension = "symmetric"
) {

  # Validation: Minimum of 8 ensures numerical
  # stability when the window finally fills.
  if (window_size < 8) stop("window_size must be at least 8.")

  # Internal state (Closure)
  buffer = numeric(0)

  processor = function(new_sample) {
    # 1. Input Validation
    if (length(new_sample) != 1) {
      stop("Stream processor accepts only one sample at a time.")
    }
    if (is.na(new_sample) || is.infinite(new_sample)) {
      warning("Invalid sample (NA or Inf) received. Returning as is.")
      return(new_sample)
    }

    # 2. Update buffer
    buffer <<- c(buffer, new_sample)
    n_curr = length(buffer)

    # Maintain sliding window size
    if (n_curr > window_size) {
      excess = n_curr - window_size
      buffer <<- buffer[(excess + 1):n_curr]
      n_curr = window_size
    }

    # 3. Filling Phase (Liu et al. Implementation)
    # "At the first stage... when sample data are not long enough...
    # we keep the data as such." (Liu et al., 2014)
    # We only process when the buffer is FULL (n_curr == window_size).
    if (n_curr < window_size) {
      return(new_sample)
    }

    # 4. Processing (Protected)
    result = tryCatch({
      # Forward LWT
      # Now we are guaranteed to have 'window_size' samples.
      lwt_obj = lwt(buffer, scheme, levels = 1, extension = extension)

      # Dynamic Threshold Calculation
      lambdas = compute_adaptive_threshold(lwt_obj, alpha, beta)

      # Apply Thresholding
      for (lvl_name in names(lambdas)) {
        coeffs_old = lwt_obj$coeffs[[lvl_name]]
        lam = lambdas[[lvl_name]]
        coeffs_new = threshold(coeffs_old, lam, method)
        lwt_obj$coeffs[[lvl_name]] = coeffs_new
      }

      # Reconstruction
      rec = ilwt(lwt_obj)

      # Return only the last sample (Causal)
      val = rec[length(rec)]

      if (is.na(val)) return(new_sample)

      return(val)

    }, error = function(e) {
      # Fallback
      return(new_sample)
    })

    return(result)
  }

  return(processor)
}

#' Causal Batch Denoising (High Performance C++)
#'
#' Processes a complete signal simulating the sequential arrival of data.
#' Optimized with a C++ loop to avoid R overhead in simulations.
#'
#' @param signal Complete vector of the noisy signal.
#' @param scheme `lifting_scheme` object.
#' @param window_size Window size.
#' @param alpha Threshold decay parameter (Eq 9).
#' @param beta Threshold gain factor (Eq 9).
#' @param method Thresholding method ("soft", "hard", "semisoft").
#' @param extension Boundary treatment ('symmetric', 'periodic').
#' @param update_freq Integer. Frequency of threshold updates.
#' 1 = every sample (slowest, most accurate).
#' 10 = every 10 samples (faster, approximation).
#' Defaults to 1 for backward compatibility.
#'
#' @return Filtered vector (same length).
#' @export
denoise_signal_causal = function(
    signal,
    scheme,
    window_size = 256,
    alpha = 0.3,
    beta = 1.2,
    method = "semisoft",
    extension = "symmetric",
    update_freq = 1
) {

  ext_int = switch(extension, "symmetric" = 1L, "periodic" = 2L, "zero" = 3L, 1L)

  output = run_turbo_batch(
    as.numeric(signal),
    scheme$steps,
    as.numeric(scheme$normalization),
    as.integer(window_size),
    as.numeric(alpha),
    as.numeric(beta),
    as.character(method),
    as.integer(ext_int),
    as.integer(update_freq)
  )

  return(output)
}
