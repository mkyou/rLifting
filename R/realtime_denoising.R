#' Create an Adaptive Wavelet Stream Processor (C++ Core)
#'
#' Generates a stateful function backed by a high-performance C++ engine.
#' It implements Sliding Window + Lifting Decomposition +
#' Adaptive Thresholding.
#'
#' @param scheme A `lifting_scheme` object.
#' @param window_size Window size (W). Must be > 8.
#' @param levels Decomposition levels (default 1).
#' @param alpha Threshold decay parameter (Eq 9).
#' @param beta Threshold gain factor (Eq 9).
#' @param method Shrinkage method: "hard", "soft", "semisoft".
#' @param extension Boundary handling ('symmetric', 'periodic', 'zero').
#' @param update_freq How often to recompute threshold stats (default 1).
#'
#' @return A closure function `processor(new_sample)` that accepts a single
#' numeric value and returns the filtered value.
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

  # Mapeia extensao para inteiro (C++)
  ext_int = switch(
    extension,
    "symmetric" = 1L, "periodic" = 2L, "zero" = 3L, 1L
  )

  # 1. Cria o Motor C++ (Ponteiro XPtr)
  # O estado (buffer, workspaces) vive na memoria do C++
  engine_ptr = create_engine_cpp(
    scheme$steps,
    as.numeric(scheme$normalization),
    as.integer(levels),
    as.integer(window_size),
    as.integer(ext_int)
  )

  # Contador interno para controlar update_freq (passado ao C++)
  step_iter = 0

  # 2. Closure (Wrapper R)
  processor = function(new_sample) {
    # Validacao de entrada (Robustez)
    if (length(new_sample) != 1) {
      # Mantemos retorno safe, mas avisamos se o teste exigir ou
      # se for input estruturalmente errado.
      stop("Stream processor accepts only one sample at a time.")
    }

    # Verifica NA ou Inf (Requisito dos testes de robustez)
    if (is.na(new_sample) || is.infinite(new_sample)) {
      warning("Invalid sample (NA or Inf) received. Returning as is.")
      return(new_sample)
    }

    # Chama o processamento de 1 amostra no C++
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
#' Uses the specialized C++ class `WaveletEngine`.
#'
#' @param signal Complete vector of the noisy signal.
#' @param scheme `lifting_scheme` object.
#' @param levels Decomposition levels.
#' @param window_size Window size.
#' @param alpha Threshold decay parameter (Eq 9).
#' @param beta Threshold gain factor (Eq 9).
#' @param method Thresholding method ("soft", "hard", "semisoft").
#' @param extension Boundary treatment ('symmetric', 'periodic').
#' @param update_freq Frequency of threshold updates.
#'
#' @return Filtered vector (same length).
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

  # Chama a versao Turbo (Batch) do C++
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
