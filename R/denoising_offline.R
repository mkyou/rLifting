#' Offline Denoising (Non-Causal / Global)
#'
#' Applies the wavelet to the entire signal at once. Uses global statistics
#' for recursive threshold calculation (Eq. 9).
#' This function is fully optimized in C++.
#'
#' @param signal Complete signal.
#' @param scheme `lifting_scheme` object.
#' @param alpha Recursive threshold parameter.
#' @param beta Threshold scale factor.
#' @param levels Decomposition levels.
#' @param method Method ('hard', 'soft', 'semisoft').
#' @param extension Extension mode ('symmetric', 'periodic', 'zero').
#'
#' @return Filtered signal.
#' @export
denoise_signal_offline = function(
    signal,
    scheme,
    alpha = 0.3,
    beta = 1.2,
    levels = 3,
    method = "semisoft",
    extension = "symmetric"
) {

  # Mapeia extensao para inteiro (C++)
  ext_int = switch(
    extension,
    "symmetric" = 1L,
    "periodic"  = 2L,
    "zero"      = 3L,
    1L
  )

  # Chama pipeline unificado
  # A funcao C++ faz LWT -> Threshold -> ILWT internamente.
  res = denoise_offline_cpp(
    as.numeric(signal),
    scheme$steps,
    as.numeric(scheme$normalization),
    as.integer(levels),
    as.numeric(alpha),
    as.numeric(beta),
    as.character(method),
    as.integer(ext_int)
  )

  return(res)
}
