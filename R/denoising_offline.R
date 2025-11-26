#' Offline Denoising (Non-Causal / Global)
#'
#' Applies the wavelet to the entire signal at once. Uses global statistics
#' for recursive threshold calculation (Eq. 9).
#'
#' @param signal Complete signal.
#' @param scheme `lifting_scheme` object.
#' @param alpha Recursive threshold parameter.
#' @param beta Threshold scale factor.
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
    method = "semisoft",
    extension = "symmetric"
    ) {

  # 1. Forward global (com extensao definida pelo usuario)
  lwt_obj = lwt(signal, scheme, levels = 3, extension = extension)

  # 2. Calcular Thresholds adaptativos (baseado na estatistica global)
  # Ref  (Eq. 9)
  lambdas = compute_adaptive_threshold(lwt_obj, alpha, beta)

  # 3. Aplicar Thresholding
  for (lvl in names(lambdas)) {
    if (!is.null(lwt_obj$coeffs[[lvl]])) {
      lwt_obj$coeffs[[lvl]] = threshold(
        lwt_obj$coeffs[[lvl]], lambdas[[lvl]], method
        )
    }
  }

  # 4. Reconstrucao
  rec = ilwt(lwt_obj)
  return(rec[1:length(signal)])
}
