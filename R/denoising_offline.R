#' Denoising Offline (Nao-Causal / Global)
#'
#' Aplica a ondaleta em todo o sinal de uma vez. Usa estatisticas globais
#' para o calculo do threshold recursivo (Eq. 9).
#'
#' @param signal Sinal completo.
#' @param scheme Objeto lifting_scheme.
#' @param alpha Parametro de threshold recursivo.
#' @param beta Fator de escala do threshold.
#' @param method Metodo ('hard', 'soft', 'semisoft').
#' @param extension Modo de extensao ('symmetric', 'periodic', 'zero').
#'
#' @return Sinal filtrado.
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
