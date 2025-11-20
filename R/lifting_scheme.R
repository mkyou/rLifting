#' Construtor para o Esquema de Lifting (Lifting Scheme)
#'
#' Cria um objeto S3 contendo os passos de predicao e atualizacao,
#' alem dos fatores de normalizacao.
#'
#' @param wavelet Nome da ondaleta (string). Ex: "haar".
#' @param custom_steps Lista de passos (opcional, para uso futuro).
#' @param custom_norm Vetor de normalizacao (opcional).
#'
#' @return Um objeto da classe `lifting_scheme`.
#' @export
#'
#' @examples
#' scheme = lifting_scheme("haar")
#' print(scheme)
lifting_scheme = function(wavelet = "haar",
                          custom_steps = NULL,
                          custom_norm = NULL) {

  steps = list()
  norm_factors = c(1, 1)

  if (wavelet == "custom") {
    if (is.null(custom_steps)) stop("Para 'custom', forneca 'steps'.")
    steps = custom_steps
    if (!is.null(custom_norm)) norm_factors = custom_norm
  } else {
    config = .get_wavelet_config(wavelet)
    steps = config$steps
    norm_factors = config$norm
  }

  structure(
    list(wavelet = wavelet, steps = steps, normalization = norm_factors),
    class = "lifting_scheme"
  )
}

#' Helper interno com os coeficientes
#' @keywords internal
.get_wavelet_config = function(name) {

  if (name == "haar") {
    steps = list(
      list(type = "predict", coeffs = c(1), start_idx = 0),
      list(type = "update",  coeffs = c(0.5), start_idx = 0)
    )
    norm = c(sqrt(2), 1/sqrt(2))
    return(list(steps = steps, norm = norm))
  }

  if (name == "db2") {
    sqrt3 = sqrt(3)
    steps = list(
      # Predict: d[k] - P. Queremos prever usando a[k-1] e a[k].
      # Com start_idx = -1:
      # coeff[1] pega a[k-1] -> deve ser -(sqrt3-2)/4
      # coeff[2] pega a[k]   -> deve ser -sqrt3/4
      list(
        type = "predict",
        coeffs = c(-(sqrt3-2)/4, -sqrt3/4),
        start_idx = -1
        ),
      # Update: depende apenas do proprio d
      list(type = "update", coeffs = c(sqrt3), start_idx = 0)
    )
    norm = c((sqrt3 + 1) / sqrt(2), (sqrt3 - 1) / sqrt(2))
    return(list(steps = steps, norm = norm))
  }

  if (name == "cdf97" || name == "bior4.4") {
    alpha = -1.586134342
    beta  = -0.05298011854
    gamma = 0.8829110762
    delta = 0.4435068522
    zeta  = 1.149604398

    steps = list(
      # P1: s[k] + s[k+1]. start_idx=0 (pega k e k+1)
      list(type = "predict", coeffs = c(-alpha, -alpha), start_idx = 0),

      # U1: d[k] + d[k-1]. start_idx=-1 (pega k-1 e k)
      list(type = "update", coeffs = c(beta, beta), start_idx = -1),

      # P2: s[k] + s[k+1]
      list(type = "predict", coeffs = c(-gamma, -gamma), start_idx = 0),

      # U2: d[k] + d[k-1]
      list(type = "update", coeffs = c(delta, delta), start_idx = -1)
    )
    norm = c(zeta, 1/zeta)
    return(list(steps = steps, norm = norm))
  }
  stop(paste("Ondaleta", name, "nao encontrada."))
}

#' @export
print.lifting_scheme = function(x, ...) {
  cat(sprintf("Lifting Scheme: %s\n", x$wavelet))
  cat(sprintf("Steps: %d\n", length(x$steps)))
  cat(sprintf("Norm:  Approx=%.4f, Detail=%.4f\n",
              x$normalization[1], x$normalization[2]))
}
