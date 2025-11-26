#' Lifting Scheme Constructor
#'
#' Creates an S3 object containing prediction and update steps.
#'
#' @param wavelet Wavelet name (string). E.g., "haar", "db2".
#' @param custom_steps List of steps (if provided, ignores internal lookup).
#' @param custom_norm Normalization vector (optional).
#'
#' @return An object of class `lifting_scheme`.
#' @export
lifting_scheme = function(
    wavelet = "haar",
    custom_steps = NULL,
    custom_norm = NULL
    ) {

  steps = list()
  norm_factors = c(1, 1)

  # Prioridade para custom_steps. Se eles existirem, usamos eles.
  # O nome 'wavelet' vira apenas um label.
  if (!is.null(custom_steps)) {
    steps = custom_steps
    if (!is.null(custom_norm)) norm_factors = custom_norm
  } else {
    # Se nao tem passos customizados, tenta carregar do banco interno
    config = .get_wavelet_config(wavelet)
    steps = config$steps
    norm_factors = config$norm
  }

  structure(
    list(wavelet = wavelet, steps = steps, normalization = norm_factors),
    class = "lifting_scheme"
  )
}

#' Implements factorizations based on Daubechies & Sweldens (1998).
#' @keywords internal
.get_wavelet_config = function(name) {
  # Baseline
  # Apenas split, sem processamento.
  # Util para debug e ensino. Ref: Daubechies & Sweldens (1998) Sec 3.
  if (name == "lazy") {
    return(list(steps = list(), norm = c(1, 1)))
  }

  if (name == "haar") {
    # Haar simples: d = odd - even; a = even + 0.5*d
    steps = list(
      list(type = "predict", coeffs = c(1), start_idx = 0),
      list(type = "update",  coeffs = c(0.5), start_idx = 0)
    )
    norm = c(sqrt(2), 1/sqrt(2))
    return(list(steps = steps, norm = norm))
  }

  # --- DB2 (D4 - 4 taps) ---
  if (name == "db2") {
    sqrt3 = sqrt(3)
    steps = list(
      list(type = "predict", coeffs = c(sqrt3), start_idx = 0),
      list(
        type = "update",
        coeffs = c(sqrt3/4, (sqrt3 - 2)/4), start_idx = 0
        ),
      list(type = "predict", coeffs = c(-1), start_idx = -1)
    )
    norm = c((sqrt3 + 1) / sqrt(2), (sqrt3 - 1) / sqrt(2))
    return(list(steps = steps, norm = norm))
  }

  # --- CDF 5/3 (LeGall) ---
  # Inteira, Rapida, Simetrica. O "Must Have" para real-time leve.
  # P: d[k] -= 0.5 * (s[k] + s[k+1])
  # U: s[k] += 0.25 * (d[k] + d[k-1])
  if (name == "cdf53" || name == "bior2.2") {
    steps = list(
      # Predict: Coeffs -0.5, -0.5.
      # (Engine faz 'odd - P', entao -0.5 vira soma)
      # Porem, CDF 5/3 padrao e d = d - floor(0.5*(s0+s1)).
      # Se assumirmos float, coeficientes sao 0.5.
      # Engine: d = d - P. P = 0.5 * s.
      list(type = "predict", coeffs = c(0.5, 0.5), start_idx = 0),

      # Update: s = s + U. U = 0.25 * d.
      list(type = "update", coeffs = c(0.25, 0.25), start_idx = -1)
    )
    # Biorrtogonal pura (sem normalizacao irracional)
    # Para Denoising (preservar energia), usamos sqrt(2)
    norm = c(sqrt(2), 1/sqrt(2))
    # Obs: Para compressao lossless inteira, norm seria c(1,1).
    return(list(steps = steps, norm = norm))
  }

  if (name == "cdf97" || name == "bior4.4") {
    # Coeficientes padrao (Cohen-Daubechies-Feauveau 9/7)
    # Valores de Daubechies & Sweldens (1998), Secao 7.7
    alpha = -1.586134342
    beta  = -0.05298011854
    gamma = 0.8829110762
    delta = 0.4435068522
    zeta  = 1.149604398

    # CORRECAO DE SINAIS:
    # A engine lwt implementa PREDICT como SUBTRACAO (odd - P).
    # A ondaleta CDF 9/7 define os passos como ADICOES sequenciais.
    # Portanto, para passos PREDICT, invertemos o sinal do coeficiente.
    # Para passos UPDATE, mantemos o sinal (pois engine usa soma).

    steps = list(
      # P1: d <- d + alpha * (s[k] + s[k+1])
      # Engine: d <- d - P. Logo P deve ser -alpha.
      list(type = "predict", coeffs = c(-alpha, -alpha), start_idx = 0),

      # U1: s <- s + beta * (d[k] + d[k-1])
      # Engine: s <- s + U. Logo U deve ser beta.
      list(type = "update",  coeffs = c(beta, beta), start_idx = -1),

      # P2: d <- d + gamma * (s[k] + s[k+1])
      # Engine: d <- d - P. Logo P deve ser -gamma.
      list(type = "predict", coeffs = c(-gamma, -gamma), start_idx = 0),

      # U2: s <- s + delta * (d[k] + d[k-1])
      # Engine: s <- s + U. Logo U deve ser delta.
      list(type = "update",  coeffs = c(delta, delta), start_idx = -1)
    )
    norm = c(zeta, 1/zeta)
    return(list(steps = steps, norm = norm))
  }

  # Interpolacao suave usando 4 vizinhos.
  # Coeficientes: -1/16, 9/16, 9/16, -1/16.
  # Ref: Daubechies & Sweldens (1998) Sec 7.4.
  if (name == "dd4" || name == "interp4") {
    p_coeffs = c(-1/16, 9/16, 9/16, -1/16)
    # Update para manter momentos nulos (simetrico ao predict / 2)
    u_coeffs = p_coeffs / 2

    steps = list(
      list(type = "predict", coeffs = p_coeffs, start_idx = -1),
      list(type = "update",  coeffs = u_coeffs, start_idx = -1)
    )
    # Ondaletas interpoladoras geralmente nao sao ortogonais
    # Normalizacao para conservacao de energia aproximada
    norm = c(sqrt(2), 1/sqrt(2))
    return(list(steps = steps, norm = norm))
  }

  stop(paste("Ondaleta", name, "nao suportada."))
}

#' Print method
#' @param x lifting_scheme object.
#' @param ... additional arguments.
#' @export
print.lifting_scheme = function(x, ...) {
  cat(sprintf("Lifting Scheme: %s\n", x$wavelet))
  cat(sprintf("Steps: %d\n", length(x$steps)))
  cat(sprintf("Norm:  Approx=%.4f, Detail=%.4f\n",
              x$normalization[1], x$normalization[2]))
}


#' Create an individual Lifting Step
#'
#' Helper function to create prediction (P) or update (U) steps,
#' abstracting the complexity of index management.
#'
#' @param type Step type: "predict" (P) or "update" (U).
#' @param coeffs Numeric vector containing the filter coefficients.
#' @param start_idx (Optional) Manual start index. If provided, ignores the
#'        `position` parameter. Use this for fine-grained control.
#' @param position Automatic index adjustment (used only if `start_idx` is NULL):
#'        \itemize{
#'          \item "center": Centers the filter (default).
#'          \item "left": Causal filter (looks into the past).
#'          \item "right": Anti-causal filter (looks into the future).
#'        }
#'
#' @return A list formatted for the internal lifting engine.
#' @export
lift_step = function(
    type = c("predict", "update"),
    coeffs,
    start_idx = NULL,
    position = "center"
) {

  type = match.arg(type)
  n = length(coeffs)

  # Logica de alinhamento automatico se start_idx nao for passado
  if (is.null(start_idx)) {
    if (position == "center") {
      # Centraliza o filtro.
      # Ex: n=3 (-1, 0, 1) -> start_idx = -1
      # Ex: n=2 (0, 1) [forward] -> start_idx = 0
      start_idx = -floor((n - 1) / 2)
    } else if (position == "left") {
      # Filtro causal estrito (apenas passado)
      # Ex: n=2 usando k-1, k -> start_idx = -1
      start_idx = -n + 1
    } else {
      # position == "right"
      # Filtro anti-causal (apenas futuro)
      start_idx = 0
    }
  }

  list(type = type, coeffs = coeffs, start_idx = start_idx)
}

#' Create a custom wavelet
#'
#' Wrapper to create a `lifting_scheme` object from manual steps.
#'
#' @param name Identifier name for the wavelet.
#' @param steps List of steps created via \code{lift_step}.
#' @param norm Normalization vector c(K, 1/K).
#'
#' @return An object of class `lifting_scheme`.
#' @export
#' @examples
#' p1 = lift_step("predict", c(1), position="center")
#' u1 = lift_step("update", c(0.5), position="center")
#' w = custom_wavelet("HaarManual", list(p1, u1), c(1.41, 0.707))
custom_wavelet = function(name, steps, norm = c(1, 1)) {
  structure(
    list(wavelet = name, steps = steps, normalization = norm),
    class = "lifting_scheme"
  )
}
