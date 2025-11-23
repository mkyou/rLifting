#' Hard Thresholding
#'
#' Zera coeficientes abaixo do limiar, mantem os demais inalterados.
#' Equivalente a "keep or kill".
#'
#' @param x Vetor de coeficientes (detalhes).
#' @param lambda Limiar (threshold) positivo.
#'
#' @return Vetor processado.
#' @export
threshold_hard = function(x, lambda) {
  # Formula: x se |x| >= lambda, senao 0
  # Implementacao vetorizada rapida (masking)
  x * (abs(x) >= lambda)
}

#' Soft Thresholding
#'
#' Zera coeficientes abaixo do limiar e contrai os demais em direcao a zero.
#' Reduz o ruido mas introduz bias (vies) na amplitude.
#'
#' @param x Vetor de coeficientes.
#' @param lambda Limiar positivo.
#'
#' @return Vetor processado.
#' @export
threshold_soft = function(x, lambda) {
  # Formula: sign(x) * max(0, |x| - lambda)
  # Em R, pmax e muito rapido
  sign(x) * pmax(0, abs(x) - lambda)
}

#' Semisoft Shrinkage (Hyperbolic)
#'
#' Implementacao baseada em Liu et al. (2014).
#' Combina a estabilidade do soft com a precisao de amplitude do Hard.
#' Funcao: sign(x) * sqrt(x^2 - lambda^2) para |x| > lambda.
#'
#' @param x Vetor de coeficientes.
#' @param lambda Limiar positivo.
#'
#' @return Vetor processado.
#' @export
threshold_semisoft = function(x, lambda) {
  # Formula Eq(12) do artigo:
  # Se |x| >= lambda -> sign(x) * sqrt(x^2 - lambda^2)
  # Senao -> 0

  abs_x = abs(x)
  mask = abs_x >= lambda

  y = numeric(length(x))

  # Calculamos apenas onde a mascara e TRUE para evitar sqrt de negativo
  if (any(mask)) {
    vals = x[mask]
    y[mask] = sign(vals) * sqrt(vals^2 - lambda^2)
  }

  return(y)
}

#' Wrapper Geral de Thresholding
#'
#' @param x Vetor de entrada.
#' @param lambda Limiar.
#' @param method "hard", "soft" ou "semisoft".
#' @export
threshold = function(x, lambda, method = "soft") {
  switch(method,
         hard = threshold_hard(x, lambda),
         soft = threshold_soft(x, lambda),
         semisoft = threshold_semisoft(x, lambda),
         stop("Metodo de threshold desconhecido")
  )
}
