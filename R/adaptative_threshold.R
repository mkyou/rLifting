#' Calcula limiar adaptativo (Recursive Thresholding)
#'
#' Estima o limiar de ruido baseado na estatistica da janela atual,
#' seguindo o metodo recursivo de Liu et al. 2014 (Eq. 9).
#'
#' @param lwt_obj Objeto retornado pela funcao `lwt()`.
#' @param alpha Parametro de ajuste recursivo (Eq. 9).
#' @param beta Fator de escala do limiar inicial (Eq. 9).
#'
#' @return Lista de limiares para cada nivel de decomposicao (d1, d2...).
#' @export
compute_adaptive_threshold = function(lwt_obj, alpha = 0.3, beta = 1.2) {

  # 1. Estimar o ruido base (Sigma) usando o Nivel 1 (Detalhes mais finos)
  # Ref: Eq (7) e Eq (9)
  # sigma = median(|d1|) / 0.6745
  d1 = lwt_obj$coeffs$d1

  # Protecao contra vetor vazio ou de zeros
  if (length(d1) == 0) return(list(d1 = 0))

  mad_val = median(abs(d1))
  sigma = mad_val / 0.6745

  # Se o sinal for constante, sigma e zero -> lambda e zero (sem denoising)
  if (sigma < 1e-15) {
    return(
      lapply(
        lwt_obj$coeffs[grep("^d", names(lwt_obj$coeffs))], function(x) 0
        )
      )
  }

  # 2. Calcular Lambda do Nivel 1
  # N1 = Tamanho da amostra no nivel 1 (Comprimento de d1)
  # Ref Eq (9): lambda_1 = beta * sigma * sqrt(2 * log(N1))
  N1 = length(d1)
  lambda_1 = beta * sigma * sqrt(2 * log(N1))

  lambdas = list()
  lambdas[["d1"]] = lambda_1

  # 3. Calcular Lambdas Recursivos para Niveis Superiores (i > 1)
  # Ref Eq (9): lambda_i = lambda_{i-1} * (i-1) / (i + alpha - 1)

  # Identifica quantos niveis existem (d1, d2, d3...)
  det_names = grep("^d[0-9]+", names(lwt_obj$coeffs), value = TRUE)
  max_level = length(det_names)

  if (max_level > 1) {
    for (i in 2:max_level) {
      prev_lambda = lambdas[[paste0("d", i-1)]]

      # Formula recursiva
      # Nota: O artigo usa 'i' como indice do nivel atual
      factor = (i - 1) / (i + alpha - 1)
      curr_lambda = prev_lambda * factor

      lambdas[[paste0("d", i)]] = curr_lambda
    }
  }

  return(lambdas)
}
