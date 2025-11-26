#' Calculate Adaptive Threshold (Recursive Thresholding)
#'
#' Estimates the noise threshold based on current window statistics.
#' Accelerated with C++.
#'
#' @param lwt_obj Object returned by `lwt()` function.
#' @param alpha Recursive adjustment parameter (Eq. 9).
#' @param beta Initial threshold scale factor (Eq. 9).
#'
#' @return List of thresholds for each decomposition level (d1, d2...).
#' @export
compute_adaptive_threshold = function(lwt_obj, alpha = 0.3, beta = 1.2) {

  d1 = lwt_obj$coeffs$d1
  if (length(d1) == 0) return(list(d1 = 0))

  # Conta quantos niveis de detalhe existem
  det_names = grep("^d[0-9]+", names(lwt_obj$coeffs), value = TRUE)
  max_level = length(det_names)

  # Chama C++
  lambdas_vec = compute_thresholds_cpp(d1, max_level, alpha, beta)

  # Converte de volta para lista nomeada (para compatibilidade)
  lambdas_list = list()
  for(i in 1:max_level) {
    lambdas_list[[paste0("d", i)]] = lambdas_vec[i]
  }

  return(lambdas_list)
}
