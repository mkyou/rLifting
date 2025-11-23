#' Lifting Wavelet Transform (Forward)
#'
#' Executa a transformada wavelet via esquema de lifting.
#' Suporta decomposicao multinivel.
#'
#' @param signal Vetor numerico contendo o sinal.
#' @param scheme Objeto da classe 'lifting_scheme'.
#' @param levels Inteiro. Numero de niveis de decomposicao.
#' @param extension Modo de tratamento de borda: "symmetric" (padrao),
#' "periodic", "zero".
#'
#' @return Objeto S3 da classe 'lwt'. Contem:
#' \item{coeffs}{Lista com detalhes (d1..dn) e aproximacao (an).}
#' \item{scheme}{O esquema utilizado.}
#' @export
#'
#' @examples
#' data = c(1, 2, 3, 4, 5, 6, 7, 8)
#' sch = lifting_scheme("haar")
#' res = lwt(data, sch, levels = 2)
#' print(res)
lwt = function(signal, scheme, levels = 1, extension = "symmetric") {

  # Validacoes
  if (!inherits(scheme, "lifting_scheme")) {
    stop("Argumento 'scheme' invalido.")
  }
  n = length(signal)
  if (n < 2^levels) {
    stop("Sinal muito curto para o numero de niveis solicitado.")
  }

  # --- WARNING DE SUPORTE ---
  # Se o sinal no ultimo nivel for muito curto (< 4 amostras),
  # filtros complexos (CDF9/7, DB4) vao sofrer distorcao severa.
  final_len = n / (2^levels)
  if (final_len < 4) {
    warning(sprintf(
      "Sinal residual no nivel %d tem apenas %.1f amostras. %s",
      levels, final_len, "Artefatos de borda podem dominar."
    ))
  }

  # Inicializacao
  current_app = signal
  details_list = list()
  boundary_mode = match.arg(extension, c("symmetric", "periodic", "zero"))

  for (j in 1:levels) {
    n_curr = length(current_app)

    # 1. Lazy Wavelet (Split)
    # No R (indice 1-based):
    # Indices impares (1, 3...) sao as amostras PARES originais (0, 2...)
    # Indices pares (2, 4...) sao as amostras IMPARES originais (1, 3...)
    even = current_app[seq(1, n_curr, by = 2)]
    odd  = current_app[seq(2, n_curr, by = 2)]

    # 2. Executa os passos do esquema (Predict/Update)
    for (step in scheme$steps) {
      if (step$type == "predict") {
        # d = odd - P(even)
        # step$coeffs contem os coeficientes de P
        prediction = .apply_filter_lifting(
          even, step$coeffs,
          step$start_idx, extension
        )

        # Garante dimensoes iguais (pode ocorrer diferenca de 1 em odd/even)
        len_min = min(length(odd), length(prediction))
        odd = odd[1:len_min] - prediction[1:len_min]

      } else if (step$type == "update") {
        # a = even + U(d)
        # step$coeffs contem os coeficientes de U
        update_val = .apply_filter_lifting(
          odd, step$coeffs,
          step$start_idx, extension
        )

        len_min = min(length(even), length(update_val))
        even = even[1:len_min] + update_val[1:len_min]
      }
    }

    # 3. Normalizacao
    # scheme$normalization = c(fator_aprox, fator_detalhe)
    approx_c = even * scheme$normalization[1]
    detail_c = odd  * scheme$normalization[2]

    # Armazena detalhe e prepara aproximacao para proximo nivel
    details_list[[paste0("d", j)]] = detail_c
    current_app = approx_c
  }

  # Adiciona a ultima aproximacao na lista de coeficientes
  # A estrutura usual e: d1, d2, ..., an
  coeffs_out = details_list
  coeffs_out[[paste0("a", levels)]] = current_app

  structure(
    list(
      coeffs = coeffs_out,
      scheme = scheme,
      levels = levels,
      original_len = n,
      extension = boundary_mode # Guardamos o modo usado
    ),
    class = "lwt"
  )
}

#' Print method para LWT
#' @param x Objeto da classe lwt.
#' @param ... Argumentos adicionais (nao utilizados).
#' @export
print.lwt = function(x, ...) {
  cat("--- LWT Decomposition ---\n")
  cat(sprintf("Levels: %d\n", x$levels))
  cat(sprintf("Wavelet: %s\n", x$scheme$wavelet))
  cat("Coefficients:\n")
  for (name in names(x$coeffs)) {
    cat(sprintf("  %s: length %d\n", name, length(x$coeffs[[name]])))
  }
}
