#' Lifting Wavelet Transform (Forward)
#'
#' Executa a transformada wavelet via esquema de lifting.
#' Otimizado com backend C++.
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

  # Validacoes Basicas
  if (!inherits(scheme, "lifting_scheme")) {
    stop("Argumento 'scheme' invalido.")
  }
  n = length(signal)
  if (n < 2^levels) {
    stop("Sinal muito curto para o numero de niveis solicitado.")
  }

  # Warning de sinal curto (mantido para consistencia)
  final_len = n / (2^levels)
  if (final_len < 4) {
    warning(sprintf(
      "Sinal residual no nivel %d tem apenas %.1f amostras.",
      levels, final_len
    ))
  }

  # Mapeia extensao para inteiro (C++)
  # 1=symmetric, 2=periodic, 3=zero
  ext_int = switch(extension,
                   "symmetric" = 1L,
                   "periodic"  = 2L,
                   "zero"      = 3L,
                   1L)

  # --- CHAMADA AO CORE C++ ---
  # Passamos a lista de steps e os parametros brutos.
  coeffs_list = lwt_cpp(
    as.numeric(signal),
    scheme$steps,
    as.numeric(scheme$normalization),
    as.integer(levels),
    as.integer(ext_int)
  )

  structure(
    list(
      coeffs = coeffs_list,
      scheme = scheme,
      levels = levels,
      original_len = n,
      extension = extension
    ),
    class = "lwt"
  )
}

#' Print method para LWT
#' @param x Objeto da classe lwt.
#' @param ... Argumentos adicionais (nao utilizados).
#' @export
print.lwt = function(x, ...) {
  cat("--- LWT Decomposition (C++ Accelerated) ---\n")
  cat(sprintf("Levels: %d\n", x$levels))
  cat(sprintf("Wavelet: %s\n", x$scheme$wavelet))
  cat("Coefficients:\n")
  # A ordem da lista que vem do C++
  # pode nao estar ordenada bonitinha (Hash map),
  # entao ordenamos para printar (d1, d2... an)
  nms = names(x$coeffs)
  # Logica simples de ordenacao visual
  for (name in sort(nms)) {
    cat(sprintf("  %s: length %d\n", name, length(x$coeffs[[name]])))
  }
}
