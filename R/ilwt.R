#' Inverse Lifting Wavelet Transform
#'
#' Reconstroi o sinal original a partir dos coeficientes de ondaleta.
#'
#' @param lwt_obj Objeto da classe 'lwt' gerado pela funcao lwt().
#' @param scheme (Opcional) Objeto lifting_scheme. Se NULL, usa o do lwt_obj.
#'
#' @return Vetor numerico com o sinal reconstruido.
#' @export
#'
#' @examples
#' s = c(1, 2, 3, 4)
#' sch = lifting_scheme("haar")
#' fwd = lwt(s, sch)
#' rec = ilwt(fwd)
#' print(rec) # Deve ser igual a s
ilwt = function(lwt_obj, scheme = NULL) {

  if (!inherits(lwt_obj, "lwt")) stop("Entrada deve ser um objeto 'lwt'.")

  # Se o usuario nao fornecer esquema, usa o que esta gravado no objeto
  if (is.null(scheme)) scheme = lwt_obj$scheme

  levels = lwt_obj$levels
  coeffs = lwt_obj$coeffs

  # Comecamos pela aproximacao mais grosseira (ultimo nivel)
  current_app = coeffs[[paste0("a", levels)]]

  # Loop reverso: do nivel mais profundo (levels) ate 1
  for (j in levels:1) {

    # Recupera os detalhes deste nivel
    current_det = coeffs[[paste0("d", j)]]

    # 1. Inverter Normalizacao
    # Na ida: a = even * K1, d = odd * K2
    # Na volta: even = a / K1, odd = d / K2
    even = current_app / scheme$normalization[1]
    odd  = current_det / scheme$normalization[2]

    # 2. Inverter Passos de Lifting (Ordem Reversa)
    steps_rev = rev(scheme$steps)

    for (step in steps_rev) {
      # Importante: O input da filtragem (.apply_filter) nao muda de sinal,
      # apenas a operacao de atualizacao (soma/subtracao).

      if (step$type == "update") {
        # Ida: a = even + U(d)
        # Volta: even = a - U(d)
        # Nota: aqui 'odd' representa o 'd' da formula
        update_val = .apply_filter_lifting(odd, step$coeffs)

        # Ajuste de tamanho (robustez)
        len_min = min(length(even), length(update_val))
        even = even[1:len_min] - update_val[1:len_min]

      } else if (step$type == "predict") {
        # Ida: d = odd - P(even)
        # Volta: odd = d + P(even)
        pred_val = .apply_filter_lifting(even, step$coeffs)

        len_min = min(length(odd), length(pred_val))
        odd = odd[1:len_min] + pred_val[1:len_min]
      }
    }

    # 3. Inverse Lazy Wavelet (Merge)
    # Intercalar even e odd
    # even -> indices impares do R (1, 3, 5...)
    # odd  -> indices pares do R (2, 4, 6...)

    n_total = length(even) + length(odd)
    merged = numeric(n_total)

    merged[seq(1, n_total, by = 2)] = even
    merged[seq(2, n_total, by = 2)] = odd

    # Atualiza a aproximacao para o proximo nivel (acima)
    current_app = merged
  }

  # Retorna o sinal reconstruido
  # (Opcional: cortar para lwt_obj$original_len se houver padding extra)
  return(current_app[1:lwt_obj$original_len])
}
