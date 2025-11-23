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
  if (is.null(scheme)) scheme = lwt_obj$scheme

  levels = lwt_obj$levels
  coeffs = lwt_obj$coeffs
  ext_mode = if(!is.null(lwt_obj$extension)) {
    lwt_obj$extension
  } else {
    "symmetric"
  }

  current_app = coeffs[[paste0("a", levels)]]

  for (j in levels:1) {
    current_det = coeffs[[paste0("d", j)]]

    # 1. Inverter Normalizacao
    even = current_app / scheme$normalization[1]
    odd  = current_det / scheme$normalization[2]

    # 2. Inverter Passos de Lifting
    steps_rev = rev(scheme$steps)

    for (step in steps_rev) {
      if (step$type == "update") {
        update_val = .apply_filter_lifting(
          odd, step$coeffs,
          step$start_idx, ext_mode
        )

        # --- WARNING DE SEGURANCA ---
        if (length(even) != length(update_val)) {
          warning(sprintf(
            "Descasamento de tamanho no Nivel %d (Update). %d vs %d",
            j, length(even), length(update_val)
          ))
        }
        # ----------------------------

        len_min = min(length(even), length(update_val))
        even = even[1:len_min] - update_val[1:len_min]

      } else if (step$type == "predict") {
        pred_val = .apply_filter_lifting(
          even, step$coeffs,
          step$start_idx, ext_mode
        )

        # --- WARNING DE SEGURANCA ---
        if (length(odd) != length(pred_val)) {
          warning(sprintf(
            "Descasamento de tamanho no Nivel %d (Predict). %d vs %d",
            j, length(odd), length(pred_val)
          ))
        }
        # ----------------------------

        len_min = min(length(odd), length(pred_val))
        odd = odd[1:len_min] + pred_val[1:len_min]
      }
    }

    # 3. Merge
    n_total = length(even) + length(odd)
    merged = numeric(n_total)
    merged[seq(1, n_total, by = 2)] = even
    merged[seq(2, n_total, by = 2)] = odd

    current_app = merged
  }

  return(current_app[1:lwt_obj$original_len])
}
