#' Inverse Lifting Wavelet Transform (C++ Accelerated)
#'
#' Reconstroi o sinal original a partir dos coeficientes de ondaleta.
#' Otimizado com backend C++.
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

  # Mapeamento de extensao
  ext_mode = if(!is.null(lwt_obj$extension)) lwt_obj$extension else "symmetric"
  ext_int = switch(
    ext_mode,
    "symmetric" = 1L,
    "periodic" = 2L,
    "zero" = 3L,
    1L
    )

  # Chamada ao Core C++
  # Passamos a lista de coeficientes inteira
  res = ilwt_cpp(
    lwt_obj$coeffs,
    scheme$steps,
    as.numeric(scheme$normalization),
    as.integer(lwt_obj$levels),
    as.integer(ext_int),
    as.integer(lwt_obj$original_len)
  )

  return(res)
}
