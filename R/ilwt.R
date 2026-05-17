#' Inverse Lifting Wavelet Transform ('C++' Accelerated)
#'
#' Reconstructs the original signal from wavelet coefficients.
#' Optimized with 'C++' backend.
#'
#' @param lwt_obj Object of class 'lwt' returned by `lwt()`.
#' @param scheme (Optional) `lifting_scheme` object.
#'  If NULL, uses the one from `lwt_obj`.
#'
#' @return Numeric vector containing the reconstructed signal.
#' @export
#'
#' @examples
#' s = c(1, 2, 3, 4)
#' sch = lifting_scheme("haar")
#' fwd = lwt(s, sch)
#' rec = ilwt(fwd)
#' print(rec) # Should match s
ilwt = function(lwt_obj, scheme = NULL) {

  if (!inherits(lwt_obj, "lwt")) stop("Input must be an 'lwt' object.")
  if (is.null(scheme)) scheme = lwt_obj$scheme

  ext_mode = if (!is.null(lwt_obj$extension)) lwt_obj$extension else "symmetric"
  ext_int = switch(
    ext_mode,
    "symmetric"    = 1L,
    "periodic"     = 2L,
    "zero"         = 3L,
    "local_linear" = 4L,
    "one_sided"    = 5L,
    1L
  )

  t_cpp  = if (!is.null(lwt_obj$t)) as.numeric(lwt_obj$t) else numeric(0)
  ll_k   = if (!is.null(lwt_obj$ll_k)) as.integer(lwt_obj$ll_k) else 2L

  res = ilwt_cpp(
    lwt_obj$coeffs,
    scheme$steps,
    as.numeric(scheme$normalization),
    as.integer(lwt_obj$levels),
    as.integer(ext_int),
    as.integer(lwt_obj$original_len),
    t_cpp,
    ll_k
  )

  return(res)
}
