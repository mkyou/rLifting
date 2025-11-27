#' Lifting Wavelet Transform (Forward)
#'
#' Performs the Forward Wavelet Transform using the Lifting Scheme.
#' Optimized with C++ backend.
#'
#' @param signal Numeric vector containing the input signal.
#' @param scheme A `lifting_scheme` object.
#' @param levels Integer. Number of decomposition levels.
#' @param extension Boundary extension mode: "symmetric" (default),
#' "periodic", or "zero".
#'
#' @return An object of class 'lwt' containing:
#' \item{coeffs}{List of details (d1..dn) and approximation (an).}
#' \item{scheme}{The scheme used.}
#' @export
#'
#' @examples
#' data = c(1, 2, 3, 4, 5, 6, 7, 8)
#' sch = lifting_scheme("haar")
#' res = lwt(data, sch, levels = 2)
#' print(res)
lwt = function(signal, scheme, levels = 1, extension = "symmetric") {

  if (!inherits(scheme, "lifting_scheme")) {
    stop("Invalid 'scheme' argument.")
  }
  n = length(signal)
  if (n < 2^levels) {
    stop("Signal is too short for the requested number of levels.")
  }

  final_len = n / (2^levels)
  if (final_len < 4) {
    warning(sprintf(
      "Residual signal at level %d has only %.1f samples.",
      levels, final_len
    ))
  }

  ext_int = switch(extension,
                   "symmetric" = 1L,
                   "periodic"  = 2L,
                   "zero"      = 3L,
                   1L)

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

#' Print method for LWT
#' @param x An object of class lwt.
#' @param ... Additional arguments.
#' @export
print.lwt = function(x, ...) {
  cat("--- LWT Decomposition (C++ Accelerated) ---\n")
  cat(sprintf("Levels: %d\n", x$levels))
  cat(sprintf("Wavelet: %s\n", x$scheme$wavelet))
  cat("Coefficients:\n")
  nms = names(x$coeffs)
  for (name in sort(nms)) {
    cat(sprintf("  %s: length %d\n", name, length(x$coeffs[[name]])))
  }
}
