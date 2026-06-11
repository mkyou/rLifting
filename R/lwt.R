#' Lifting Wavelet Transform (Forward)
#'
#' Performs the Forward Wavelet Transform using the Lifting Scheme.
#' Optimized with 'C++' backend.
#'
#' @param signal Numeric vector containing the input signal.
#' @param scheme A \code{lifting_scheme} object.
#' @param levels Integer. Number of decomposition levels.
#' @param extension Boundary extension mode: \code{"symmetric"} (default),
#'   \code{"periodic"}, \code{"zero"}, \code{"local_linear"} (linear
#'   extrapolation from boundary samples), or \code{"one_sided"} (asymmetric
#'   filter renormalisation at the boundary).
#' @param t Optional numeric vector of sample positions for irregular grids.
#'   Must be sorted and have the same length as \code{signal}. When supplied,
#'   irregular-grid Lagrange interpolation is applied in the predict steps
#'   and \code{lwt_obj$t} is stored for use by \code{ilwt()}. Ignored by
#'   \code{extension = "one_sided"} (with a warning).
#' @param ll_k Local-linear neighbourhood size, used only when
#'   \code{extension = "local_linear"}. Default 4L; minimum 2; clamped to the
#'   signal length if larger.
#'
#' @return An object of class \code{lwt}. It is a list containing
#'   \code{coeffs} (list of details d1..dn and approximation an),
#'   \code{scheme} (the scheme object used), \code{levels},
#'   \code{original_len}, \code{extension}, \code{ll_k}, and \code{t}.
#' @export
#'
#' @examples
#' data = c(1, 2, 3, 4, 5, 6, 7, 8)
#' sch = lifting_scheme("haar")
#' res = lwt(data, sch, levels = 2)
#' print(res)
lwt = function(signal, scheme, levels = 1, extension = "symmetric", t = NULL,
               ll_k = 4L) {

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
    "periodic" = 2L,
    "zero" = 3L,
    "local_linear" = 4L,
    "one_sided" = 5L,
    1L)

  if (extension == "local_linear" && ll_k > n)
    warning(sprintf(
      "ll_k (%d) > signal length (%d): clamped to n.", ll_k, n
    ))

  if (!is.null(t)) {
    if (length(t) != n) stop("'t' must have the same length as 'signal'.")
    if (is.unsorted(t))  stop("'t' must be sorted in increasing order.")
    if (extension == "one_sided")
      warning(paste0(
        "extension 'one_sided' ignores irregular grid positions: ",
        "Lagrange interpolation will not be applied. ",
        "Use 'symmetric' or 'local_linear' for irregular-grid processing."
      ))
    .check_irregular_scheme(scheme)
  }

  t_cpp = if (is.null(t)) numeric(0) else as.numeric(t)

  coeffs_list = lwt_cpp(
    as.numeric(signal),
    scheme$steps,
    as.numeric(scheme$normalization),
    as.integer(levels),
    as.integer(ext_int),
    t_cpp,
    as.integer(ll_k)
  )

  structure(
    list(
      coeffs = coeffs_list,
      scheme = scheme,
      levels = levels,
      original_len = n,
      extension = extension,
      ll_k = ll_k,
      t = t
    ),
    class = "lwt"
  )
}

#' Print method for LWT
#' @param x An object of class lwt.
#' @param ... Additional arguments.
#' @return Invisibly returns \code{NULL}. Called for side effects (printing).
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

#' Plot method for LWT Decomposition
#'
#' @param x An object of class \code{lwt}.
#' @param ... Additional arguments.
#' @return Invisibly returns \code{NULL}.
#' @export
plot.lwt = function(x, ...) {
  oldpar = par(no.readonly = TRUE)
  on.exit(par(oldpar))

  n_plots = x$levels + 1
  par(mfrow = c(n_plots, 1), mar = c(2, 4, 2, 1))

  for (i in 1:x$levels) {
    name = paste0("d", i)
    data = x$coeffs[[name]]
    ts.plot(data, main = paste("Detail Level", i), ylab = "Amp", col = "blue")
    grid()
  }

  approx_name = paste0("a", x$levels)
  ts.plot(x$coeffs[[approx_name]],
    main = paste("Approximation Level", x$levels),
    ylab = "Amp", col = "red")
  grid()
}
