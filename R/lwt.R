#' Lifting Wavelet Transform (Forward)
#'
#' Performs the Forward Wavelet Transform using the Lifting Scheme.
#' Optimized with 'C++' backend.
#'
#' @param signal Numeric vector containing the input signal.
#' @param scheme A \code{lifting_scheme} object.
#' @param levels Integer. Number of decomposition levels.
#' @param extension Boundary extension mode: "symmetric" (default),
#' "periodic", or "zero".
#'
#' @return An object of class 'lwt'. It is a list containing
#' 'coeffs' (list of details d1..dn and approximation an) and
#' 'scheme' (the scheme object used).
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
  # Setup layout: 1 row per level (+1 for approximation)
  oldpar = par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  n_plots = x$levels + 1
  par(mfrow = c(n_plots, 1), mar = c(2, 4, 2, 1))
  
  # Plot Details (d1 to dn)
  for (i in 1:x$levels) {
    name = paste0("d", i)
    data = x$coeffs[[name]]
    ts.plot(data, main = paste("Detail Level", i), ylab = "Amp", col = "blue")
    grid()
  }
  
  # Plot Approximation (an)
  approx_name = paste0("a", x$levels)
  ts.plot(x$coeffs[[approx_name]], main = paste("Approximation Level", x$levels), 
          ylab = "Amp", col = "red")
  grid()
}

