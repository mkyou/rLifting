#' Calculate Adaptive Threshold (Recursive)
#'
#' Estimates the optimal noise threshold based on current window statistics.
#' Implements the recursive formula from Liu et al. (2014).
#' Accelerated with 'C++'.
#'
#' @param lwt_obj Object returned by \code{lwt()}.
#' @param alpha Recursive adjustment parameter (Eq. 9).
#' @param beta Initial threshold scale factor (Eq. 9).
#'
#' @return Object of class \code{adaptive_thresholds} (a list of thresholds).
#' @export
compute_adaptive_threshold = function(lwt_obj, alpha = 0.3, beta = 1.2) {

  d1 = lwt_obj$coeffs$d1
  if (length(d1) == 0) {
    res = list(d1 = 0)
    class(res) = "adaptive_thresholds"
    return(res)
  }

  det_names = grep("^d[0-9]+", names(lwt_obj$coeffs), value = TRUE)
  max_level = length(det_names)

  lambdas_vec = compute_thresholds_cpp(d1, max_level, alpha, beta)

  lambdas_list = list()
  for (i in 1:max_level) {
    lambdas_list[[paste0("d", i)]] = lambdas_vec[i]
  }

  class(lambdas_list) = "adaptive_thresholds"
  return(lambdas_list)
}

#' Print method for Adaptive Thresholds
#'
#' @param x Object of class \code{adaptive_thresholds}.
#' @param ... Additional arguments.
#' @return Invisibly returns \code{x}.
#' @export
print.adaptive_thresholds = function(x, ...) {
  cat("--- Adaptive Thresholds (Recursive) ---\n")
  for (name in names(x)) {
    cat(sprintf("  %s: %.6f\n", name, x[[name]]))
  }
  invisible(x)
}

#' Plot method for Adaptive Thresholds
#'
#' @param x Object of class \code{adaptive_thresholds}.
#' @param ... Additional arguments.
#' @return Invisibly returns \code{NULL}.
#' @export
plot.adaptive_thresholds = function(x, ...) {
  vals = unlist(x)
  barplot(vals, main = "Adaptive Thresholds per Level",
          ylab = "Threshold Value", col = "steelblue", border = NA, ...)
  grid(nx = NA, ny = NULL)
}
