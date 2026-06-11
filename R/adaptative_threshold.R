#' Calculate Adaptive Threshold (Universal / Recursive)
#'
#' Estimates the per-level noise threshold from the finest-level detail
#' coefficients and applies the recursive Liu et al. (2014) decay across
#' levels. This is the step-by-step entry point for the universal threshold
#' rule.
#'
#' To use SureShrink instead, call \code{denoise_signal_offline()} or the
#' causal/stream functions with \code{threshold_method = "sure"} — the SURE
#' branch lives in the C++ engine and is not exposed as a standalone R
#' routine. For automatic selection of \code{alpha} and \code{beta},
#' see \code{\link{tune_alpha_beta}}.
#'
#' @param lwt_obj Object returned by \code{lwt()}.
#' @param alpha Recursive adjustment parameter (Eq. 9 of Liu et al., 2014).
#' @param beta Initial threshold scale factor (Eq. 9 of Liu et al., 2014).
#'
#' @return Object of class \code{adaptive_thresholds} (a list of thresholds).
#'
#' @references
#' Donoho, D. L., & Johnstone, I. M. (1994). Ideal spatial adaptation by
#' wavelet shrinkage. \emph{Biometrika}, 81(3), 425--455.
#'
#' Liu, Z., Mi, Y., & Mao, Y. (2014). Improved real-time denoising method
#' based on lifting wavelet transform. \emph{Measurement Science Review},
#' 14(3), 152--159. \doi{10.2478/msr-2014-0020}
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
