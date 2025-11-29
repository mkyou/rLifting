#' Lifting Scheme Constructor
#'
#' Creates an S3 object containing the prediction (P) and update (U) steps
#' required for the Lifting Transform.
#'
#' @param wavelet Wavelet name (string). Options:
#'  "haar", "db2", "cdf53", "cdf97", "dd4", "lazy".
#' @param custom_steps List of custom steps (optional).
#'  If provided, ignores internal lookup.
#' @param custom_norm Normalization vector (optional).
#'
#' @return An object of class \code{lifting_scheme}.
#' @export
lifting_scheme = function(
    wavelet = "haar",
    custom_steps = NULL,
    custom_norm = NULL
) {

  steps = list()
  norm_factors = c(1, 1)

  if (!is.null(custom_steps)) {
    steps = custom_steps
    if (!is.null(custom_norm)) norm_factors = custom_norm
  } else {
    config = .get_wavelet_config(wavelet)
    steps = config$steps
    norm_factors = config$norm
  }

  structure(
    list(wavelet = wavelet, steps = steps, normalization = norm_factors),
    class = "lifting_scheme"
  )
}

#' Implements factorizations based on Daubechies and Sweldens (1998).
#' @keywords internal
.get_wavelet_config = function(name) {
  if (name == "lazy") {
    return(list(steps = list(), norm = c(1, 1)))
  }

  if (name == "haar") {
    steps = list(
      list(type = "predict", coeffs = c(1), start_idx = 0),
      list(type = "update",  coeffs = c(0.5), start_idx = 0)
    )
    norm = c(sqrt(2), 1/sqrt(2))
    return(list(steps = steps, norm = norm))
  }

  if (name == "db2") {
    sqrt3 = sqrt(3)
    steps = list(
      list(type = "predict", coeffs = c(sqrt3), start_idx = 0),
      list(
        type = "update",
        coeffs = c(sqrt3/4, (sqrt3 - 2)/4), start_idx = 0
      ),
      list(type = "predict", coeffs = c(-1), start_idx = -1)
    )
    norm = c((sqrt3 + 1) / sqrt(2), (sqrt3 - 1) / sqrt(2))
    return(list(steps = steps, norm = norm))
  }

  if (name == "cdf53" || name == "bior2.2") {
    steps = list(
      list(type = "predict", coeffs = c(0.5, 0.5), start_idx = 0),
      list(type = "update", coeffs = c(0.25, 0.25), start_idx = -1)
    )
    norm = c(sqrt(2), 1/sqrt(2))
    return(list(steps = steps, norm = norm))
  }

  if (name == "cdf97" || name == "bior4.4") {
    alpha = -1.586134342
    beta  = -0.05298011854
    gamma = 0.8829110762
    delta = 0.4435068522
    zeta  = 1.149604398

    steps = list(
      list(type = "predict", coeffs = c(-alpha, -alpha), start_idx = 0),
      list(type = "update",  coeffs = c(beta, beta), start_idx = -1),
      list(type = "predict", coeffs = c(-gamma, -gamma), start_idx = 0),
      list(type = "update",  coeffs = c(delta, delta), start_idx = -1)
    )
    norm = c(zeta, 1/zeta)
    return(list(steps = steps, norm = norm))
  }

  if (name == "dd4" || name == "interp4") {
    p_coeffs = c(-1/16, 9/16, 9/16, -1/16)
    u_coeffs = p_coeffs / 2

    steps = list(
      list(type = "predict", coeffs = p_coeffs, start_idx = -1),
      list(type = "update",  coeffs = u_coeffs, start_idx = -1)
    )
    norm = c(sqrt(2), 1/sqrt(2))
    return(list(steps = steps, norm = norm))
  }

  stop(paste("Wavelet", name, "not supported."))
}

#' Print method
#' @param x object of class \code{lifting_scheme}.
#' @param ... additional arguments.
#' @export
print.lifting_scheme = function(x, ...) {
  cat(sprintf("Lifting Scheme: %s\n", x$wavelet))
  cat(sprintf("Steps: %d\n", length(x$steps)))
  cat(sprintf("Norm:  Approx=%.4f, Detail=%.4f\n",
              x$normalization[1], x$normalization[2]))
}


#' Create an individual Lifting Step
#'
#' Helper function to create prediction (P) or update (U) steps,
#' abstracting the complexity of index management.
#'
#' @param type Step type: "predict" (P) or "update" (U).
#' @param coeffs Numeric vector containing the filter coefficients.
#' @param start_idx (Optional) Manual start index. If provided, ignores the
#'        \code{position} parameter. Use this for fine-grained control.
#' @param position Automatic index adjustment
#'  (used only if \code{start_idx} is NULL):
#' \itemize{
#'   \item "center": Centers the filter (default).
#'   \item "left": Causal filter (looks into the past).
#'   \item "right": Anti-causal filter (looks into the future).
#' }
#'
#' @return A list formatted for the internal lifting engine.
#' @export
lift_step = function(
    type = c("predict", "update"),
    coeffs,
    start_idx = NULL,
    position = "center"
) {

  type = match.arg(type)
  n = length(coeffs)

  if (is.null(start_idx)) {
    if (position == "center") {
      start_idx = -floor((n - 1) / 2)
    } else if (position == "left") {
      start_idx = -n + 1
    } else {
      start_idx = 0
    }
  }

  list(type = type, coeffs = coeffs, start_idx = start_idx)
}

#' Create a custom wavelet
#'
#' Wrapper to create a \code{lifting_scheme} object from manual steps.
#'
#' @param name Identifier name for the wavelet.
#' @param steps List of steps created via \code{lift_step}.
#' @param norm Normalization vector c(K, 1/K).
#'
#' @return An object of class \code{lifting_scheme}.
#' @export
#' @examples
#' p1 = lift_step("predict", c(1), position="center")
#' u1 = lift_step("update", c(0.5), position="center")
#' w = custom_wavelet("HaarManual", list(p1, u1), c(1.41, 0.707))
custom_wavelet = function(name, steps, norm = c(1, 1)) {
  structure(
    list(wavelet = name, steps = steps, normalization = norm),
    class = "lifting_scheme"
  )
}
