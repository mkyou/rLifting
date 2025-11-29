#' Hard Thresholding
#'
#' Sets coefficients below the threshold to zero, keeping others unchanged.
#' Known as the "keep or kill" policy.
#'
#' @param x Vector of coefficients (details).
#' @param lambda Positive threshold value.
#'
#' @return Processed vector.
#' @export
threshold_hard = function(x, lambda) {
  threshold_hard_cpp(as.numeric(x), as.numeric(lambda))
}

#' Soft Thresholding
#'
#' Sets coefficients below the threshold to zero and shrinks others
#' towards zero.
#' Reduces noise variance but introduces amplitude bias.
#'
#' @param x Vector of coefficients.
#' @param lambda Positive threshold value.
#'
#' @return Processed vector.
#' @export
threshold_soft = function(x, lambda) {
  threshold_soft_cpp(as.numeric(x), as.numeric(lambda))
}

#' Semisoft Shrinkage (Hyperbolic)
#'
#' Implementation based on Liu et al. (2014).
#' Combines the stability of Soft Thresholding with the amplitude
#' precision of Hard Thresholding.
#' Function: \code{sign(x) * sqrt(x^2 - lambda^2)} for values above lambda.
#'
#' @param x Vector of coefficients.
#' @param lambda Positive threshold value.
#'
#' @return Processed vector.
#' @export
threshold_semisoft = function(x, lambda) {
  threshold_semisoft_cpp(as.numeric(x), as.numeric(lambda))
}

#' General Thresholding Wrapper
#'
#' @param x Input vector.
#' @param lambda Threshold value.
#' @param method Method: "hard", "soft" or "semisoft".
#' @export
threshold = function(x, lambda, method = "soft") {
  switch(
    method,
    hard = threshold_hard(x, lambda),
    soft = threshold_soft(x, lambda),
    semisoft = threshold_semisoft(x, lambda),
    stop("Unknown threshold method")
  )
}
