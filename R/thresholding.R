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

#' SCAD Shrinkage (Antoniadis & Fan, 2001)
#'
#' Smoothly Clipped Absolute Deviation shrinkage. Three-region rule with
#' continuity at lambda, 2*lambda and a*lambda. Identity above a*lambda
#' eliminates the bias that soft thresholding imposes on large coefficients,
#' while keeping sparsity in the zero region.
#'
#' @param x Vector of coefficients.
#' @param lambda Positive threshold value.
#' @param a Shape parameter, a > 2. Default 3.7 (Fan-Li canonical).
#'
#' @return Processed vector.
#' @export
threshold_scad = function(x, lambda, a = 3.7) {
  threshold_scad_cpp(as.numeric(x), as.numeric(lambda), as.numeric(a))
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
#'
#' @return Numeric vector of the same length as \code{x} with thresholded coefficients.
#' @export
threshold = function(x, lambda, method = "soft", a = 3.7) {
  switch(
    method,
    hard = threshold_hard(x, lambda),
    soft = threshold_soft(x, lambda),
    semisoft = threshold_semisoft(x, lambda),
    scad = threshold_scad(x, lambda, a = a),
    stop("Unknown threshold method")
  )
}
