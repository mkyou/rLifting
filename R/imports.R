#' rLifting: Real-Time Signal Denoising via Lifting Scheme
#'
#' A package for wavelet transforms using the Lifting Scheme, focused on
#' real-time (causal) applications, adaptive thresholding, and custom
#' wavelet construction.
#'
#' @docType package
#' @name rLifting
#' @importFrom stats median rnorm cor sd
#' @importFrom utils setTxtProgressBar tail txtProgressBar head
#' @importFrom graphics par ts.plot grid lines legend title
#' @importFrom Rcpp sourceCpp
#' @useDynLib rLifting, .registration = TRUE
NULL
