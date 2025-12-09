#' rLifting: High-Performance Wavelet Lifting Transforms
#'
#' A unified framework for Wavelet Transforms using the Lifting Scheme.
#' It provides robust tools for offline signal analysis and functional data
#' analysis (FDA), while also enabling high-performance causal processing
#' for real-time applications via a specialized 'C++' core.
#'
#' @name rLifting
#' @importFrom stats median rnorm cor sd ts.plot
#' @importFrom utils setTxtProgressBar tail txtProgressBar head
#' @importFrom graphics par grid lines legend title
#' @importFrom Rcpp sourceCpp
#' @useDynLib rLifting, .registration = TRUE
"_PACKAGE"
