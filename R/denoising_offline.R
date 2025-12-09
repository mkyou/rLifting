#' Offline Denoising (Global Batch)
#'
#' Performs denoising on the entire signal at once using a non-causal approach.
#' Uses global statistics for recursive threshold calculation (Eq. 9).
#' This function is fully optimized in 'C++' (Zero-Allocation).
#'
#' @param signal Numeric vector containing the complete signal.
#' @param scheme A \code{lifting_scheme} object.
#' @param alpha Recursive threshold parameter.
#' @param beta Threshold scale factor.
#' @param levels Number of decomposition levels.
#' @param method Thresholding method ("hard", "soft", "semisoft").
#' @param extension Extension mode ("symmetric", "periodic", "zero").
#'
#' @return Filtered numeric vector (same length as input).
#' @export
denoise_signal_offline = function(
    signal,
    scheme,
    alpha = 0.3,
    beta = 1.2,
    levels = 3,
    method = "semisoft",
    extension = "symmetric"
) {

  ext_int = switch(
    extension,
    "symmetric" = 1L,
    "periodic"  = 2L,
    "zero"      = 3L,
    1L
  )

  res = denoise_offline_cpp(
    as.numeric(signal),
    scheme$steps,
    as.numeric(scheme$normalization),
    as.integer(levels),
    as.numeric(alpha),
    as.numeric(beta),
    as.character(method),
    as.integer(ext_int)
  )

  return(res)
}
