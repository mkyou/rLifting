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
#' @param threshold_method Threshold-selection rule. Currently only
#'   `"universal"` (Donoho-Johnstone universal threshold with the recursive
#'   per-level decay parameterised by \code{alpha} and \code{beta}).
#' @param shrinkage Shrinkage rule applied above the threshold:
#'   `"hard"`, `"soft"`, or `"semisoft"`.
#' @param method Deprecated. Use \code{shrinkage} instead. If provided, takes
#'   precedence over \code{shrinkage} with a deprecation warning.
#' @param extension Extension mode ("symmetric", "periodic", "zero",
#'   "local_linear").
#'
#' @return Filtered numeric vector (same length as input).
#' @export
denoise_signal_offline = function(
  signal,
  scheme,
  alpha = 0.3,
  beta = 1.2,
  levels = 3,
  threshold_method = "universal",
  shrinkage = NULL,
  method = NULL,
  extension = "symmetric",
  t = NULL,
  ll_k = 4L
) {

  resolved = .resolve_shrinkage_args(method, shrinkage, threshold_method)

  if (extension == "local_linear" && ll_k > length(signal))
    warning(sprintf(
      "ll_k (%d) > signal length (%d): clamped to n.", ll_k, length(signal)
    ))

  if (!is.null(t)) {
    if (length(t) != length(signal))
      stop("'t' must have the same length as 'signal'.")
    if (is.unsorted(t))
      stop("'t' must be sorted in increasing order.")
    if (extension == "one_sided")
      warning(paste0(
        "extension 'one_sided' ignores irregular grid positions: ",
        "Lagrange interpolation will not be applied. ",
        "Use 'symmetric' or 'local_linear' for irregular-grid processing."
      ))
    .check_irregular_scheme(scheme)
  }

  ext_int = switch(
    extension,
    "symmetric" = 1L,
    "periodic" = 2L,
    "zero" = 3L,
    "local_linear" = 4L,
    "one_sided" = 5L,
    1L
  )

  t_cpp = if (is.null(t)) numeric(0) else as.numeric(t)

  res = denoise_offline_cpp(
    as.numeric(signal),
    scheme$steps,
    as.numeric(scheme$normalization),
    as.integer(levels),
    as.numeric(alpha),
    as.numeric(beta),
    as.character(resolved$shrinkage),
    as.integer(ext_int),
    t_cpp,
    as.integer(ll_k),
    as.character(resolved$threshold_method)
  )

  return(res)
}
