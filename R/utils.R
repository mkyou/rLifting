#' Signal padding/extension handling
#'
#' Kept for compatibility with legacy tests and visualization,
#' but the C++ core now performs virtual padding on-the-fly.
#'
#' @param x Input signal.
#' @param n_left Samples to the left.
#' @param n_right Samples to the right.
#' @param mode Extension mode: "symmetric", "periodic", "zero".
#'
#' @keywords internal
.pad_signal = function(x, n_left, n_right, mode = "symmetric") {
  n = length(x)
  if (n == 0) return(x)

  left_pad = numeric(0)
  right_pad = numeric(0)

  if (mode == "zero") {
    if (n_left > 0) left_pad = rep(0, n_left)
    if (n_right > 0) right_pad = rep(0, n_right)
  } else if (mode == "periodic") {
    if (n_left > 0) {
      full_idx = rep(1:n, ceiling(n_left/n) + 1)
      left_pad = x[utils::tail(full_idx, n_left)]
    }
    if (n_right > 0) {
      idx = (0:(n_right-1)) %% n + 1
      right_pad = x[idx]
    }
  } else if (mode == "symmetric") {
    if (n_left > 0) {
      nl = min(n_left, n)
      left_pad = x[seq(nl, 1, by = -1)]
      if (n_left > n) left_pad = c(rep(x[n], n_left - n), left_pad)
    }
    if (n_right > 0) {
      nr = min(n_right, n)
      right_pad = x[seq(n, n - nr + 1, by = -1)]
      if (n_right > n) right_pad = c(right_pad, rep(x[1], n_right - n))
    }
  }
  c(left_pad, x, right_pad)
}

#' Applies lifting filter with boundary support (C++ Wrapper)
#'
#' This function bridges the R logic to the compiled C++ core for speed.
#'
#' @param x Input signal.
#' @param coeffs Coefficients of filter P or U.
#' @param start_idx Relative offset.
#' @param extension Extension mode.
#'
#' @keywords internal
.apply_filter_lifting = function(
    x,
    coeffs,
    start_idx = 0,
    extension = "symmetric"
) {

  # Map string to integer for C++ switch-case
  mode_int = switch(extension,
                    "symmetric" = 1L,
                    "periodic"  = 2L,
                    "zero"      = 3L,
                    1L) # Default symmetric

  # Call C++
  # Note: R indexing is 1-based, C++ is 0-based.
  # The old R logic used physical padding which shifted indices.
  # The C++ logic uses virtual padding.
  # Tests indicate start_idx might need +1 or -1 adjustment depending on
  # exact filter alignment. Starting with direct mapping.

  return(apply_filter_cpp(x, coeffs, start_idx, mode_int))
}
