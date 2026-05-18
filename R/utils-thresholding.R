.VALID_THRESHOLD_METHODS = c("universal", "sure")

.sure_optimal_lambda_level = function(d, sigma) {
  n = length(d)
  if (n == 0L || sigma < 1e-15) return(0)
  abs_d = sort(abs(d))
  sigma_sq = sigma * sigma
  universal_cap = sigma * sqrt(2 * log(n))
  cumsum_sq = c(0, cumsum(abs_d * abs_d))
  n_sigma_sq = n * sigma_sq

  best_sure = n_sigma_sq
  best_lambda = 0

  for (k in seq_len(n)) {
    lam = abs_d[k]
    sure = n_sigma_sq + cumsum_sq[k + 1L] +
      lam * lam * (n - k) - 2 * sigma_sq * k
    if (sure < best_sure) {
      best_sure = sure
      best_lambda = lam
    }
  }
  min(best_lambda, universal_cap)
}
.VALID_SHRINKAGE = c("hard", "soft", "semisoft", "scad")

.resolve_shrinkage_args = function(method, shrinkage, threshold_method) {
  if (!is.null(method) && !is.null(shrinkage)) {
    stop("Pass either `method` (deprecated) or `shrinkage`, not both.",
         call. = FALSE)
  }
  if (!is.null(method)) {
    warning(
      "Argument `method` is deprecated; use `shrinkage` instead.",
      call. = FALSE
    )
    shrinkage = method
  }
  if (is.null(shrinkage)) shrinkage = "semisoft"

  if (!is.character(threshold_method) || length(threshold_method) != 1L ||
        !(threshold_method %in% .VALID_THRESHOLD_METHODS)) {
    stop(sprintf(
      "Invalid `threshold_method`: must be one of %s.",
      paste(shQuote(.VALID_THRESHOLD_METHODS), collapse = ", ")
    ), call. = FALSE)
  }
  if (!is.character(shrinkage) || length(shrinkage) != 1L ||
        !(shrinkage %in% .VALID_SHRINKAGE)) {
    stop(sprintf(
      "Invalid `shrinkage`: must be one of %s.",
      paste(shQuote(.VALID_SHRINKAGE), collapse = ", ")
    ), call. = FALSE)
  }

  list(threshold_method = threshold_method, shrinkage = shrinkage)
}
