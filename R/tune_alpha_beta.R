#' Tune Adaptive Threshold Parameters via SURE
#'
#' Selects the recursive-threshold parameters \code{alpha} and \code{beta}
#' that minimise Stein's Unbiased Risk Estimate (SURE) for soft thresholding
#' applied to the wavelet detail coefficients of the supplied signal.
#'
#' SURE is computed under the soft-threshold estimator assumption (Donoho &
#' Johnstone, 1995). The chosen parameters are typically usable with
#' \code{shrinkage = "soft"}, \code{"semisoft"}, \code{"hard"}, or
#' \code{"scad"}, since all four share the same threshold location.
#'
#' @param signal Numeric vector.
#' @param scheme A \code{lifting_scheme} object.
#' @param levels Decomposition depth.
#' @param extension Boundary mode (passed to \code{lwt}).
#' @param ll_k Local-linear neighborhood size (passed to \code{lwt}).
#' @param alpha_range Bounds for \code{alpha}; default \code{c(0, 10)}.
#' @param beta_range Bounds for \code{beta}; default \code{c(0.5, 3.0)}.
#'
#' @return A list with components \code{alpha}, \code{beta}, \code{sure}
#'   (the minimised SURE value), and \code{converged} (logical).
#' @export
tune_alpha_beta = function(signal, scheme, levels = 3,
                           extension = "symmetric", ll_k = 4L,
                           alpha_range = c(0, 10),
                           beta_range = c(0.5, 3.0)) {

  lw = lwt(signal, scheme, levels = levels,
           extension = extension, ll_k = ll_k)
  details = lw$coeffs[seq_len(levels)]

  d1 = details[[1L]]
  sigma_hat = stats::mad(d1, constant = 1.4826)
  n_finest = length(d1)

  objective = function(par) {
    .sure_from_details(details, sigma_hat, n_finest,
                       alpha = par[1L], beta = par[2L])
  }

  # Phase 1: coarse grid (handles non-smooth SURE objective).
  alpha_grid = seq(alpha_range[1L], alpha_range[2L], length.out = 21L)
  beta_grid = seq(beta_range[1L], beta_range[2L], length.out = 21L)
  best_val = Inf
  best_par = c(alpha_grid[1L], beta_grid[1L])
  for (a in alpha_grid) {
    for (b in beta_grid) {
      v = objective(c(a, b))
      if (v < best_val) {
        best_val = v
        best_par = c(a, b)
      }
    }
  }

  # Phase 2: Nelder-Mead refinement from grid minimum.
  opt = tryCatch(
    stats::optim(
      par = best_par,
      fn = objective,
      method = "Nelder-Mead",
      control = list(reltol = 1e-8)
    ),
    error = function(e) list(par = best_par, value = best_val, convergence = 99L)
  )

  # Clip refined value to bounds (Nelder-Mead is unconstrained).
  alpha_opt = max(alpha_range[1L], min(alpha_range[2L], opt$par[1L]))
  beta_opt  = max(beta_range[1L],  min(beta_range[2L],  opt$par[2L]))
  # Use refined value only if it's still inside bounds and improves on grid.
  refined_val = if (alpha_opt == opt$par[1L] && beta_opt == opt$par[2L]) {
    opt$value
  } else {
    objective(c(alpha_opt, beta_opt))
  }
  if (refined_val > best_val) {
    alpha_opt = best_par[1L]
    beta_opt = best_par[2L]
    refined_val = best_val
  }

  list(
    alpha = unname(alpha_opt),
    beta = unname(beta_opt),
    sure = refined_val,
    converged = opt$convergence == 0L
  )
}

.lambda_recursion = function(sigma_hat, n_finest, levels, alpha, beta) {
  lam = numeric(levels)
  lam[1L] = beta * sigma_hat * sqrt(2 * log(n_finest))
  if (levels >= 2L) {
    for (k in 2:levels) {
      lam[k] = lam[k - 1L] * (k - 1L) / (k + alpha - 1L)
    }
  }
  lam
}

.sure_from_details = function(details, sigma_hat, n_finest,
                              alpha, beta) {
  levels = length(details)
  lam = .lambda_recursion(sigma_hat, n_finest, levels, alpha, beta)

  total = 0
  for (j in seq_len(levels)) {
    d = details[[j]]
    n_j = length(d)
    # Per-level sigma_j for risk evaluation (the pipeline still uses
    # sigma_hat from d1 to set lam, but SURE measures the actual risk).
    sigma_j = stats::mad(d, constant = 1.4826)
    if (sigma_j < 1e-15) sigma_j = sigma_hat
    sigma_j_sq = sigma_j * sigma_j

    d_sq = d * d
    below = abs(d) <= lam[j]
    total = total +
      n_j * sigma_j_sq +
      sum(pmin(d_sq, lam[j] * lam[j])) -
      2 * sigma_j_sq * sum(below)
  }
  total
}

.sure_alpha_beta = function(signal, scheme, levels = 3,
                            extension = "symmetric", ll_k = 4L,
                            alpha = 0.3, beta = 1.2) {
  lw = lwt(signal, scheme, levels = levels,
           extension = extension, ll_k = ll_k)
  details = lw$coeffs[seq_len(levels)]
  sigma_hat = stats::mad(details[[1L]], constant = 1.4826)
  n_finest = length(details[[1L]])
  .sure_from_details(details, sigma_hat, n_finest, alpha, beta)
}
