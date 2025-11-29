#' Generate standard signals for wavelet testing (Donoho-Johnstone Benchmark)
#'
#' Creates synthetic classic signals used for wavelet validation.
#' Formulas are based on Donoho and Johnstone and Liu et al.
#'
#' @param type Signal type: "const", "ramp", "poly2", "poly3", "random",
#'             "impulse", "sine", "doppler", "heavisine", "bumps".
#' @param n Signal length (default 512).
#'
#' @return Numeric vector.
#' @keywords internal
.generate_signal = function(type, n = 512) {
  t = seq(0, 1, length.out = n)

  if (type == "const")   return(rep(100, n))
  if (type == "ramp")    return(seq(1, 100, length.out = n))
  if (type == "poly2")   return(100 * t^2)
  if (type == "poly3")   return(100 * t^3)
  if (type == "random")  return(rnorm(n))
  if (type == "impulse") {
    x = rep(0, n)
    x[floor(n/2)] = 1
    return(x)
  }
  if (type == "sine")    return(sin(4 * pi * t))

  # Advanced Benchmarks
  if (type == "doppler") {
    eps = 0.05
    return(sqrt(t * (1 - t)) * sin((2 * pi * 1.05) / (t + eps)))
  }

  if (type == "heavisine") {
    return(4 * sin(4 * pi * t) - sign(t - 0.3) - sign(0.72 - t))
  }

  if (type == "bumps") {
    pos = c(0.1, 0.13, 0.15, 0.23, 0.25, 0.40, 0.44, 0.65, 0.76, 0.78, 0.81)
    h   = c(4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 5.1, -4.2)
    w   = c(0.005, 0.005, 0.006, 0.01, 0.01, 0.03, 0.01, 0.01, 0.005,
            0.008, 0.005)

    x = numeric(n)
    for (j in seq_along(pos)) {
      x = x + h[j] * (1 + abs((t - pos[j]) / w[j])^4)^(-1)
    }
    return(x)
  }

  stop(sprintf("Unknown signal type '%s'.", type))
}

#' Validate Perfect Reconstruction (Stress Test)
#'
#' Verifies wavelet invertibility against a battery of signals.
#'
#' @param scheme Object of class \code{lifting_scheme}.
#' @param tol Numerical error tolerance (default 1e-9).
#'
#' @return List with global status and maximum error found.
#' @export
validate_perfect_reconstruction = function(scheme, tol = 1e-9) {
  signals = c("random", "ramp", "sine", "doppler", "heavisine", "bumps")
  max_err_global = 0
  failed_signals = c()

  for (sig_type in signals) {
    x = .generate_signal(sig_type, n = 512)

    # Use periodic to focus on math, not boundary
    res = lwt(x, scheme, extension = "periodic")
    rec = ilwt(res)

    err = max(abs(x - rec))
    max_err_global = max(max_err_global, err)

    if (err >= tol) {
      failed_signals = c(failed_signals, sig_type)
    }
  }

  passed = length(failed_signals) == 0
  msg = if (passed) {
    sprintf("Max Error (Global): %.2e", max_err_global)
  } else {
    sprintf("FAILED on: %s. Max Error: %.2e",
            paste(failed_signals, collapse=", "), max_err_global)
  }

  list(
    name = "Perfect Reconstruction (Stress Test)",
    passed = passed,
    metric = max_err_global,
    msg = msg
  )
}

#' Validate Vanishing Moments
#'
#' Verifies if the wavelet cancels polynomials of a specific degree.
#'
#' @param scheme Object of class \code{lifting_scheme}.
#' @param degree Polynomial degree (0=Constant, 1=Ramp, 2=Parabola...).
#' @param tol Residual energy tolerance (default 1e-9).
#'
#' @return List with status and residual energy.
#' @export
validate_vanishing_moments = function(scheme, degree = 0, tol = 1e-9) {
  n = 128
  type_map = c("0"="const", "1"="ramp", "2"="poly2", "3"="poly3")
  sig_type = type_map[as.character(degree)]

  if (is.na(sig_type)) return(
    list(name = "VM", passed = FALSE, msg = "Invalid Degree")
  )

  x = .generate_signal(sig_type, n)

  res = lwt(x, scheme, extension = "periodic")
  d1 = res$coeffs$d1

  cut = floor(n * 0.15)
  d_core = d1[(cut+1):(length(d1)-cut)]

  energy = sum(d_core^2)
  passed = energy < tol

  list(
    name = sprintf("Vanishing Moments (Degree %d)", degree),
    passed = passed,
    metric = energy,
    msg = sprintf("Residual Energy: %.2e (Tol: %.2e)", energy, tol)
  )
}

#' Validate Orthogonality (Energy Conservation)
#'
#' Verifies Parseval's Theorem. Only applicable for orthogonal wavelets.
#'
#' @param scheme Object of class \code{lifting_scheme}.
#' @param expected Boolean. If TRUE, expects orthogonality.
#' @param tol Tolerance (default 1e-9).
#'
#' @return List with status and energy ratio (Out/In).
#' @export
validate_orthogonality = function(scheme, expected = TRUE, tol = 1e-9) {
  x = .generate_signal("random", n = 512)
  energy_in = sum(x^2)

  res = lwt(x, scheme, extension = "periodic")
  energy_out = sum(res$coeffs$a1^2) + sum(res$coeffs$d1^2)

  ratio = energy_out / energy_in
  is_ortho = abs(ratio - 1) < tol

  passed = (expected && is_ortho) || (!expected)
  if (expected && !is_ortho) passed = FALSE

  status_str = if (is_ortho) "Orthogonal" else "Non-Orthogonal"

  list(
    name = "Orthogonality (Energy Conservation)",
    passed = passed,
    metric = ratio,
    msg = sprintf("Ratio: %.6f (%s) [Exp: %s]",
                  ratio, status_str, expected)
  )
}

#' Validate Compact Support (FIR Compliance)
#'
#' Verifies if the impulse response is finite (FIR Filter).
#'
#' @param scheme Object of class \code{lifting_scheme}.
#' @param max_width Maximum expected width (number of taps).
#'
#' @return List with status and number of active taps.
#' @export
validate_compact_support = function(scheme, max_width) {
  n = 256
  coeffs_mock = list(a1 = rep(0, n/2), d1 = rep(0, n/2))
  coeffs_mock$d1[n/4] = 1

  mock_lwt = structure(
    list(
      coeffs = coeffs_mock, scheme = scheme,
      levels = 1, original_len = n, extension = "zero"
    ),
    class = "lwt"
  )

  psi = ilwt(mock_lwt)
  nonzero = sum(abs(psi) > 1e-10)

  passed = nonzero <= (max_width + 2) && nonzero >= 1

  list(
    name = "Compact Support (FIR)",
    passed = passed,
    metric = nonzero,
    msg = sprintf("Active Taps: %d (Max Expected: %d)", nonzero, max_width)
  )
}

#' Validate Shift Sensitivity (Shift Variance)
#'
#' Decimated wavelets are not translation invariant.
#' This test quantifies the variation in detail energy when
#' shifting the input signal by 1 sample.
#'
#' @param scheme Object of class \code{lifting_scheme}.
#'
#' @return List with status and percentage variation.
#' @export
validate_shift_sensitivity = function(scheme) {
  n = 128
  x = rep(0, n)
  x[n/2] = 1

  res1 = lwt(x, scheme, levels = 1, extension = "periodic")
  energy1 = sum(res1$coeffs$d1^2)

  x_shift = c(x[n], x[1:(n-1)])
  res2 = lwt(x_shift, scheme, levels = 1, extension = "periodic")
  energy2 = sum(res2$coeffs$d1^2)

  diff_pct = 0
  if (energy1 > 1e-15) {
    diff_pct = abs(energy1 - energy2) / energy1 * 100
  }

  msg = sprintf("Energy Variation (Shift 1): %.2f%%", diff_pct)

  list(
    name = "Shift Sensitivity",
    passed = TRUE,
    metric = diff_pct,
    msg = msg
  )
}

#' Visualize Basis Functions (Scaling and Wavelet)
#'
#' Plots the waveform by iterating the reconstruction over several levels.
#'
#' @param scheme Object of class \code{lifting_scheme}.
#' @param plot Boolean.
#' @param levels Number of cascade levels.
#'
#' @export
visualize_wavelet_basis = function(scheme, plot = TRUE, levels = 8) {
  n = 2^levels

  # 1. Wavelet Function (Psi)
  coeffs_psi = list()
  for (j in 1:levels) coeffs_psi[[paste0("d", j)]] = numeric(n / 2^j)
  coeffs_psi[[paste0("a", levels)]] = numeric(n / 2^levels)
  coeffs_psi[[paste0("d", levels)]][1] = 1

  lwt_psi = structure(
    list(coeffs = coeffs_psi, scheme = scheme, levels = levels,
         original_len = n, extension = "zero"),
    class = "lwt"
  )
  psi = ilwt(lwt_psi)

  # 2. Scaling Function (Phi)
  coeffs_phi = list()
  for (j in 1:levels) coeffs_phi[[paste0("d", j)]] = numeric(n / 2^j)
  coeffs_phi[[paste0("a", levels)]] = numeric(n / 2^levels)
  coeffs_phi[[paste0("a", levels)]][1] = 1

  lwt_phi = structure(
    list(coeffs = coeffs_phi, scheme = scheme, levels = levels,
         original_len = n, extension = "zero"),
    class = "lwt"
  )
  phi = ilwt(lwt_phi)

  if (plot) {
    par(mfrow = c(1, 2))
    ts.plot(phi, main = paste("Scaling (Phi):", scheme$wavelet),
            ylab = "Amp", col = "blue", lwd = 1.5)
    grid()
    ts.plot(psi, main = paste("Wavelet (Psi):", scheme$wavelet),
            ylab = "Amp", col = "red", lwd = 1.5)
    grid()
    par(mfrow = c(1, 1))
  }
}

#' Complete Wavelet Diagnosis
#'
#' Runs a battery of physical and mathematical tests on a wavelet.
#'
#' @param wavelet_name Name string or a \code{lifting_scheme} object.
#' @param config Configuration list (is_ortho, vm_degrees, max_taps).
#' @param verbose Print results to console?
#'
#' @return (Invisible) A list containing the results of each test
#' (Perfect Reconstruction, Orthogonality, etc).
#' Each item contains: \code{passed} (TRUE/FALSE),
#' \code{metric} (numeric), and \code{msg}.
#' @export
diagnose_wavelet = function(wavelet_name, config, verbose = TRUE) {

  if (is.character(wavelet_name)) {
    sch = tryCatch(lifting_scheme(wavelet_name), error = function(e) NULL)
  } else {
    sch = wavelet_name
    wavelet_name = sch$wavelet
  }

  if (is.null(sch)) stop("Invalid Wavelet.")

  if (verbose) {
    cat(sprintf("\n=== DIAGNOSIS: %s ===\n", toupper(wavelet_name)))
  }

  tests = list()
  tests[[1]] = validate_perfect_reconstruction(sch)
  tests[[2]] = validate_orthogonality(sch, expected = config$is_ortho)

  for (deg in config$vm_degrees) {
    idx = length(tests) + 1
    tests[[idx]] = validate_vanishing_moments(sch, degree = deg)
  }

  tests[[length(tests)+1]] = validate_compact_support(
    sch, max_width = config$max_taps
  )

  tests[[length(tests)+1]] = validate_shift_sensitivity(sch)

  if (verbose) {
    cat("Generating basis visualization (check plot window)...\n")
    try(visualize_wavelet_basis(sch, plot = TRUE))
  }

  if (verbose) {
    for (res in tests) {
      status = if(res$passed) "[PASS]" else "[FAIL]"
      if (res$name == "Shift Sensitivity") status = "[INFO]"

      cat(sprintf("%-6s %-35s | %s\n", status, res$name, res$msg))
    }
  }
  invisible(tests)
}
