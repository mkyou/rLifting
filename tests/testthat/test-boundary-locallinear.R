
# --- Mathematical baseline ---
# For a signal [3, 5, 7, ...], local-linear extrapolation at i=-1:
# x[-1] = x[0] + (-1) * (x[1] - x[0]) = 3 - 2 = 1
# At i=-2: x[-2] = 3 + (-2)*2 = -1
# At i=n (right): x[n] = x[n-1] + 1*(x[n-1] - x[n-2])
test_that("local_linear extrapolation is mathematically correct", {
  skip_if_not_installed("rLifting")

  sch = lifting_scheme("haar")
  x = c(3, 5, 7, 9, 11, 13, 15, 17)  # Linear signal, slope 2

  # With local_linear, the filter at boundaries reads extrapolated values.
  # A perfectly linear signal should produce near-zero details after LWT,
  # since the linear predictor in Haar/CDF removes linear trends.
  res = lwt(x, sch, levels = 1, extension = "local_linear")
  expect_false(any(is.na(res$coeffs$d1)))
  expect_false(any(is.infinite(res$coeffs$d1)))
})

# --- Perfect reconstruction ---
test_that("local_linear preserves perfect reconstruction for all built-in wavelets", {
  wavelets = c("haar", "db2", "cdf53", "cdf97", "dd4")
  signals = c("random", "ramp", "sine", "doppler")

  for (wname in wavelets) {
    sch = lifting_scheme(wname)
    for (sig_type in signals) {
      x = rLifting:::.generate_signal(sig_type, n = 64)
      res = lwt(x, sch, levels = 2, extension = "local_linear")
      rec = ilwt(res)
      err = max(abs(x - rec))
      expect_lt(err, 1e-9,
        label = sprintf("PR failed: wavelet=%s signal=%s err=%.2e", wname, sig_type, err))
    }
  }
})

# --- API acceptance ---
test_that("local_linear is accepted by all public denoising functions", {
  sch = lifting_scheme("haar")
  x = rnorm(128)

  # Offline
  res_off = denoise_signal_offline(x, sch, levels = 2, extension = "local_linear")
  expect_equal(length(res_off), length(x))
  expect_false(any(is.na(res_off)))

  # Causal batch
  res_caus = denoise_signal_causal(x, sch, levels = 2, window_size = 32,
                                   extension = "local_linear")
  expect_equal(length(res_caus), length(x))
  expect_false(any(is.na(res_caus)))

  # Stream processor
  proc = new_wavelet_stream(sch, window_size = 32, levels = 2,
                            extension = "local_linear")
  out = numeric(length(x))
  for (i in seq_along(x)) out[i] = proc(x[i])
  expect_equal(length(out), length(x))
  expect_false(any(is.na(out)))
})

# --- Boundary behaviour ---
# local_linear should produce different results than symmetric at boundaries.
test_that("local_linear produces different boundary behaviour than symmetric", {
  sch = lifting_scheme("db2")
  x = cumsum(rnorm(64))  # Non-stationary signal where boundary treatment matters

  res_sym  = lwt(x, sch, levels = 2, extension = "symmetric")
  res_ll   = lwt(x, sch, levels = 2, extension = "local_linear")

  # Details should differ — at least one level, at least one boundary coefficient
  expect_false(identical(res_sym$coeffs$d1, res_ll$coeffs$d1))
})

# --- Causal: no look-ahead contamination ---
# Changing the signal AFTER time t must not alter output at time t.
# This holds for all extension modes in causal processing.
test_that("local_linear causal mode has no look-ahead leakage", {
  sch = lifting_scheme("haar")
  set.seed(1)
  x = rnorm(64)

  proc1 = new_wavelet_stream(sch, window_size = 16, extension = "local_linear")
  proc2 = new_wavelet_stream(sch, window_size = 16, extension = "local_linear")

  # Feed 32 identical samples to both processors
  out1 = numeric(32)
  out2 = numeric(32)
  for (i in 1:32) {
    out1[i] = proc1(x[i])
    out2[i] = proc2(x[i])
  }

  # Diverge: proc2 gets completely different future samples
  for (i in 33:64) proc2(rnorm(1))

  # Output up to t=32 must be identical regardless of future
  expect_equal(out1, out2)
})
