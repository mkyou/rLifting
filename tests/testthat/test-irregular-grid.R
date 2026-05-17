
# ── Helpers ───────────────────────────────────────────────────────────────────

make_irregular = function(n, seed = 42) {
  set.seed(seed)
  cumsum(c(0, abs(rnorm(n - 1, mean = 1, sd = 0.4))))
}

make_regular = function(n) seq(0, n - 1, by = 1)

# ── Degree auto-inference ─────────────────────────────────────────────────────

test_that("lift_step infers degree automatically for interpolating predict", {
  s_haar  = lift_step("predict", coeffs = c(1),         start_idx = 0)
  s_lin   = lift_step("predict", coeffs = c(0.5, 0.5),  start_idx = 0)
  s_cubic = lift_step("predict", coeffs = c(-1/16, 9/16, 9/16, -1/16), start_idx = -1)
  s_upd   = lift_step("update",  coeffs = c(0.5),       start_idx = 0)
  s_db2   = lift_step("predict", coeffs = c(sqrt(3)),   start_idx = 0)

  expect_equal(s_haar$degree,  0L)
  expect_equal(s_lin$degree,   1L)
  expect_equal(s_cubic$degree, 3L)
  expect_equal(s_upd$degree,  -1L)  # update steps: never inferred
  expect_equal(s_db2$degree,  -1L)  # sum != 1
})

test_that("lift_step respects explicit degree override", {
  s = lift_step("predict", coeffs = c(0.5, 0.5), start_idx = 0, degree = 0L)
  expect_equal(s$degree, 0L)
})

# ── Regular grid equivalence ──────────────────────────────────────────────────

test_that("lwt with equispaced t gives identical result to lwt without t", {
  for (wname in c("haar", "cdf53", "dd4")) {
    sch = lifting_scheme(wname)
    set.seed(1)
    x = rnorm(64)
    t = make_regular(64)

    r_notl = lwt(x, sch, levels = 2)
    r_witht = lwt(x, sch, levels = 2, t = t)

    expect_equal(r_notl$coeffs, r_witht$coeffs,
      tolerance = 1e-12,
      label = paste("regular equivalence failed:", wname))
  }
})

# ── Perfect reconstruction with irregular t ───────────────────────────────────

test_that("ilwt(lwt(x, t=irr)) == x for interpolating wavelets", {
  for (wname in c("haar", "cdf53", "dd4")) {
    sch = lifting_scheme(wname)
    for (n in c(32, 64)) {
      set.seed(7)
      x = rnorm(n)
      t = make_irregular(n)

      res = lwt(x, sch, levels = 2, t = t)
      rec = ilwt(res)
      err = max(abs(x - rec))
      expect_lt(err, 1e-9,
        label = sprintf("PR failed: %s n=%d err=%.2e", wname, n, err))
    }
  }
})

# ── Warning for non-interpolating wavelets ────────────────────────────────────

test_that("lwt warns when t provided for non-interpolating wavelets", {
  for (wname in c("db2", "cdf97")) {
    sch = lifting_scheme(wname)
    x   = rnorm(64)
    t   = make_irregular(64)
    expect_warning(lwt(x, sch, levels = 2, t = t),
      regexp = "irregular",
      label = paste("expected warning for", wname))
  }
})

test_that("lwt with t on non-interpolating wavelet still gives perfect reconstruction", {
  for (wname in c("db2", "cdf97")) {
    sch = lifting_scheme(wname)
    set.seed(3); x = rnorm(64)
    t = make_irregular(64)
    suppressWarnings({
      res = lwt(x, sch, levels = 2, t = t)
      rec = ilwt(res)
    })
    err = max(abs(x - rec))
    expect_lt(err, 1e-9,
      label = sprintf("PR failed (fallback): %s err=%.2e", wname, err))
  }
})

# ── Irregular grid improves denoising vs ignoring positions ──────────────────

test_that("linear signal: position-aware CDF 5/3 reduces spurious detail energy", {
  # CDF 5/3 has 1 vanishing moment: annihilates linear polynomials.
  # With correct positions, interior detail coefficients should be ~0.
  # With fixed (0.5,0.5) on irregular grid, interior coefficients are NOT ~0.
  # Boundary positions may still have small non-zero values due to boundary
  # handling — we test the ratio, not an absolute threshold.
  n   = 256
  t   = make_irregular(n, seed = 7)
  t_n = (t - min(t)) / (max(t) - min(t))
  x   = 3 * t_n + 1   # purely linear

  sch = lifting_scheme("cdf53")

  res_correct = lwt(x, sch, levels = 2, t = t)
  res_ignore  = lwt(x, sch, levels = 2)

  # Interior coefficients (trim 2 boundary positions each side)
  d1_correct = res_correct$coeffs$d1
  d1_ignore  = res_ignore$coeffs$d1
  k = 2L
  d1_correct_int = d1_correct[(k+1):(length(d1_correct)-k)]
  d1_ignore_int  = d1_ignore[(k+1):(length(d1_ignore)-k)]

  energy_correct = sum(d1_correct_int^2)
  energy_ignore  = sum(d1_ignore_int^2)

  # Position-aware → interior energy near-machine-precision
  expect_lt(energy_correct, 1e-20)
  # Ignoring positions → significant spurious energy
  expect_gt(energy_ignore,  1e-4)
})

test_that("denoising with correct t outperforms ignoring positions (irregular grid)", {
  set.seed(42)
  n   = 256
  t   = make_irregular(n, seed = 42)
  t_n = (t - min(t)) / (max(t) - min(t))

  # Piecewise linear signal — clear advantage for position-aware CDF 5/3
  pure  = ifelse(t_n < 0.5, 4 * t_n, 4 - 4 * t_n) + 0.5 * t_n
  noisy = pure + rnorm(n, sd = 0.3)

  sch = lifting_scheme("cdf53")
  lvl = 4L

  r_correct = denoise_signal_offline(noisy, sch, levels = lvl, t = t)
  r_ignore  = denoise_signal_offline(noisy, sch, levels = lvl)

  mse_correct = mean((pure - r_correct)^2)
  mse_ignore  = mean((pure - r_ignore)^2)

  expect_lt(mse_correct, mse_ignore,
    label = sprintf("MSE correct=%.5f should be < ignore=%.5f",
                    mse_correct, mse_ignore))
})

# ── Custom wavelet with explicit degree ───────────────────────────────────────

test_that("custom wavelet with explicit degree works on irregular grid", {
  p = lift_step("predict", coeffs = c(0.5, 0.5), start_idx = 0, degree = 1L)
  u = lift_step("update",  coeffs = c(0.25, 0.25), start_idx = -1)
  w = custom_wavelet("CustomLinear", list(p, u), c(sqrt(2), 1/sqrt(2)))

  set.seed(5); x = rnorm(64)
  t = make_irregular(64)

  res = lwt(x, w, levels = 2, t = t)
  rec = ilwt(res)
  expect_lt(max(abs(x - rec)), 1e-9)
})

# ── API: t parameter accepted by all public functions ─────────────────────────

test_that("t parameter accepted by denoise_signal_offline", {
  sch = lifting_scheme("cdf53")
  x   = rnorm(128)
  t   = make_irregular(128)
  res = denoise_signal_offline(x, sch, levels = 4, t = t)
  expect_equal(length(res), 128)
  expect_false(any(is.na(res)))
})

# ── Causal mode with irregular t ──────────────────────────────────────────────

test_that("denoise_signal_causal accepts t for irregular grid", {
  sch = lifting_scheme("cdf53")
  set.seed(9); x = rnorm(128)
  t = make_irregular(128)
  res = denoise_signal_causal(x, sch, levels = 2, window_size = 32, t = t)
  expect_equal(length(res), 128)
  expect_false(any(is.na(res)))
})

test_that("new_wavelet_stream closure accepts t_val per sample (irregular mode)", {
  sch  = lifting_scheme("cdf53")
  proc = new_wavelet_stream(sch, window_size = 32, levels = 2, irregular = TRUE)
  set.seed(11); x = rnorm(64)
  t = make_irregular(64)
  res = numeric(64)
  for (i in seq_along(x)) res[i] = proc(x[i], t[i])
  expect_equal(length(res), 64)
  expect_false(any(is.na(res)))
})

test_that("causal with equispaced t gives identical result to causal without t", {
  for (wname in c("haar", "cdf53")) {
    sch = lifting_scheme(wname)
    set.seed(2); x = rnorm(128)
    t   = make_regular(128)

    # Use odd window_size so both paths use the same size (t auto-adjusts even→odd)
    r_notl  = denoise_signal_causal(x, sch, levels = 2, window_size = 33)
    r_witht = denoise_signal_causal(x, sch, levels = 2, window_size = 33, t = t)

    expect_equal(r_notl, r_witht, tolerance = 1e-12,
      label = paste("causal regular equivalence:", wname))
  }
})

test_that("causal position-aware CDF5/3 nearly perfectly reconstructs noiseless linear signal", {
  # With correct positions, detail coefficients of a linear signal are ~0, so
  # the adaptive threshold lambda is ~0 and nothing is shrunk → near-perfect PR.
  # With fixed coefficients, details carry the position-induced residual (~0.003),
  # which the adaptive threshold then zeroes out → signal distortion.
  n   = 256
  t   = make_irregular(n, seed = 42)
  t_n = (t - min(t)) / (max(t) - min(t))
  x   = 3 * t_n + 1   # purely linear, NO noise

  sch = lifting_scheme("cdf53")
  warmup = 63L
  r_correct = denoise_signal_causal(x, sch, levels = 1, window_size = 64, t = t)
  r_ignore  = denoise_signal_causal(x, sch, levels = 1, window_size = 64)

  # Exclude warmup (raw pass-through identical for both)
  x_post      = x[-(1:warmup)]
  mse_correct = mean((x_post - r_correct[-(1:warmup)])^2)
  mse_ignore  = mean((x_post - r_ignore[-(1:warmup)])^2)

  expect_lt(mse_correct, mse_ignore,
    label = sprintf("causal noiseless MSE: correct=%.2e vs ignore=%.2e",
                    mse_correct, mse_ignore))
  expect_lt(mse_correct, 1e-6,
    label = sprintf("causal near-PR: mse_correct=%.2e", mse_correct))
})
