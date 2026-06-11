# SCAD shrinkage depends on `a` only when detail coefficients fall in the
# transitional region (2*lambda, a*lambda]. The signals below are sized so
# that this region is non-empty for both shrinkage methods compared.

test_that("denoise_signal_causal honours SCAD `a` argument", {
  scheme <- lifting_scheme("haar")
  set.seed(1)
  signal <- 3 * sin(seq(0, 8 * pi, length.out = 256)) +
    rnorm(256, sd = 0.5)
  out_narrow <- denoise_signal_causal(
    signal, scheme,
    shrinkage = "scad", a = 2.1, beta = 0.4,
    window_size = 33, update_freq = 1
  )
  out_wide <- denoise_signal_causal(
    signal, scheme,
    shrinkage = "scad", a = 50, beta = 0.4,
    window_size = 33, update_freq = 1
  )
  expect_false(isTRUE(all.equal(out_narrow, out_wide)))
})

test_that("new_wavelet_stream honours SCAD `a` argument", {
  scheme <- lifting_scheme("haar")
  set.seed(1)
  signal <- 3 * sin(seq(0, 8 * pi, length.out = 256)) +
    rnorm(256, sd = 0.5)
  stream_narrow <- new_wavelet_stream(
    scheme,
    shrinkage = "scad", a = 2.1, beta = 0.4,
    window_size = 33, update_freq = 1
  )
  stream_wide <- new_wavelet_stream(
    scheme,
    shrinkage = "scad", a = 50, beta = 0.4,
    window_size = 33, update_freq = 1
  )
  out_narrow <- vapply(signal, stream_narrow, numeric(1))
  out_wide <- vapply(signal, stream_wide, numeric(1))
  expect_false(isTRUE(all.equal(out_narrow, out_wide)))
})

test_that("denoise_signal_offline honours SCAD `a` argument", {
  scheme <- lifting_scheme("haar")
  set.seed(1)
  signal <- 3 * sin(seq(0, 8 * pi, length.out = 256)) +
    rnorm(256, sd = 0.5)
  out_narrow <- denoise_signal_offline(
    signal, scheme,
    shrinkage = "scad", a = 2.1, beta = 0.4
  )
  out_wide <- denoise_signal_offline(
    signal, scheme,
    shrinkage = "scad", a = 50, beta = 0.4
  )
  expect_false(isTRUE(all.equal(out_narrow, out_wide)))
})

test_that("default SCAD `a` preserves prior behaviour", {
  scheme <- lifting_scheme("haar")
  set.seed(2)
  signal <- rnorm(64)
  out_implicit <- denoise_signal_causal(
    signal, scheme,
    shrinkage = "scad", window_size = 33, update_freq = 1
  )
  out_explicit <- denoise_signal_causal(
    signal, scheme,
    shrinkage = "scad", a = 3.7, window_size = 33, update_freq = 1
  )
  expect_equal(out_implicit, out_explicit)
})
