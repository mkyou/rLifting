test_that("new_wavelet_stream warns on update_freq = 0", {
  scheme <- lifting_scheme("haar")
  expect_warning(
    new_wavelet_stream(scheme, window_size = 33, update_freq = 0),
    "update_freq = 0"
  )
})

test_that("denoise_signal_causal warns on update_freq = 0", {
  scheme <- lifting_scheme("haar")
  signal <- rnorm(64)
  expect_warning(
    denoise_signal_causal(signal, scheme, window_size = 33, update_freq = 0),
    "update_freq = 0"
  )
})

test_that("update_freq < 0 raises an error", {
  scheme <- lifting_scheme("haar")
  expect_error(
    new_wavelet_stream(scheme, window_size = 33, update_freq = -1),
    "update_freq"
  )
  expect_error(
    denoise_signal_causal(
      rnorm(64), scheme,
      window_size = 33, update_freq = -1
    ),
    "update_freq"
  )
})

test_that("update_freq = 0 produces finite output without crash (causal batch)", {
  scheme <- lifting_scheme("haar")
  set.seed(3)
  signal <- rnorm(128, sd = 0.5) +
    sin(seq(0, 4 * pi, length.out = 128))
  out <- suppressWarnings(
    denoise_signal_causal(
      signal, scheme,
      window_size = 33, update_freq = 0
    )
  )
  expect_length(out, 128)
  expect_true(all(is.finite(out)))
})

test_that("update_freq = 0 freezes thresholds after warm-up", {
  # With update_freq = 0, thresholds are computed once at warm-up and reused
  # for all subsequent samples; with update_freq = 1 they are recomputed at
  # every step. On a signal whose variance changes over time the two outputs
  # must differ.
  scheme <- lifting_scheme("haar")
  set.seed(4)
  signal <- c(rnorm(64, sd = 0.2), rnorm(64, sd = 2.0))
  out_frozen <- suppressWarnings(
    denoise_signal_causal(
      signal, scheme,
      window_size = 33, update_freq = 0
    )
  )
  out_adaptive <- denoise_signal_causal(
    signal, scheme,
    window_size = 33, update_freq = 1
  )
  expect_false(isTRUE(all.equal(out_frozen, out_adaptive)))
})
