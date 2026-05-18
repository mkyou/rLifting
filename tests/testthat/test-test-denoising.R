test_that("Offline Denoising returns correct dimensions", {
  x = rnorm(512)
  sch = lifting_scheme("db2")

  for (met in c("hard", "soft", "semisoft")) {
    res = denoise_signal_offline(x, sch, method = met)
    expect_equal(length(res), length(x))
    expect_false(any(is.na(res)))
  }
})

test_that("Causal Denoising runs correctly via stream processor", {
  x = rnorm(100)
  sch = lifting_scheme("haar")

  proc = new_wavelet_stream(sch, window_size = 16)

  out = numeric(100)
  for (i in 1:100) {
    val = proc(x[i])
    expect_false(is.na(val))
    out[i] = val
  }

  expect_equal(length(out), 100)
})

test_that("Thresholding functions math check", {
  expect_equal(threshold_soft(3, 1), 2)
  expect_equal(threshold_soft(0.5, 1), 0)

  expect_equal(threshold_hard(3, 1), 3)
  expect_equal(threshold_hard(0.5, 1), 0)
})
