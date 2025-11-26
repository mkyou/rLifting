test_that("Offline Denoising returns correct dimensions", {
  x = rnorm(512)
  sch = lifting_scheme("db2")

  # Testa Hard, Soft e Semisoft
  for (met in c("hard", "soft", "semisoft")) {
    res = denoise_signal_offline(x, sch, method = met)
    expect_equal(length(res), length(x))
    expect_false(any(is.na(res)))
  }
})

test_that("Causal Denoising runs correctly via stream processor", {
  x = rnorm(100)
  sch = lifting_scheme("haar")

  # Inicializa processador
  proc = new_wavelet_stream(sch, window_size = 16)

  # Simula stream
  out = numeric(100)
  for(i in 1:100) {
    val = proc(x[i])
    expect_false(is.na(val))
    out[i] = val
  }

  expect_equal(length(out), 100)
})

test_that("Thresholding functions math check", {
  # Soft Thresholding: |3| - 1 = 2. sign(3) = 1. Result 2.
  expect_equal(threshold_soft(3, 1), 2)
  # Soft: |0.5| - 1 < 0. Result 0.
  expect_equal(threshold_soft(0.5, 1), 0)

  # Hard: |3| > 1. Result 3.
  expect_equal(threshold_hard(3, 1), 3)
  # Hard: |0.5| < 1. Result 0.
  expect_equal(threshold_hard(0.5, 1), 0)
})
