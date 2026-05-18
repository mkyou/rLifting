test_that("threshold_scad zera coeficientes |d| <= lambda", {
  x = c(0, 0.5, -0.5, 1.0, -1.0)
  expect_equal(threshold_scad(x, lambda = 1.0), rep(0, length(x)))
})

test_that("threshold_scad: regiao soft em lambda < |d| <= 2*lambda", {
  # lambda=1, a=3.7. d=1.5 -> sign(d)*(|d|-lambda) = 0.5
  expect_equal(threshold_scad(1.5, lambda = 1.0), 0.5)
  expect_equal(threshold_scad(-1.5, lambda = 1.0), -0.5)
  # boundary 2*lambda: continuidade com proxima regiao
  expect_equal(threshold_scad(2.0, lambda = 1.0), 1.0)
})

test_that("threshold_scad: regiao SCAD em 2*lambda < |d| <= a*lambda", {
  # a=3.7 (default). lambda=1, d=2.5
  # ((a-1)*d - sign(d)*a*lambda) / (a-2)
  expected = ((3.7 - 1) * 2.5 - 3.7) / (3.7 - 2)
  expect_equal(threshold_scad(2.5, lambda = 1.0), expected)
  expect_equal(threshold_scad(-2.5, lambda = 1.0), -expected)
})

test_that("threshold_scad: regiao identidade em |d| > a*lambda", {
  # a=3.7, lambda=1, d=5 -> 5 (sem viés)
  expect_equal(threshold_scad(5.0, lambda = 1.0), 5.0)
  expect_equal(threshold_scad(-5.0, lambda = 1.0), -5.0)
  # boundary a*lambda: continuidade
  expect_equal(threshold_scad(3.7, lambda = 1.0), 3.7)
})

test_that("threshold_scad: parametro a customizavel", {
  # a=4, lambda=1, d=2.5 -> ((4-1)*2.5 - 4)/(4-2) = (7.5-4)/2 = 1.75
  expect_equal(threshold_scad(2.5, lambda = 1.0, a = 4.0), 1.75)
})

test_that("threshold_scad: a <= 2 da erro", {
  expect_error(threshold_scad(1.0, lambda = 1.0, a = 2.0), "a")
  expect_error(threshold_scad(1.0, lambda = 1.0, a = 1.5), "a")
})

test_that("threshold dispatcher aceita method='scad'", {
  x = c(0, 0.5, 1.5, 2.5, 5.0)
  out = threshold(x, lambda = 1.0, method = "scad")
  expected = c(0, 0, 0.5, ((3.7-1)*2.5 - 3.7)/(3.7-2), 5.0)
  expect_equal(out, expected)
})

test_that("denoise_signal_offline aceita shrinkage='scad'", {
  set.seed(101)
  x = rnorm(128)
  sch = lifting_scheme("haar")
  expect_no_error(
    denoise_signal_offline(x, sch, levels = 2, shrinkage = "scad")
  )
})

test_that("denoise_signal_causal aceita shrinkage='scad'", {
  set.seed(102)
  x = rnorm(256)
  sch = lifting_scheme("haar")
  expect_no_error(
    denoise_signal_causal(x, sch, levels = 2, window_size = 64,
                          shrinkage = "scad")
  )
})

test_that("new_wavelet_stream aceita shrinkage='scad'", {
  sch = lifting_scheme("haar")
  expect_no_error(
    new_wavelet_stream(sch, window_size = 64, levels = 2,
                       shrinkage = "scad")
  )
})

test_that("SCAD reduz viés vs soft em coeficientes grandes (offline)", {
  set.seed(103)
  # Sinal com coeficientes grandes e ruído pequeno
  n = 512
  pure = sin(seq(0, 4 * pi, length.out = n)) * 3
  noisy = pure + rnorm(n, sd = 0.2)
  sch = lifting_scheme("haar")

  out_soft = denoise_signal_offline(noisy, sch, levels = 3, shrinkage = "soft")
  out_scad = denoise_signal_offline(noisy, sch, levels = 3, shrinkage = "scad")

  # SCAD não encolhe coeficientes grandes -> MSE menor ou igual em sinal com sinal forte
  mse_soft = mean((pure - out_soft)^2)
  mse_scad = mean((pure - out_scad)^2)
  expect_lt(mse_scad, mse_soft * 1.1)  # margem
})
