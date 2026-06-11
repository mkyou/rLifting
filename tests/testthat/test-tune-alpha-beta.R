test_that("tune_alpha_beta retorna lista com alpha, beta, sure, converged", {
  set.seed(201)
  x = rnorm(512)
  sch = lifting_scheme("haar")
  res = tune_alpha_beta(x, sch, levels = 3)

  expect_type(res, "list")
  expect_named(res, c("alpha", "beta", "sure", "converged"),
               ignore.order = TRUE)
  expect_true(is.numeric(res$alpha) && length(res$alpha) == 1L)
  expect_true(is.numeric(res$beta) && length(res$beta) == 1L)
  expect_true(is.finite(res$alpha))
  expect_true(is.finite(res$beta))
})

test_that("tune_alpha_beta respeita bounds informados", {
  set.seed(202)
  x = rnorm(512)
  sch = lifting_scheme("haar")
  res = tune_alpha_beta(x, sch, levels = 3,
                        alpha_range = c(0, 5),
                        beta_range = c(0.8, 2.0))

  expect_gte(res$alpha, 0)
  expect_lte(res$alpha, 5)
  expect_gte(res$beta, 0.8)
  expect_lte(res$beta, 2.0)
})

test_that("tune_alpha_beta reduz SURE vs defaults em sinal com estrutura", {
  set.seed(203)
  n = 1024
  pure = c(rep(0, n/2), rep(3, n/2))  # piecewise constant
  noisy = pure + rnorm(n, sd = 0.3)
  sch = lifting_scheme("haar")

  # SURE com defaults
  sure_default = rLifting:::.sure_alpha_beta(noisy, sch, levels = 3,
                                             alpha = 0.3, beta = 1.2)
  # SURE tunada
  tuned = tune_alpha_beta(noisy, sch, levels = 3)

  expect_lt(tuned$sure, sure_default)
})

test_that("tune_alpha_beta produz MSE menor ou comparável ao default em sinal real", {
  set.seed(204)
  n = 1024
  pure = sin(seq(0, 6 * pi, length.out = n))
  noisy = pure + rnorm(n, sd = 0.3)
  sch = lifting_scheme("cdf53")

  tuned = tune_alpha_beta(noisy, sch, levels = 4)
  out_def = denoise_signal_offline(noisy, sch, levels = 4,
                                   alpha = 0.3, beta = 1.2,
                                   shrinkage = "soft")
  out_tun = denoise_signal_offline(noisy, sch, levels = 4,
                                   alpha = tuned$alpha, beta = tuned$beta,
                                   shrinkage = "soft")
  mse_def = mean((pure - out_def)^2)
  mse_tun = mean((pure - out_tun)^2)
  # Tuned não deve piorar MSE em mais de 5% (margem para variância amostral)
  expect_lt(mse_tun, mse_def * 1.05)
})

test_that("tune_alpha_beta emite warning quando optim não converge", {
  local_mocked_bindings(
    optim = function(...) {
      list(par = c(0.3, 1.2), value = 1.0, convergence = 1L)
    },
    .package = "stats"
  )
  set.seed(206)
  x = rnorm(256)
  sch = lifting_scheme("haar")
  expect_warning(tune_alpha_beta(x, sch, levels = 2), "did not converge")
})

test_that("tune_alpha_beta não emite warning quando converge", {
  set.seed(207)
  x = rnorm(512)
  sch = lifting_scheme("haar")
  expect_no_warning(tune_alpha_beta(x, sch, levels = 3))
})

test_that(".sure_alpha_beta é determinístico", {
  set.seed(205)
  x = rnorm(256)
  sch = lifting_scheme("haar")
  s1 = rLifting:::.sure_alpha_beta(x, sch, levels = 3, alpha = 0.5, beta = 1.0)
  s2 = rLifting:::.sure_alpha_beta(x, sch, levels = 3, alpha = 0.5, beta = 1.0)
  expect_equal(s1, s2)
})
