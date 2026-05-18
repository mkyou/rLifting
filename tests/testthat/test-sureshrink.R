test_that("denoise_signal_offline aceita threshold_method='sure'", {
  set.seed(301)
  x = rnorm(256)
  sch = lifting_scheme("haar")
  expect_no_error(
    denoise_signal_offline(x, sch, levels = 3,
                           threshold_method = "sure",
                           shrinkage = "soft")
  )
})

test_that("denoise_signal_offline com sure combina com qualquer shrinkage", {
  set.seed(302)
  x = rnorm(256)
  sch = lifting_scheme("haar")
  for (sh in c("hard", "soft", "semisoft", "scad")) {
    expect_no_error(
      denoise_signal_offline(x, sch, levels = 3,
                             threshold_method = "sure",
                             shrinkage = sh),
      message = sh
    )
  }
})

test_that("SureShrink em puro ruído Gaussiano produz reconstrução próxima de zero", {
  set.seed(303)
  x = rnorm(1024, sd = 1)
  sch = lifting_scheme("haar")

  out_sure = denoise_signal_offline(x, sch, levels = 4,
                                    threshold_method = "sure",
                                    shrinkage = "soft")
  # Sob H_0 (sinal zero, só ruído), threshold deveria zerar quase tudo
  expect_lt(mean(out_sure^2), 0.4)  # MSE pequeno vs variância do ruído (1)
})

test_that("SureShrink produz output válido em sinal com estrutura", {
  set.seed(304)
  n = 1024
  pure = sin(seq(0, 6 * pi, length.out = n))
  noisy = pure + rnorm(n, sd = 0.3)
  sch = lifting_scheme("haar")

  out_sure = denoise_signal_offline(noisy, sch, levels = 4,
                                    threshold_method = "sure",
                                    shrinkage = "soft")
  # Output finito e dimensão correta
  expect_length(out_sure, n)
  expect_true(all(is.finite(out_sure)))
  # MSE reduz vs sinal ruidoso (denoising acontece de fato)
  mse_noisy = mean((pure - noisy)^2)
  mse_sure = mean((pure - out_sure)^2)
  expect_lt(mse_sure, mse_noisy)
})

test_that(".sure_optimal_lambda_level retorna lambda finito não-negativo", {
  set.seed(305)
  d = rnorm(128, sd = 0.5)
  lam = .sure_optimal_lambda_level(d, sigma = 0.5)
  expect_true(is.finite(lam))
  expect_gte(lam, 0)
})

test_that(".sure_optimal_lambda_level: lambda <= max(|d|) sempre", {
  set.seed(306)
  d = rnorm(128, sd = 0.5)
  lam = .sure_optimal_lambda_level(d, sigma = 0.5)
  expect_lte(lam, max(abs(d)))
})

test_that(".sure_optimal_lambda_level: minimiza SURE entre candidatos", {
  set.seed(307)
  d = c(rnorm(120, sd = 0.5), 3, -3.2, 4.1, 2.8)  # ruído + sinais grandes
  sigma = 0.5
  lam_star = .sure_optimal_lambda_level(d, sigma = sigma)

  sure_fn = function(lambda) {
    sigma_sq = sigma * sigma
    length(d) * sigma_sq +
      sum(pmin(d * d, lambda * lambda)) -
      2 * sigma_sq * sum(abs(d) <= lambda)
  }
  sure_star = sure_fn(lam_star)
  # Testa em algumas candidates fora do ótimo
  for (alt in c(0, sigma, sigma * sqrt(2 * log(length(d))), max(abs(d)))) {
    expect_gte(sure_fn(alt), sure_star - 1e-9,
               label = sprintf("SURE(%.3f) vs SURE*", alt))
  }
})

test_that("threshold_method='sure' com hard shrinkage em sinal sparse", {
  set.seed(308)
  n = 512
  pure = c(rep(0, n - 5), 5, -5, 4, 3, 2)  # very sparse signal
  noisy = pure + rnorm(n, sd = 0.3)
  sch = lifting_scheme("haar")
  out = denoise_signal_offline(noisy, sch, levels = 3,
                               threshold_method = "sure",
                               shrinkage = "hard")
  expect_true(is.numeric(out))
  expect_length(out, n)
})

test_that("denoise_signal_causal aceita threshold_method='sure'", {
  set.seed(309)
  x = rnorm(512)
  sch = lifting_scheme("haar")
  expect_no_error(
    denoise_signal_causal(x, sch, levels = 2, window_size = 128,
                          threshold_method = "sure",
                          shrinkage = "soft")
  )
})

test_that("new_wavelet_stream aceita threshold_method='sure'", {
  sch = lifting_scheme("haar")
  expect_no_error(
    new_wavelet_stream(sch, window_size = 128, levels = 2,
                       threshold_method = "sure",
                       shrinkage = "soft")
  )
})
