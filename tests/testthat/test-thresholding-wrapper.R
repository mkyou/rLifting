
test_that("General threshold wrapper dispatches correctly", {
  x = c(-2, -0.5, 0, 0.5, 2)
  lambda = 1

  # Soft
  res_soft = threshold(x, lambda, method = "soft")
  exp_soft = threshold_soft(x, lambda)
  expect_equal(res_soft, exp_soft)

  # Hard
  res_hard = threshold(x, lambda, method = "hard")
  exp_hard = threshold_hard(x, lambda)
  expect_equal(res_hard, exp_hard)

  # Semisoft
  res_semisoft = threshold(x, lambda, method = "semisoft")
  exp_semisoft = threshold_semisoft(x, lambda)
  expect_equal(res_semisoft, exp_semisoft)

  # Error handling
  expect_error(threshold(x, lambda, method = "invalid_method"))
})
