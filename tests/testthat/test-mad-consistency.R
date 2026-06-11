test_that("compute_mad_cpp returns centered MAD for even n", {
  # median(c(-1,-2,-3,-4)) = -2.5; |deviations| = c(1.5, 0.5, 0.5, 1.5)
  # sorted: 0.5, 0.5, 1.5, 1.5 -> median = (0.5+1.5)/2 = 1.0
  expect_equal(rLifting:::compute_mad_cpp(c(-1, -2, -3, -4)), 1.0)
  expect_equal(rLifting:::compute_mad_cpp(c(1, 2, 3, 4)), 1.0)
})

test_that("compute_mad_cpp returns centered MAD for odd n", {
  # median(c(1,2,3,4,5)) = 3; |deviations| = c(2,1,0,1,2); median = 1
  expect_equal(rLifting:::compute_mad_cpp(c(1, 2, 3, 4, 5)), 1.0)
})

test_that("compute_mad_cpp returns 0 for empty input", {
  expect_equal(rLifting:::compute_mad_cpp(numeric(0)), 0)
})

test_that("compute_mad_cpp matches mad(x, constant=1) exactly", {
  set.seed(42)
  x <- rnorm(99, mean = 5, sd = 2)
  expect_equal(rLifting:::compute_mad_cpp(x), mad(x, constant = 1))
})

test_that("compute_mad_cpp handles non-zero-centered input correctly", {
  # Key case: old implementation would return median(|x|)=12; new returns MAD=1
  x <- c(10, 11, 12, 13, 14)
  expect_equal(rLifting:::compute_mad_cpp(x), 1.0)
  expect_equal(rLifting:::compute_mad_cpp(x), mad(x, constant = 1))
})

test_that("compute_thresholds_cpp uses centered MAD on d1", {
  d1 <- c(-1, -2, 3, 4)
  # median(d1) = 1; |deviations| = c(2,3,2,3); median = 2.5
  expected_sigma <- 2.5 / 0.6745
  expected_lambda1 <- expected_sigma * sqrt(2 * log(length(d1)))
  lambdas <- rLifting:::compute_thresholds_cpp(d1, 1L, 0, 1)
  expect_equal(lambdas[1], expected_lambda1)
})

test_that("denoise_signal_offline sigma matches centered MAD on d1", {
  # Lazy wavelet: d1 == odd subsamples = c(1,2,3,4)
  # centered MAD: median=2.5, |deviations|=c(1.5,0.5,0.5,1.5), MAD=1.0
  scheme <- lifting_scheme("lazy")
  signal <- c(0, 1, 0, 2, 0, 3, 0, 4)
  d1 <- signal[seq(2, length(signal), by = 2)]
  expected_sigma <- mad(d1, constant = 1) / 0.6745
  lambdas_universal <- rLifting:::compute_thresholds_cpp(d1, 1L, 0, 1)
  expect_equal(lambdas_universal[1],
               expected_sigma * sqrt(2 * log(length(d1))))
})
