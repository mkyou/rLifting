test_that("Haar native verification", {
  sch = lifting_scheme("haar")

  pr = validate_perfect_reconstruction(sch)
  expect_true(pr$passed, label = pr$msg)

  vm0 = validate_vanishing_moments(sch, degree = 0)
  expect_true(vm0$passed, label = vm0$msg)

  ortho = validate_orthogonality(sch, expected = TRUE)
  expect_true(ortho$passed, label = ortho$msg)
})

test_that("CDF 5/3 (LeGall) verification", {
  sch = lifting_scheme("cdf53")

  expect_true(validate_perfect_reconstruction(sch)$passed)

  expect_true(validate_vanishing_moments(sch, degree = 0)$passed)
  expect_true(validate_vanishing_moments(sch, degree = 1)$passed)

  expect_false(validate_vanishing_moments(sch, degree = 2)$passed)

  expect_true(validate_orthogonality(sch, expected = FALSE)$passed)
})

test_that("CDF 9/7 (Cohen-Daubechies-Feauveau) verification", {
  sch = lifting_scheme("cdf97")

  expect_true(validate_perfect_reconstruction(sch)$passed)

  for (d in 0:3) {
    res = validate_vanishing_moments(sch, degree = d)
    expect_true(res$passed, label = paste("Falha no grau", d, "-", res$msg))
  }
})

test_that("DD4 (Interpolating Cubic) verification", {
  sch = lifting_scheme("dd4")
  expect_true(validate_perfect_reconstruction(sch)$passed)

  expect_true(validate_vanishing_moments(sch, degree = 3)$passed)
})

test_that("Lazy Wavelet verification", {
  sch = lifting_scheme("lazy")
  expect_true(validate_perfect_reconstruction(sch)$passed)
  expect_false(validate_vanishing_moments(sch, degree = 0)$passed)
})
