test_that("Haar native verification", {
  sch = lifting_scheme("haar")

  # 1. Reconstrucao Perfeita
  pr = validate_perfect_reconstruction(sch)
  expect_true(pr$passed, label = pr$msg)

  # 2. Momentos Nulos (Haar tem 1 momento nulo: Grau 0)
  vm0 = validate_vanishing_moments(sch, degree = 0)
  expect_true(vm0$passed, label = vm0$msg)

  # 3. Ortogonalidade
  ortho = validate_orthogonality(sch, expected = TRUE)
  expect_true(ortho$passed, label = ortho$msg)
})

test_that("CDF 5/3 (LeGall) verification", {
  sch = lifting_scheme("cdf53")

  # Reconstrucao Perfeita
  expect_true(validate_perfect_reconstruction(sch)$passed)

  # Momentos Nulos (Deve passar no Grau 0 e 1)
  expect_true(validate_vanishing_moments(sch, degree = 0)$passed)
  expect_true(validate_vanishing_moments(sch, degree = 1)$passed)

  # Deve FALHAR no Grau 2 (Spline Linear nao faz curva)
  expect_false(validate_vanishing_moments(sch, degree = 2)$passed)

  # Nao deve ser ortogonal (Biorrtogonal)
  expect_true(validate_orthogonality(sch, expected = FALSE)$passed)
})

test_that("CDF 9/7 (Cohen-Daubechies-Feauveau) verification", {
  sch = lifting_scheme("cdf97")

  # Stress test matematico
  expect_true(validate_perfect_reconstruction(sch)$passed)

  # Deve ter 4 momentos nulos (0 a 3)
  for(d in 0:3) {
    res = validate_vanishing_moments(sch, degree = d)
    expect_true(res$passed, label = paste("Falha no grau", d, "-", res$msg))
  }
})

test_that("DD4 (Interpolating Cubic) verification", {
  sch = lifting_scheme("dd4")
  expect_true(validate_perfect_reconstruction(sch)$passed)

  # Interpoladora cubica tem 4 momentos nulos
  expect_true(validate_vanishing_moments(sch, degree = 3)$passed)
})

test_that("Lazy Wavelet verification", {
  sch = lifting_scheme("lazy")
  expect_true(validate_perfect_reconstruction(sch)$passed)
  # Lazy nao deve ter momentos nulos
  expect_false(validate_vanishing_moments(sch, degree = 0)$passed)
})
