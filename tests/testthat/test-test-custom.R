test_that("Custom Wavelet API works", {
  p1 = lift_step("predict", c(1), position = "center")
  u1 = lift_step("update", c(0.5), position = "center")

  custom_haar = custom_wavelet(
    "ManualHaar",
    list(p1, u1),
    c(sqrt(2), 1 / sqrt(2))
  )

  expect_true(validate_perfect_reconstruction(custom_haar)$passed)

  x = c(1, 2, 3, 4, 5, 6, 7, 8)
  native_haar = lifting_scheme("haar")

  res_custom = lwt(x, custom_haar, levels = 1)
  res_native = lwt(x, native_haar, levels = 1)

  expect_equal(res_custom$coeffs$d1, res_native$coeffs$d1)
  expect_equal(res_custom$coeffs$a1, res_native$coeffs$a1)
})

test_that("lift_step handles index logic correctly", {
  s1 = lift_step("predict", c(1, 1, 1), position = "center")
  expect_equal(s1$start_idx, -1)

  s2 = lift_step("predict", c(1, 1), position = "left")
  expect_equal(s2$start_idx, -1)

  s3 = lift_step("predict", c(1), start_idx = 100)
  expect_equal(s3$start_idx, 100)
})
