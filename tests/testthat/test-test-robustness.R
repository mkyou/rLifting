test_that("Stream processor handles NA and Inf gracefully", {
  sch = lifting_scheme("haar")
  proc = new_wavelet_stream(sch, window_size = 16)

  v1 = proc(1)
  expect_equal(v1, 1)

  expect_warning({ val_na = proc(NA) })
  expect_true(is.na(val_na))

  expect_warning({ val_inf = proc(Inf) })
  expect_true(is.infinite(val_inf))
})

test_that("Stream processor rejects bad window sizes", {
  sch = lifting_scheme("haar")
  expect_error(new_wavelet_stream(sch, window_size = 7))
})
