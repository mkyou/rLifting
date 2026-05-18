test_that("denoise_signal_offline aceita shrinkage e threshold_method", {
  set.seed(1)
  x = rnorm(128)
  sch = lifting_scheme("haar")

  expect_no_error(
    denoise_signal_offline(x, sch, levels = 2,
                           threshold_method = "universal",
                           shrinkage = "semisoft")
  )
  expect_no_error(
    denoise_signal_offline(x, sch, levels = 2, shrinkage = "hard")
  )
})

test_that("denoise_signal_offline: shrinkage produz mesmo output que method (deprecated)", {
  set.seed(2)
  x = rnorm(128)
  sch = lifting_scheme("haar")

  for (sh in c("hard", "soft", "semisoft")) {
    new_out = denoise_signal_offline(x, sch, levels = 2, shrinkage = sh)
    old_out = suppressWarnings(
      denoise_signal_offline(x, sch, levels = 2, method = sh)
    )
    expect_equal(new_out, old_out, info = sh)
  }
})

test_that("denoise_signal_offline emite warning quando method é usado", {
  set.seed(3)
  x = rnorm(64)
  sch = lifting_scheme("haar")
  expect_warning(
    denoise_signal_offline(x, sch, levels = 2, method = "soft"),
    "deprecated|shrinkage"
  )
})

test_that("denoise_signal_offline: passar method e shrinkage juntos é erro", {
  set.seed(4)
  x = rnorm(64)
  sch = lifting_scheme("haar")
  expect_error(
    denoise_signal_offline(x, sch, method = "soft", shrinkage = "hard"),
    "method.*shrinkage|both"
  )
})

test_that("denoise_signal_offline: threshold_method inválido é erro", {
  x = rnorm(64); sch = lifting_scheme("haar")
  expect_error(
    denoise_signal_offline(x, sch, threshold_method = "nonsense"),
    "threshold_method"
  )
})

test_that("denoise_signal_offline: shrinkage inválido é erro", {
  x = rnorm(64); sch = lifting_scheme("haar")
  expect_error(
    denoise_signal_offline(x, sch, shrinkage = "nonsense"),
    "shrinkage"
  )
})

test_that("denoise_signal_causal aceita shrinkage e threshold_method", {
  set.seed(5)
  x = rnorm(256)
  sch = lifting_scheme("haar")

  expect_no_error(
    denoise_signal_causal(x, sch, levels = 2, window_size = 64,
                          threshold_method = "universal",
                          shrinkage = "soft")
  )
})

test_that("denoise_signal_causal: shrinkage equivalente a method (deprecated)", {
  set.seed(6)
  x = rnorm(256)
  sch = lifting_scheme("haar")
  for (sh in c("hard", "soft", "semisoft")) {
    new_out = denoise_signal_causal(x, sch, levels = 2, window_size = 64,
                                    shrinkage = sh)
    old_out = suppressWarnings(
      denoise_signal_causal(x, sch, levels = 2, window_size = 64, method = sh)
    )
    expect_equal(new_out, old_out, info = sh)
  }
})

test_that("new_wavelet_stream aceita shrinkage e threshold_method", {
  sch = lifting_scheme("haar")
  expect_no_error(
    new_wavelet_stream(sch, window_size = 64, levels = 2,
                       threshold_method = "universal",
                       shrinkage = "soft")
  )
})

test_that("new_wavelet_stream: shrinkage equivalente a method", {
  set.seed(7)
  x = rnorm(256)
  sch = lifting_scheme("haar")
  for (sh in c("hard", "soft", "semisoft")) {
    p_new = new_wavelet_stream(sch, window_size = 64, levels = 2,
                               shrinkage = sh)
    p_old = suppressWarnings(
      new_wavelet_stream(sch, window_size = 64, levels = 2, method = sh)
    )
    out_new = vapply(x, p_new, numeric(1))
    out_old = vapply(x, p_old, numeric(1))
    expect_equal(out_new, out_old, info = sh)
  }
})
