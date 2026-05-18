
test_that(".generate_signal produces all D-J signals with correct length", {
  dj_signals = c("doppler", "heavisine", "bumps", "blocks")
  for (sig in dj_signals) {
    x = rLifting:::.generate_signal(sig, n = 256)
    expect_equal(length(x), 256)
    expect_false(any(is.na(x)))
  }
  x_blocks = rLifting:::.generate_signal("blocks", n = 1024)
  diffs = diff(x_blocks)
  expect_gt(mean(diffs == 0), 0.95)
})

test_that("compute_adaptive_threshold returns correct list structure", {
  mock_lwt = list(
    coeffs = list(
      d1 = rnorm(128),
      d2 = rnorm(64),
      d3 = rnorm(32),
      a3 = rnorm(32)
    )
  )
  class(mock_lwt) = "lwt"

  res = compute_adaptive_threshold(mock_lwt, alpha = 0.5, beta = 1.0)

  expect_type(res, "list")
  expect_named(res, c("d1", "d2", "d3"))
  expect_true(is.numeric(res$d1))
  expect_true(length(res$d1) == 1)
})

test_that("validate_compact_support identifies FIR constraints", {
  sch = lifting_scheme("haar")
  res = validate_compact_support(sch, max_width = 2)
  expect_true(res$passed)
  expect_equal(res$metric, 2)
})

test_that("diagnose_wavelet runs without error (smoke test)", {
  sch = lifting_scheme("db2")
  config = list(is_ortho = TRUE, vm_degrees = c(0, 1), max_taps = 4)

  expect_silent(diagnose_wavelet(sch, config, verbose = FALSE))

  expect_output(diagnose_wavelet(sch, config, verbose = TRUE), "DIAGNOSIS")
})

test_that('S3 Plot methods run without error', {
  sch = lifting_scheme('haar')
  tmp = tempfile()
  pdf(file = tmp)
  on.exit({dev.off(); unlink(tmp)}, add = TRUE)

  expect_error(plot(sch), NA)

  x = rnorm(64)
  res = lwt(x, sch, levels = 2)
  expect_error(plot(res), NA)

})
