
test_that(
  "one_sided preserves perfect reconstruction for all built-in wavelets",
  {
    wavelets = c("haar", "db2", "cdf53", "cdf97", "dd4")
    signals = c("random", "ramp", "sine", "doppler", "heavisine")
    for (wname in wavelets) {
      sch = lifting_scheme(wname)
      for (sig in signals) {
        x = rLifting:::.generate_signal(sig, n = 64)
        res = lwt(x, sch, levels = 2, extension = "one_sided")
        err = max(abs(x - ilwt(res)))
        expect_lt(err, 1e-9,
          label = sprintf("PR failed: %s + %s, err=%.2e", wname, sig, err))
      }
    }
  }
)

test_that("one_sided is accepted by all public functions", {
  sch = lifting_scheme("haar")
  x = rnorm(128)

  res = denoise_signal_offline(x, sch, levels = 2, extension = "one_sided")
  expect_equal(length(res), 128)
  expect_false(any(is.na(res)))

  res = denoise_signal_causal(x, sch, levels = 2, window_size = 32,
    extension = "one_sided")
  expect_equal(length(res), 128)
  expect_false(any(is.na(res)))

  proc = new_wavelet_stream(sch, window_size = 32, levels = 2,
    extension = "one_sided")
  out = numeric(128)
  for (i in seq_along(x)) out[i] = proc(x[i])
  expect_equal(length(out), 128)
  expect_false(any(is.na(out)))
})

test_that("one_sided produces different boundary behaviour than symmetric", {
  sch = lifting_scheme("cdf97")
  set.seed(7)
  x = cumsum(rnorm(64))

  r_sym = lwt(x, sch, levels = 2, extension = "symmetric")
  r_os = lwt(x, sch, levels = 2, extension = "one_sided")
  expect_false(identical(r_sym$coeffs$d1, r_os$coeffs$d1))

  set.seed(7); x2 = cumsum(rnorm(128))
  off_sym = denoise_signal_offline(x2, sch, levels = 2,
    extension = "symmetric")
  off_os = denoise_signal_offline(x2, sch, levels = 2,
    extension = "one_sided")
  expect_false(identical(off_sym, off_os))

  set.seed(7); x3 = cumsum(rnorm(128))
  caus_sym = denoise_signal_causal(x3, sch, levels = 2, window_size = 32,
    extension = "symmetric")
  caus_os = denoise_signal_causal(x3, sch, levels = 2, window_size = 32,
    extension = "one_sided")
  expect_false(identical(caus_sym, caus_os))
})

test_that(
  "one_sided is more noise-robust at boundary than 2-point local_linear",
  {
    set.seed(123)
    sch = lifting_scheme("cdf97")
    n = 128
    x_clean = seq(0, 1, length.out = n)

    x_noisy = x_clean
    x_noisy[2] = x_noisy[2] + 100

    res_ll = denoise_signal_offline(x_noisy, sch, levels = 2,
      extension = "local_linear", ll_k = 2L)
    res_os = denoise_signal_offline(x_noisy, sch, levels = 2,
      extension = "one_sided")

    mse = function(a, b) mean((a - b)^2)
    truth = denoise_signal_offline(x_clean, sch, levels = 2,
      extension = "one_sided")
    expect_lt(mse(res_os[1:5], truth[1:5]),
              mse(res_ll[1:5], truth[1:5]))
  }
)

test_that("one_sided causal mode has no look-ahead leakage", {
  sch = lifting_scheme("haar")
  set.seed(1)
  x = rnorm(64)

  p1 = new_wavelet_stream(sch, window_size = 16, extension = "one_sided")
  p2 = new_wavelet_stream(sch, window_size = 16, extension = "one_sided")

  out1 = numeric(32); out2 = numeric(32)
  for (i in 1:32) { out1[i] = p1(x[i]); out2[i] = p2(x[i]) }
  for (i in 33:64) p2(rnorm(1))

  expect_equal(out1, out2)
})

test_that("one_sided warns when t is supplied (irregular grid)", {
  sch = lifting_scheme("cdf53")
  n = 64
  x = rnorm(n)
  t = cumsum(c(0, abs(rnorm(n - 1, mean = 1, sd = 0.4))))

  expect_warning(lwt(x, sch, extension = "one_sided", t = t),
    "one_sided")
  expect_warning(denoise_signal_offline(x, sch, extension = "one_sided",
    t = t), "one_sided")
  expect_warning(denoise_signal_causal(x, sch, extension = "one_sided",
    t = t, window_size = 31), "one_sided")
  expect_warning(new_wavelet_stream(sch, extension = "one_sided",
    irregular = TRUE), "one_sided")
})

test_that("one_sided does not warn on regular grid (no t)", {
  sch = lifting_scheme("cdf53")
  x = rnorm(64)

  expect_no_warning(lwt(x, sch, extension = "one_sided"))
  expect_no_warning(denoise_signal_offline(x, sch, extension = "one_sided"))
  expect_no_warning(new_wavelet_stream(sch, extension = "one_sided",
    irregular = FALSE))
})

test_that("one_sided overhead vs symmetric is < 15% for N=1024", {
  skip_if_not_installed("microbenchmark")
  library(microbenchmark)

  sch = lifting_scheme("cdf97")
  x = rnorm(1024)
  lvls = 4

  mb = microbenchmark(
    sym = denoise_signal_offline(x, sch, levels = lvls,
      extension = "symmetric"),
    os = denoise_signal_offline(x, sch, levels = lvls,
      extension = "one_sided"),
    times = 50L
  )
  t_sym = median(mb$time[mb$expr == "sym"])
  t_os = median(mb$time[mb$expr == "os"])

  overhead_pct = (t_os / t_sym - 1) * 100
  expect_lt(overhead_pct, 15,
    label = sprintf("overhead %.1f%% exceeds 15%%", overhead_pct))
})
