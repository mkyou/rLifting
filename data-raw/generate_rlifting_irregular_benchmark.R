# Irregular-grid benchmark for rLifting.
#
# Sweep: 12 method configs x 6 wavelets x 5 boundaries x 7 signals x 3 modes
#   = 7560 jobs
# Methods mirror the regular v2 grid (universal/sure x hard/soft/semisoft/scad,
# default and tuned alpha/beta where applicable).
# Signals: 3 physically motivated (linear_phys, trend_events, blocks_gapped)
# plus 4 Donoho-Johnstone classics re-sampled in irregular grids
# (blocks_dj_irr, bumps_dj_irr, doppler_dj_irr, heavisine_dj_irr).
#
# Offline mode produces both position-aware ("pos", t = t_phys passed in)
# and uniform-grid-ignoring ("fix", t = NULL) MSE/time to expose the value
# of irregular handling. Causal/stream produce only position-aware results.

suppressPackageStartupMessages({
  library(parallel)
  library(microbenchmark)
  library(rLifting)
})

N_SIM = 1000L
N_PTS = 1024L
LEVELS_OFF = 4L
LEVELS_CAUSAL = 3L
WINDOW_SIZE = 255L
WARMUP = WINDOW_SIZE - 1L

DJ_PHYS = c("linear_phys", "trend_events", "blocks_gapped")
DJ_IRR  = c("blocks_dj_irr", "bumps_dj_irr",
            "doppler_dj_irr", "heavisine_dj_irr")
DJ_SIGNALS = c(DJ_PHYS, DJ_IRR)

NOISE_SD = c(linear_phys = 0.15, trend_events = 0.15, blocks_gapped = 0.50,
             blocks_dj_irr = 0.30, bumps_dj_irr = 0.30,
             doppler_dj_irr = 0.30, heavisine_dj_irr = 0.30)

WAVELETS   = c("haar", "db2", "cdf53", "cdf97", "dd4", "lazy")
BOUNDARIES = c("symmetric", "periodic", "zero", "local_linear", "one_sided")
MODES      = c("offline", "causal", "stream")

NEW_METHODS = list(
  list(label = "universal_hard",          threshold = "universal", shrinkage = "hard",     tuned = FALSE),
  list(label = "universal_soft",          threshold = "universal", shrinkage = "soft",     tuned = FALSE),
  list(label = "universal_semisoft",      threshold = "universal", shrinkage = "semisoft", tuned = FALSE),
  list(label = "universal_scad",          threshold = "universal", shrinkage = "scad",     tuned = FALSE),
  list(label = "universal_tuned_hard",    threshold = "universal", shrinkage = "hard",     tuned = TRUE),
  list(label = "universal_tuned_soft",    threshold = "universal", shrinkage = "soft",     tuned = TRUE),
  list(label = "universal_tuned_semisoft",threshold = "universal", shrinkage = "semisoft", tuned = TRUE),
  list(label = "universal_tuned_scad",    threshold = "universal", shrinkage = "scad",     tuned = TRUE),
  list(label = "sure_hard",               threshold = "sure",      shrinkage = "hard",     tuned = FALSE),
  list(label = "sure_soft",               threshold = "sure",      shrinkage = "soft",     tuned = FALSE),
  list(label = "sure_semisoft",           threshold = "sure",      shrinkage = "semisoft", tuned = FALSE),
  list(label = "sure_scad",               threshold = "sure",      shrinkage = "scad",     tuned = FALSE)
)

TMPDIR   = "data/tmp_rlifting_irr"
CACHEDIR = "data/tmp_rlifting_irr_v2"
dir.create(TMPDIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(CACHEDIR, showWarnings = FALSE, recursive = TRUE)

make_t_phys = function(signal_name) {
  seeds = c(linear_phys = 101L, trend_events = 102L, blocks_gapped = 103L,
            blocks_dj_irr = 201L, bumps_dj_irr = 202L,
            doppler_dj_irr = 203L, heavisine_dj_irr = 204L)
  set.seed(seeds[[signal_name]])
  if (signal_name %in% c("linear_phys", "trend_events")) {
    steps = abs(rnorm(N_PTS - 1L, mean = 1, sd = 0.9))
  } else if (signal_name == "blocks_gapped") {
    steps = abs(rnorm(N_PTS - 1L, mean = 0.5, sd = 0.3))
    big = sample.int(N_PTS - 1L, 12L)
    steps[big] = runif(12L, min = 6, max = 18)
  } else {
    # DJ-irregular: moderate jitter on top of uniform spacing
    steps = abs(rnorm(N_PTS - 1L, mean = 1, sd = 0.4))
  }
  cumsum(c(0, steps))
}

gen_pure = function(signal_name, t_phys) {
  n = length(t_phys)
  t_n = (t_phys - min(t_phys)) / (max(t_phys) - min(t_phys))
  if (signal_name == "linear_phys") return(0.5 * t_phys)
  if (signal_name == "trend_events") {
    trend = 0.5 * t_phys
    events = Reduce("+", lapply(c(0.25, 0.55, 0.80), function(p) {
      3 * exp(-((t_phys - p * max(t_phys)) / 5)^2)
    }))
    return(trend + events)
  }
  if (signal_name == "blocks_gapped") {
    pos = c(0.12, 0.28, 0.42, 0.58, 0.72, 0.87)
    h = c(5, -7, 4, -5, 6, -4)
    x = numeric(n)
    for (j in seq_along(pos)) x = x + h[j] * (1 + sign(t_n - pos[j])) / 2
    return(x)
  }
  if (signal_name == "doppler_dj_irr") {
    eps = 0.05
    return(sqrt(t_n * (1 - t_n)) * sin((2 * pi * 1.05) / (t_n + eps)))
  }
  if (signal_name == "heavisine_dj_irr") {
    return(4 * sin(4 * pi * t_n) - sign(t_n - 0.3) - sign(0.72 - t_n))
  }
  if (signal_name == "blocks_dj_irr") {
    pos = c(0.1, 0.13, 0.15, 0.23, 0.25, 0.40, 0.44, 0.65, 0.76, 0.78, 0.81)
    h = c(4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 5.1, -4.2)
    x = numeric(n)
    for (j in seq_along(pos)) x = x + h[j] * (1 + sign(t_n - pos[j])) / 2
    return(x)
  }
  if (signal_name == "bumps_dj_irr") {
    pos = c(0.1, 0.13, 0.15, 0.23, 0.25, 0.40, 0.44, 0.65, 0.76, 0.78, 0.81)
    h = c(4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 5.1, -4.2)
    w = c(0.005, 0.005, 0.006, 0.01, 0.01, 0.03, 0.01, 0.01, 0.005, 0.008, 0.005)
    x = numeric(n)
    for (j in seq_along(pos)) {
      x = x + h[j] * (1 + abs((t_n - pos[j]) / w[j])^4)^(-1)
    }
    return(x)
  }
  stop(sprintf("Unknown signal: '%s'", signal_name))
}

col_stats = function(x, prefix) {
  x = x[!is.na(x)]
  if (length(x) == 0) {
    s = rep(NA_real_, 7)
  } else {
    s = c(min(x), quantile(x, 0.25, names = FALSE), median(x),
          mean(x), quantile(x, 0.75, names = FALSE), max(x),
          sd(x) / sqrt(length(x)))
  }
  setNames(as.list(s),
           paste0(prefix, c("min", "q1", "median", "mean", "q3", "max", "se")))
}

# Tune alpha/beta per (signal, wavelet) on a seeded clean realisation.
cache_path = file.path(CACHEDIR, "alpha_beta_cache.rds")
ab_cache = if (file.exists(cache_path)) readRDS(cache_path) else list()
for (sig in DJ_SIGNALS) for (wav in WAVELETS) {
  key = paste(sig, wav, sep = "__")
  if (!is.null(ab_cache[[key]])) next
  set.seed(2025L + sum(utf8ToInt(key)))
  t_phys = make_t_phys(sig)
  pure   = gen_pure(sig, t_phys)
  noisy  = pure + rnorm(N_PTS, sd = NOISE_SD[[sig]])
  sch    = lifting_scheme(wav)
  tuned  = tune_alpha_beta(noisy, sch, levels = LEVELS_OFF)
  ab_cache[[key]] = list(alpha = tuned$alpha, beta = tuned$beta)
  message(sprintf("Tuned %s: alpha=%.3f beta=%.3f", key,
                  tuned$alpha, tuned$beta))
}
saveRDS(ab_cache, cache_path)

CONFIGS = expand.grid(
  Signal      = DJ_SIGNALS,
  Wavelet     = WAVELETS,
  Boundary    = BOUNDARIES,
  MethodLabel = vapply(NEW_METHODS, `[[`, character(1), "label"),
  Mode        = MODES,
  stringsAsFactors = FALSE
)
method_lookup = setNames(NEW_METHODS,
                         vapply(NEW_METHODS, `[[`, character(1), "label"))

jobs = vector("list", nrow(CONFIGS))
for (i in seq_len(nrow(CONFIGS))) {
  cfg = CONFIGS[i, ]
  m   = method_lookup[[cfg$MethodLabel]]
  key = paste(cfg$Signal, cfg$Wavelet, cfg$Boundary,
              cfg$MethodLabel, cfg$Mode, sep = "__")
  fout = file.path(TMPDIR, paste0(key, ".rds"))
  jobs[[i]] = list(cfg = cfg, m = m, fout = fout,
                   done = file.exists(fout))
}

pending = jobs[!vapply(jobs, `[[`, logical(1), "done")]
message(sprintf(
  "Configs: %d | Pending: %d | Skipping: %d already done",
  length(jobs), length(pending), length(jobs) - length(pending)
))

run_job = function(job) {
  cfg = job$cfg
  m   = job$m
  sig = cfg$Signal; wav = cfg$Wavelet
  bnd = cfg$Boundary; mode = cfg$Mode

  sch    = lifting_scheme(wav)
  t_phys = make_t_phys(sig)
  pure   = gen_pure(sig, t_phys)
  noise_sd = NOISE_SD[[sig]]
  ab_key = paste(sig, wav, sep = "__")

  alpha = if (m$tuned) ab_cache[[ab_key]]$alpha else 0.3
  beta  = if (m$tuned) ab_cache[[ab_key]]$beta  else 1.2

  mse_pos_v  = numeric(N_SIM)
  mse_fix_v  = numeric(N_SIM)
  time_pos_v = numeric(N_SIM)
  time_fix_v = numeric(N_SIM)

  for (i in seq_len(N_SIM)) {
    set.seed(i)
    noisy = pure + rnorm(N_PTS, sd = noise_sd)

    tryCatch({
      if (mode == "offline") {
        mb_pos = microbenchmark(
          out_pos <- suppressWarnings(
            denoise_signal_offline(noisy, sch, levels = LEVELS_OFF,
                                   alpha = alpha, beta = beta,
                                   threshold_method = m$threshold,
                                   shrinkage = m$shrinkage,
                                   extension = bnd, t = t_phys)
          ), times = 3L)
        mb_fix = microbenchmark(
          out_fix <- suppressWarnings(
            denoise_signal_offline(noisy, sch, levels = LEVELS_OFF,
                                   alpha = alpha, beta = beta,
                                   threshold_method = m$threshold,
                                   shrinkage = m$shrinkage,
                                   extension = bnd)
          ), times = 3L)
        mse_pos_v[i]  = mean((pure - out_pos)^2)
        mse_fix_v[i]  = mean((pure - out_fix)^2)
        time_pos_v[i] = median(mb_pos$time) / 1e9
        time_fix_v[i] = median(mb_fix$time) / 1e9

      } else if (mode == "causal") {
        mb = microbenchmark(
          out_pos <- suppressWarnings(
            denoise_signal_causal(noisy, sch, levels = LEVELS_CAUSAL,
                                  window_size = WINDOW_SIZE,
                                  alpha = alpha, beta = beta,
                                  threshold_method = m$threshold,
                                  shrinkage = m$shrinkage,
                                  extension = bnd, t = t_phys)
          ), times = 1L)
        mse_pos_v[i]  = mean((pure - out_pos)^2)
        mse_fix_v[i]  = NA_real_
        time_pos_v[i] = mb$time[1] / 1e9
        time_fix_v[i] = NA_real_

      } else {  # stream
        proc = suppressWarnings(
          new_wavelet_stream(sch, window_size = WINDOW_SIZE,
                             levels = LEVELS_CAUSAL,
                             alpha = alpha, beta = beta,
                             threshold_method = m$threshold,
                             shrinkage = m$shrinkage,
                             extension = bnd, update_freq = 1L,
                             irregular = TRUE)
        )
        out_pos = numeric(N_PTS)
        t0 = as.numeric(Sys.time())
        for (j in seq_len(N_PTS)) out_pos[j] = proc(noisy[j], t_phys[j])
        t_total = as.numeric(Sys.time()) - t0
        mse_pos_v[i]  = mean((pure - out_pos)^2)
        mse_fix_v[i]  = NA_real_
        time_pos_v[i] = t_total
        time_fix_v[i] = NA_real_
      }
    }, error = function(e) {
      message(sprintf("ERROR in %s__%s__%s__%s__%s: %s",
                      sig, wav, bnd, m$label, mode, conditionMessage(e)))
      mse_pos_v[i]  <<- NA_real_
      mse_fix_v[i]  <<- NA_real_
      time_pos_v[i] <<- NA_real_
      time_fix_v[i] <<- NA_real_
    })
  }

  ratio = mse_fix_v / mse_pos_v
  result = as.data.frame(c(
    list(Signal = sig,
         Pkg = "rLifting",
         Mode = mode,
         Wavelet = wav,
         Boundary = bnd,
         Method = m$label,
         ThresholdMethod = m$threshold,
         Shrinkage = m$shrinkage,
         AlphaUsed = alpha,
         BetaUsed = beta,
         NoiseSd = noise_sd,
         N = sum(!is.na(mse_pos_v))),
    col_stats(mse_pos_v,  "MSEpos_"),
    col_stats(mse_fix_v,  "MSEfix_"),
    col_stats(ratio,      "Ratio_"),
    col_stats(time_pos_v, "Timepos_"),
    col_stats(time_fix_v, "Timefix_")
  ), stringsAsFactors = FALSE)

  saveRDS(result, job$fout)
  result
}

if (length(pending) > 0) {
  RNGkind("L'Ecuyer-CMRG"); set.seed(2026)
  n_cores = max(1L, detectCores() - 1L)
  message(sprintf("Running on %d cores...", n_cores))
  t0 = proc.time()
  mclapply(pending, run_job, mc.cores = n_cores,
           mc.preschedule = FALSE, mc.set.seed = TRUE)
  elapsed = (proc.time() - t0)[["elapsed"]]
  message(sprintf("Done in %.1f s  (%.1f min)", elapsed, elapsed / 60))
}

message("Combining checkpoints...")
files = list.files(TMPDIR, pattern = "\\.rds$", full.names = TRUE)
benchmark_rlifting_irregular = do.call(rbind, lapply(files, readRDS))
benchmark_rlifting_irregular = benchmark_rlifting_irregular[
  order(benchmark_rlifting_irregular$Signal,
        benchmark_rlifting_irregular$Mode,
        benchmark_rlifting_irregular$Wavelet,
        benchmark_rlifting_irregular$Boundary,
        benchmark_rlifting_irregular$Method), ]
rownames(benchmark_rlifting_irregular) = NULL
message(sprintf("Final: %d rows x %d cols",
                nrow(benchmark_rlifting_irregular),
                ncol(benchmark_rlifting_irregular)))
save(benchmark_rlifting_irregular,
     file = "data/benchmark_rlifting_irregular.rda", compress = "xz")
message("Saved data/benchmark_rlifting_irregular.rda")
