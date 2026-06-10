suppressPackageStartupMessages({
  library(parallel)
  library(microbenchmark)
  library(adlift)
})

N_PTS = 1024L

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
    steps = abs(rnorm(N_PTS - 1L, mean = 1, sd = 0.4))
  }
  cumsum(c(0, steps))
}

gen_pure = function(signal_name, t_phys) {
  n = length(t_phys)
  t_n = (t_phys - min(t_phys)) / (max(t_phys) - min(t_phys))
  if (signal_name == "linear_phys")
    return(0.5 * t_phys)
  if (signal_name == "trend_events") {
    trend = 0.5 * t_phys
    events = Reduce("+", lapply(c(0.25, 0.55, 0.80), function(p)
      3 * exp(-((t_phys - p * max(t_phys)) / 5)^2)))
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
  stop("Unknown signal: ", signal_name)
}

NOISE_SD = c(linear_phys = 0.15, trend_events = 0.15, blocks_gapped = 0.50,
             blocks_dj_irr = 0.30, bumps_dj_irr = 0.30,
             doppler_dj_irr = 0.30, heavisine_dj_irr = 0.30)

PRED_FNS = list(LinearPred = adlift::LinearPred, QuadPred = adlift::QuadPred,
                CubicPred = adlift::CubicPred, AdaptPred = adlift::AdaptPred)

CONFIGS = do.call(rbind, lapply(names(PRED_FNS), function(pname) {
  expand.grid(pred_name = pname,
              neighbours = c(1L, 2L, 3L),
              int = c(TRUE, FALSE),
              clo = c(FALSE, TRUE),
              rule = c("median", "mean"),
              stringsAsFactors = FALSE)
}))
CONFIGS$label = with(CONFIGS, paste(
  tolower(sub("Pred", "", pred_name)), paste0("n", neighbours),
  ifelse(int, "int", "noint"), ifelse(clo, "clo", "noclo"), rule, sep = "_"))

DJ_SIGNALS = c("linear_phys", "trend_events", "blocks_gapped",
               "blocks_dj_irr", "bumps_dj_irr",
               "doppler_dj_irr", "heavisine_dj_irr")
N_SIM = 1000L

col_stats = function(x, prefix) {
  x = x[!is.na(x)]
  s = c(min = min(x), q1 = unname(quantile(x, .25)), median = median(x),
        mean = mean(x), q3 = unname(quantile(x, .75)), max = max(x),
        se = sd(x)/sqrt(length(x)))
  setNames(as.list(s),
    paste0(prefix, c("min","q1","median","mean","q3","max","se")))
}

TMPDIR = "data/tmp_adlift_irr"
dir.create(TMPDIR, showWarnings = FALSE, recursive = TRUE)

jobs = list()
for (sig in DJ_SIGNALS) {
  for (ci in seq_len(nrow(CONFIGS))) {
    cfg = CONFIGS[ci, ]
    key = sprintf("%s__%s", sig, cfg$label)
    fout = file.path(TMPDIR, paste0(key, ".rds"))
    if (!file.exists(fout))
      jobs[[length(jobs) + 1]] = list(signal = sig, cfg = cfg, fout = fout)
  }
}

total = nrow(CONFIGS) * length(DJ_SIGNALS)
message(sprintf("Configs: %d | Signals: %d | Total: %d | Jobs: %d",
  nrow(CONFIGS), length(DJ_SIGNALS), total, length(jobs)))

run_job = function(job) {
  sig = job$signal
  cfg = job$cfg
  t_phys = make_t_phys(sig)
  pure = gen_pure(sig, t_phys)
  noise_sd = NOISE_SD[[sig]]
  pred_fn = PRED_FNS[[cfg$pred_name]]
  neigh = as.integer(cfg$neighbours)

  mse_v = numeric(N_SIM)
  time_v = numeric(N_SIM)

  for (i in seq_len(N_SIM)) {
    set.seed(i)
    noisy = pure + rnorm(length(pure), sd = noise_sd)
    res = tryCatch({
      mb = microbenchmark(
        out <- as.vector(adlift::denoise(t_phys, noisy, pred = pred_fn,
          neigh = neigh, int = cfg$int, clo = cfg$clo,
          keep = 2L, rule = cfg$rule)),
        times = 3L)
      list(mse = mean((pure - out)^2), time = median(mb$time)/1e9)
    }, error = function(e) NULL)
    mse_v[i] = if (!is.null(res)) res$mse else NA_real_
    time_v[i] = if (!is.null(res)) res$time else NA_real_
  }

  result = as.data.frame(c(
    list(Signal = sig,
         Pkg = "adlift",
         Wavelet = cfg$pred_name,
         Boundary = cfg$label,
         NoiseSd = noise_sd,
         N = sum(!is.na(mse_v))),
    col_stats(mse_v, "MSE_"),
    col_stats(time_v, "Time_")
  ), stringsAsFactors = FALSE)

  saveRDS(result, job$fout)
  result
}

if (length(jobs) > 0) {
  RNGkind("L'Ecuyer-CMRG"); set.seed(2025)
  n_cores = max(1L, detectCores() - 1L)
  t0 = proc.time()
  mclapply(jobs, run_job, mc.cores = n_cores,
           mc.preschedule = FALSE, mc.set.seed = TRUE)
  message(sprintf("Done in %.1f min", (proc.time()-t0)[["elapsed"]]/60))
}

files = list.files(TMPDIR, pattern = "\\.rds$", full.names = TRUE)
benchmark_adlift_irregular = do.call(rbind, lapply(files, readRDS))
benchmark_adlift_irregular = benchmark_adlift_irregular[
  order(benchmark_adlift_irregular$Signal,
        benchmark_adlift_irregular$Wavelet,
        benchmark_adlift_irregular$Boundary), ]
rownames(benchmark_adlift_irregular) = NULL
save(benchmark_adlift_irregular,
     file = "data/benchmark_adlift_irregular.rda", compress = "xz")
message("Saved data/benchmark_adlift_irregular.rda")
