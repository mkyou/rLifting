
# Generates data/benchmark_rlifting_irregular.rda
#
# Purpose: quantify the MSE advantage of position-aware lifting (t=t_phys)
#          versus fixed-coefficient lifting on irregularly-sampled signals.
#          wavethresh excluded: no native irregular-grid support.
#
# Signals (physical time, NOT normalised to [0,1]):
#   "linear_phys"    — pure linear trend: y = 0.5 * t_phys
#   "trend_events"   — linear trend + 3 Gaussian events
#   "blocks_gapped"  — blocks with 12 injected large gaps near jumps
#
# Noise is ABSOLUTE (not relative to signal range) so the position-induced
# residual (~0.3) stays comparable to noise level and the advantage is visible.
#   linear_phys / trend_events : noise_sd = 0.15
#   blocks_gapped               : noise_sd = 0.50
#
# The irregular grid (t_phys) is FIXED per signal (seed = 101/102/103).
# Only noise varies across replications.
#
# Configs: same 90 as generate_rlifting_benchmark.R
#   Wavelets   : haar, db2, cdf53, cdf97, dd4, lazy    (6)
#   Extensions : symmetric, periodic, zero,
#                local_linear, one_sided                (5)
#   Methods    : hard, soft, semisoft                  (3)
#   Total      : 6 x 5 x 3 = 90
#
# Each row reports offline posicional/fixo AND causal posicional/fixo:
#   MSEpos_*, MSEfix_*    — offline stats over reps
#   Ratio_*               — MSEfix / MSEpos offline (>1 means posicional wins)
#   Timepos_*, Timefix_*  — offline timing
#   MSEcpos_*, MSEcfix_*  — causal stats (warmup excluded)
#   CRatio_*              — MSEcfix / MSEcpos causal
#   Timecpos_*, Timecfix_ — causal timing
# Causal: window_size=64, levels=3, warmup=63 samples excluded from MSE
#
# Reps: 1000 | n: 256 | levels: 4
# Checkpoint/resume: data/tmp_rlifting_irr/
# Runtime estimate: ~15 min on 11 cores
#
# source("data-raw/generate_rlifting_irregular_benchmark.R")

suppressPackageStartupMessages({
  library(parallel)
  library(microbenchmark)
  library(rLifting)
})

# ── Canonical irregular grids (shared with adlift / nlt scripts) ──────────────

make_t_phys = function(signal_name) {
  seeds = c(linear_phys=101L, trend_events=102L, blocks_gapped=103L)
  seed  = seeds[[signal_name]]
  n     = 256L
  set.seed(seed)
  if (signal_name %in% c("linear_phys", "trend_events")) {
    steps = abs(rnorm(n - 1L, mean = 1, sd = 0.9))
  } else {                         # blocks_gapped: mostly fine + 12 big gaps
    steps = abs(rnorm(n - 1L, mean = 0.5, sd = 0.3))
    big   = sample.int(n - 1L, 12L)
    steps[big] = runif(12L, min = 6, max = 18)
  }
  cumsum(c(0, steps))
}

gen_pure = function(signal_name, t_phys) {
  n   = length(t_phys)
  t_n = (t_phys - min(t_phys)) / (max(t_phys) - min(t_phys))
  if (signal_name == "linear_phys") {
    return(0.5 * t_phys)
  }
  if (signal_name == "trend_events") {
    trend  = 0.5 * t_phys
    events = Reduce("+", lapply(c(0.25, 0.55, 0.80), function(p) {
      3 * exp(-((t_phys - p * max(t_phys)) / 5)^2)
    }))
    return(trend + events)
  }
  if (signal_name == "blocks_gapped") {
    pos = c(0.12, 0.28, 0.42, 0.58, 0.72, 0.87)
    h   = c(5, -7, 4, -5, 6, -4)
    x   = numeric(n)
    for (j in seq_along(pos)) x = x + h[j] * (1 + sign(t_n - pos[j])) / 2
    return(x)
  }
  stop("Unknown signal: ", signal_name)
}

NOISE_SD = c(linear_phys = 0.15, trend_events = 0.15, blocks_gapped = 0.50)

# ── Config grid (same as regular benchmark) ───────────────────────────────────

CONFIGS = expand.grid(
  wavelet   = c("haar", "db2", "cdf53", "cdf97", "dd4", "lazy"),
  extension = c("symmetric", "periodic", "zero", "local_linear", "one_sided"),
  method    = c("hard", "soft", "semisoft"),
  stringsAsFactors = FALSE
)

DJ_SIGNALS = c("linear_phys", "trend_events", "blocks_gapped")
N_SIM         = 1000L
N_PTS         = 256L
LEVELS        = 4L
LEVELS_CAUSAL = 3L
WINDOW_SIZE   = 64L
WARMUP        = WINDOW_SIZE - 1L
ALPHA         = 0.3
BETA          = 1.2

# ── Helpers ───────────────────────────────────────────────────────────────────

col_stats = function(x, prefix) {
  x = x[!is.na(x)]
  s = c(min    = min(x),
        q1     = unname(quantile(x, 0.25)),
        median = median(x),
        mean   = mean(x),
        q3     = unname(quantile(x, 0.75)),
        max    = max(x),
        se     = sd(x) / sqrt(length(x)))
  setNames(as.list(s), paste0(prefix, c("min","q1","median","mean","q3","max","se")))
}

# ── Checkpoint setup ──────────────────────────────────────────────────────────

TMPDIR = "data/tmp_rlifting_irr"
dir.create(TMPDIR, showWarnings = FALSE, recursive = TRUE)

jobs = list()
for (sig in DJ_SIGNALS) {
  for (ci in seq_len(nrow(CONFIGS))) {
    cfg  = CONFIGS[ci, ]
    key  = sprintf("%s__%s__%s__%s", sig, cfg$wavelet, cfg$extension, cfg$method)
    fout = file.path(TMPDIR, paste0(key, ".rds"))
    if (!file.exists(fout))
      jobs[[length(jobs) + 1]] = list(signal=sig, cfg=cfg, ci=ci, fout=fout)
  }
}

total = nrow(CONFIGS) * length(DJ_SIGNALS)
message(sprintf("Configs: %d | Signals: %d | Total rows: %d",
                nrow(CONFIGS), length(DJ_SIGNALS), total))
message(sprintf("Jobs to run: %d  (skipping: %d already done)",
                length(jobs), total - length(jobs)))

# ── Worker ────────────────────────────────────────────────────────────────────

run_job = function(job) {
  sig      = job$signal
  cfg      = job$cfg
  sch      = lifting_scheme(cfg$wavelet)
  t_phys   = make_t_phys(sig)
  pure     = gen_pure(sig, t_phys)
  noise_sd = NOISE_SD[[sig]]

  mse_pos  = numeric(N_SIM)
  mse_fix  = numeric(N_SIM)
  time_pos = numeric(N_SIM)
  time_fix = numeric(N_SIM)

  for (i in seq_len(N_SIM)) {
    set.seed(i)
    noisy = pure + rnorm(N_PTS, sd = noise_sd)

    res = tryCatch({
      mb_pos = microbenchmark(
        out_pos <- suppressWarnings(
          denoise_signal_offline(noisy, sch, levels=LEVELS,
                                 alpha=ALPHA, beta=BETA,
                                 method=cfg$method, extension=cfg$extension,
                                 t=t_phys)),
        times = 3L)
      mb_fix = microbenchmark(
        out_fix <- denoise_signal_offline(noisy, sch, levels=LEVELS,
                                          alpha=ALPHA, beta=BETA,
                                          method=cfg$method, extension=cfg$extension),
        times = 3L)
      list(
        mse_pos  = mean((pure - out_pos)^2),
        mse_fix  = mean((pure - out_fix)^2),
        time_pos = median(mb_pos$time) / 1e9,
        time_fix = median(mb_fix$time) / 1e9
      )
    }, error = function(e) NULL)

    if (!is.null(res)) {
      mse_pos[i]  = res$mse_pos
      mse_fix[i]  = res$mse_fix
      time_pos[i] = res$time_pos
      time_fix[i] = res$time_fix
    } else {
      mse_pos[i] = mse_fix[i] = time_pos[i] = time_fix[i] = NA_real_
    }
  }

  ratio = mse_fix / mse_pos   # > 1 → posicional wins

  result = as.data.frame(c(
    list(Signal    = sig,
         Pkg       = "rLifting",
         Wavelet   = cfg$wavelet,
         Boundary  = cfg$extension,
         Method    = cfg$method,
         NoiseSd   = noise_sd,
         N         = sum(!is.na(mse_pos))),
    col_stats(mse_pos,  "MSEpos_"),
    col_stats(mse_fix,  "MSEfix_"),
    col_stats(ratio,    "Ratio_"),
    col_stats(time_pos, "Timepos_"),
    col_stats(time_fix, "Timefix_")
  ), stringsAsFactors = FALSE)

  saveRDS(result, job$fout)
  result
}

# ── Run ───────────────────────────────────────────────────────────────────────

if (length(jobs) > 0) {
  RNGkind("L'Ecuyer-CMRG"); set.seed(2025)
  n_cores = max(1L, detectCores() - 1L)
  t0 = proc.time()
  mclapply(jobs, run_job, mc.cores = n_cores,
           mc.preschedule = FALSE, mc.set.seed = TRUE)
  elapsed = (proc.time() - t0)[["elapsed"]]
  message(sprintf("Done in %.1f s (%.1f min)", elapsed, elapsed/60))
}

# ── Combine ───────────────────────────────────────────────────────────────────

files = list.files(TMPDIR, pattern = "\\.rds$", full.names = TRUE)
message(sprintf("Combining %d checkpoint files...", length(files)))
benchmark_rlifting_irregular = do.call(rbind, lapply(files, readRDS))
benchmark_rlifting_irregular = benchmark_rlifting_irregular[
  order(benchmark_rlifting_irregular$Signal,
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
message("Remove data/tmp_rlifting_irr/ when satisfied.")
