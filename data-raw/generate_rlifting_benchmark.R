
# Generates data/benchmark_rlifting.rda
#
# Design:
#   Signals    : blocks, bumps, doppler, heavisine                  (4)
#   Wavelets   : haar, db2, cdf53, cdf97, dd4                       (5)
#   Boundaries : symmetric, local_linear, one_sided, periodic, zero (5)
#   Methods    : hard, soft, semisoft                               (3)
#   Modes      : offline, causal, stream                            (3)
#   Total rows : 4 x 5 x 5 x 3 x 3 = 900
#   Reps       : 1000 | n : 1024 | sigma : 0.3
#
# Parameters:
#   Offline : levels = 4, ll_k = 4
#   Causal  : window_size = 255, levels = 3, ll_k = 4
#   Stream  : window_size = 255, levels = 3, ll_k = 4
#
# MSE columns (all modes):
#   MSE_*          : on full T_PTS samples
#   MSE_settled_*  : on (WARMUP+1):T_PTS — directly comparable across modes
#
# Time columns:
#   Time_total_*     : total wall time for full T_PTS-sample pipeline (seconds)
#   Per_sample_us_*  : per-sample time (µs)
#     offline/causal : Time_total / T_PTS * 1e6
#     stream         : settled loop time / (T_PTS - WARMUP) * 1e6
#
# Checkpoint/resume: data/tmp_rlifting/
#
# Runtime (measured on 11 cores, cdf53/symmetric/semisoft as reference):
#   Offline : ~0.07 ms/rep  →  ~0.2 s/job  →  ~5   s total
#   Causal  : ~6.8  ms/rep  →  ~6.8 s/job  →  ~3   min total
#   Stream  : ~29   ms/rep  →  ~29  s/job  →  ~13  min total
#   Estimated total: ~16–20 min on 10 cores
#
# source("data-raw/generate_rlifting_benchmark.R")

suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(microbenchmark))
suppressPackageStartupMessages(library(rLifting))

# ── Constants ─────────────────────────────────────────────────────────────────

N_SIM   = 1000L
T_PTS   = 1024L
SIGMA   = 0.3
W       = 255L    # window_size — odd, close to 256
WARMUP  = W       # samples to exclude from settled MSE / per-sample timing
L_OFF   = 4L      # levels for offline
L_CAUS  = 3L      # levels for causal / stream
LL_K    = 4L      # local_linear OLS points

DJ_SIGNALS = c("blocks", "bumps", "doppler", "heavisine")
WAVELETS   = c("haar", "db2", "cdf53", "cdf97", "dd4")
BOUNDARIES = c("symmetric", "local_linear", "one_sided", "periodic", "zero")
METHODS    = c("hard", "soft", "semisoft")
MODES      = c("offline", "causal", "stream")

TMPDIR = "data/tmp_rlifting"
dir.create(TMPDIR, showWarnings = FALSE, recursive = TRUE)

# ── Summary helper ────────────────────────────────────────────────────────────

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

# ── Job list ──────────────────────────────────────────────────────────────────

CONFIGS = expand.grid(
  Signal   = DJ_SIGNALS,
  Wavelet  = WAVELETS,
  Boundary = BOUNDARIES,
  Method   = METHODS,
  Mode     = MODES,
  stringsAsFactors = FALSE
)

jobs = vector("list", nrow(CONFIGS))
for (i in seq_len(nrow(CONFIGS))) {
  cfg  = CONFIGS[i, ]
  key  = paste(cfg$Signal, cfg$Wavelet, cfg$Boundary,
               cfg$Method, cfg$Mode, sep = "__")
  fout = file.path(TMPDIR, paste0(key, ".rds"))
  jobs[[i]] = list(cfg = cfg, fout = fout, done = file.exists(fout))
}

pending = jobs[!vapply(jobs, `[[`, logical(1), "done")]
message(sprintf(
  "Configs: %d | Pending: %d | Skipping: %d already done",
  length(jobs), length(pending), length(jobs) - length(pending)
))

# ── Worker ────────────────────────────────────────────────────────────────────

run_job = function(job) {
  cfg  = job$cfg
  sig  = cfg$Signal;   wav  = cfg$Wavelet
  bnd  = cfg$Boundary; mth  = cfg$Method; mode = cfg$Mode

  sch = lifting_scheme(wav)

  mse_full_vec    = numeric(N_SIM)
  mse_settled_vec = numeric(N_SIM)
  time_total_vec  = numeric(N_SIM)
  per_sample_vec  = numeric(N_SIM)

  for (i in seq_len(N_SIM)) {
    pure  = rLifting:::.generate_signal(sig, n = T_PTS)
    noisy = pure + rnorm(T_PTS, sd = SIGMA)

    tryCatch({

      if (mode == "offline") {

        mb = microbenchmark(
          out <- suppressWarnings(
            denoise_signal_offline(noisy, sch, levels = L_OFF,
                                   method = mth, extension = bnd, ll_k = LL_K)
          ), times = 3L)
        t_total = median(mb$time) / 1e9

        mse_full_vec[i]    = mean((pure - out)^2)
        mse_settled_vec[i] = mean((pure[(WARMUP + 1L):T_PTS] -
                                    out[(WARMUP + 1L):T_PTS])^2)
        time_total_vec[i]  = t_total
        per_sample_vec[i]  = t_total / T_PTS * 1e6

      } else if (mode == "causal") {

        mb = microbenchmark(
          out <- suppressWarnings(
            denoise_signal_causal(noisy, sch, levels = L_CAUS,
                                  window_size = W, method = mth,
                                  extension = bnd, ll_k = LL_K)
          ), times = 1L)
        t_total = mb$time[1] / 1e9

        mse_full_vec[i]    = mean((pure - out)^2)
        mse_settled_vec[i] = mean((pure[(WARMUP + 1L):T_PTS] -
                                    out[(WARMUP + 1L):T_PTS])^2)
        time_total_vec[i]  = t_total
        per_sample_vec[i]  = t_total / T_PTS * 1e6

      } else {  # stream

        proc = suppressWarnings(
          new_wavelet_stream(sch, window_size = W, levels = L_CAUS,
                             method = mth, extension = bnd,
                             update_freq = 1L, ll_k = LL_K)
        )
        out = numeric(T_PTS)

        t0 = as.numeric(Sys.time())
        for (j in seq_len(T_PTS)) out[j] = proc(noisy[j])
        t_total = as.numeric(Sys.time()) - t0

        mse_full_vec[i]    = mean((pure - out)^2)
        mse_settled_vec[i] = mean((pure[(WARMUP + 1L):T_PTS] -
                                    out[(WARMUP + 1L):T_PTS])^2)
        time_total_vec[i]  = t_total
        per_sample_vec[i]  = t_total / T_PTS * 1e6
      }

    }, error = function(e) {
      mse_full_vec[i]    <<- NA_real_
      mse_settled_vec[i] <<- NA_real_
      time_total_vec[i]  <<- NA_real_
      per_sample_vec[i]  <<- NA_real_
    })
  }

  n_valid = sum(!is.na(mse_full_vec))
  row = c(
    list(Signal   = sig,
         Pkg      = "rLifting",
         Mode     = mode,
         Wavelet  = wav,
         Boundary = bnd,
         Method   = mth,
         N        = n_valid),
    col_stats(mse_full_vec,    "MSE_"),
    col_stats(mse_settled_vec, "MSE_settled_"),
    col_stats(time_total_vec,  "Time_total_"),
    col_stats(per_sample_vec,  "Per_sample_us_")
  )
  result = as.data.frame(row, stringsAsFactors = FALSE)
  saveRDS(result, job$fout)
  result
}

# ── Parallel execution ────────────────────────────────────────────────────────

if (length(pending) > 0) {
  RNGkind("L'Ecuyer-CMRG"); set.seed(2025)
  n_cores = max(1L, detectCores() - 1L)
  message(sprintf("Running on %d cores...", n_cores))
  t0 = proc.time()
  mclapply(pending, run_job, mc.cores = n_cores,
           mc.preschedule = FALSE, mc.set.seed = TRUE)
  elapsed = (proc.time() - t0)[["elapsed"]]
  message(sprintf("Done in %.1f s  (%.1f min)", elapsed, elapsed / 60))
}

# ── Combine and save ──────────────────────────────────────────────────────────

files = list.files(TMPDIR, pattern = "\\.rds$", full.names = TRUE)
message(sprintf("Combining %d checkpoint files...", length(files)))
benchmark_rlifting = do.call(rbind, lapply(files, readRDS))
rownames(benchmark_rlifting) = NULL
message(sprintf("Final: %d rows x %d cols",
                nrow(benchmark_rlifting), ncol(benchmark_rlifting)))
save(benchmark_rlifting, file = "data/benchmark_rlifting.rda", compress = "xz")
message("Saved data/benchmark_rlifting.rda — remove data/tmp_rlifting/ when satisfied.")
