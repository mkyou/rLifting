# V2 benchmark: extends benchmark_rlifting.rda with the new threshold/shrinkage
# combinations introduced in feat/boundary-modes:
#   - (universal, scad)        : SCAD shrinkage with default alpha/beta
#   - (universal_tuned, semisoft): default semisoft with alpha/beta tuned via SURE
#   - (sure, soft)             : SureShrink with canonical soft shrinkage
#   - (sure, scad)             : SureShrink with SCAD shrinkage
#
# Existing 900 checkpoints in data/tmp_rlifting/ from the v1 run are skipped.
# Tuned alpha/beta are computed once per (Signal, Wavelet) pair on a seeded
# clean realization and cached in data/tmp_rlifting_v2/alpha_beta_cache.rds.

suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(microbenchmark))
suppressPackageStartupMessages(library(rLifting))

N_SIM = 1000L
T_PTS = 1024L
SIGMA = 0.3
W = 255L
WARMUP = W
L_OFF = 4L
L_CAUS = 3L
LL_K = 4L

DJ_SIGNALS = c("blocks", "bumps", "doppler", "heavisine")
WAVELETS = c("haar", "db2", "cdf53", "cdf97", "dd4")
BOUNDARIES = c("symmetric", "local_linear", "one_sided", "periodic", "zero")
MODES = c("offline", "causal", "stream")

# NEW method-axes only. v1 covered (universal, hard|soft|semisoft) with default
# alpha/beta — those are intentionally NOT regenerated.
NEW_METHODS = list(
  list(label = "universal_scad",   threshold = "universal", shrinkage = "scad",     tuned = FALSE),
  list(label = "universal_tuned",  threshold = "universal", shrinkage = "semisoft", tuned = TRUE),
  list(label = "sure_soft",        threshold = "sure",      shrinkage = "soft",     tuned = FALSE),
  list(label = "sure_scad",        threshold = "sure",      shrinkage = "scad",     tuned = FALSE)
)

TMPDIR = "data/tmp_rlifting"
CACHEDIR = "data/tmp_rlifting_v2"
dir.create(CACHEDIR, showWarnings = FALSE, recursive = TRUE)

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

# Precompute tuned (alpha, beta) per (Signal, Wavelet) — single realization
# at a fixed seed for reproducibility. Used for the universal_tuned variant.
cache_path = file.path(CACHEDIR, "alpha_beta_cache.rds")
if (file.exists(cache_path)) {
  ab_cache = readRDS(cache_path)
} else {
  ab_cache = list()
}
for (sig in DJ_SIGNALS) for (wav in WAVELETS) {
  key = paste(sig, wav, sep = "__")
  if (!is.null(ab_cache[[key]])) next
  set.seed(2025L + sum(utf8ToInt(key)))
  pure = rLifting:::.generate_signal(sig, n = T_PTS)
  noisy = pure + rnorm(T_PTS, sd = SIGMA)
  sch = lifting_scheme(wav)
  tuned = tune_alpha_beta(noisy, sch, levels = L_OFF)
  ab_cache[[key]] = list(alpha = tuned$alpha, beta = tuned$beta)
  message(sprintf("Tuned %s: alpha=%.3f beta=%.3f", key,
                  tuned$alpha, tuned$beta))
}
saveRDS(ab_cache, cache_path)

CONFIGS = expand.grid(
  Signal = DJ_SIGNALS,
  Wavelet = WAVELETS,
  Boundary = BOUNDARIES,
  MethodLabel = vapply(NEW_METHODS, `[[`, character(1), "label"),
  Mode = MODES,
  stringsAsFactors = FALSE
)

method_lookup = setNames(NEW_METHODS,
                         vapply(NEW_METHODS, `[[`, character(1), "label"))

jobs = vector("list", nrow(CONFIGS))
for (i in seq_len(nrow(CONFIGS))) {
  cfg = CONFIGS[i, ]
  m = method_lookup[[cfg$MethodLabel]]
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
  m = job$m
  sig = cfg$Signal; wav = cfg$Wavelet
  bnd = cfg$Boundary; mode = cfg$Mode

  sch = lifting_scheme(wav)
  ab_key = paste(sig, wav, sep = "__")
  if (m$tuned) {
    ab = ab_cache[[ab_key]]
    alpha = ab$alpha; beta = ab$beta
  } else {
    alpha = 0.3; beta = 1.2
  }

  mse_full_vec = numeric(N_SIM)
  mse_settled_vec = numeric(N_SIM)
  time_total_vec = numeric(N_SIM)
  per_sample_vec = numeric(N_SIM)

  for (i in seq_len(N_SIM)) {
    pure = rLifting:::.generate_signal(sig, n = T_PTS)
    noisy = pure + rnorm(T_PTS, sd = SIGMA)

    tryCatch({
      if (mode == "offline") {
        mb = microbenchmark(
          out <- suppressWarnings(
            denoise_signal_offline(noisy, sch, levels = L_OFF,
                                   alpha = alpha, beta = beta,
                                   threshold_method = m$threshold,
                                   shrinkage = m$shrinkage,
                                   extension = bnd, ll_k = LL_K)
          ), times = 3L)
        t_total = median(mb$time) / 1e9

      } else if (mode == "causal") {
        mb = microbenchmark(
          out <- suppressWarnings(
            denoise_signal_causal(noisy, sch, levels = L_CAUS,
                                  window_size = W,
                                  alpha = alpha, beta = beta,
                                  threshold_method = m$threshold,
                                  shrinkage = m$shrinkage,
                                  extension = bnd, ll_k = LL_K)
          ), times = 1L)
        t_total = mb$time[1] / 1e9

      } else {
        proc = suppressWarnings(
          new_wavelet_stream(sch, window_size = W, levels = L_CAUS,
                             alpha = alpha, beta = beta,
                             threshold_method = m$threshold,
                             shrinkage = m$shrinkage,
                             extension = bnd, update_freq = 1L,
                             ll_k = LL_K)
        )
        out = numeric(T_PTS)
        t0 = as.numeric(Sys.time())
        for (j in seq_len(T_PTS)) out[j] = proc(noisy[j])
        t_total = as.numeric(Sys.time()) - t0
      }

      mse_full_vec[i] = mean((pure - out)^2)
      mse_settled_vec[i] = mean((pure[(WARMUP + 1L):T_PTS] -
                                   out[(WARMUP + 1L):T_PTS])^2)
      time_total_vec[i] = t_total
      per_sample_vec[i] = t_total / T_PTS * 1e6

    }, error = function(e) {
      message(sprintf("ERROR in %s__%s__%s__%s__%s: %s",
                      sig, wav, bnd, m$label, mode,
                      conditionMessage(e)))
      mse_full_vec[i] <<- NA_real_
      mse_settled_vec[i] <<- NA_real_
      time_total_vec[i] <<- NA_real_
      per_sample_vec[i] <<- NA_real_
    })
  }

  n_valid = sum(!is.na(mse_full_vec))
  row = c(
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
         Version = "v2",
         N = n_valid),
    col_stats(mse_full_vec, "MSE_"),
    col_stats(mse_settled_vec, "MSE_settled_"),
    col_stats(time_total_vec, "Time_total_"),
    col_stats(per_sample_vec, "Per_sample_us_")
  )
  result = as.data.frame(row, stringsAsFactors = FALSE)
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

# Reassemble: combine v1 + v2 with backfilled metadata.
message("Combining v1 + v2 results...")

load("data/benchmark_rlifting.rda")  # benchmark_rlifting (v1, 900 rows)
v1 = benchmark_rlifting
if (!"Version" %in% names(v1)) v1$Version = "v1"
if (!"ThresholdMethod" %in% names(v1)) v1$ThresholdMethod = "universal"
if (!"Shrinkage" %in% names(v1)) v1$Shrinkage = v1$Method
if (!"AlphaUsed" %in% names(v1)) v1$AlphaUsed = 0.3
if (!"BetaUsed" %in% names(v1)) v1$BetaUsed = 1.2

new_files = list.files(TMPDIR, pattern = "\\.rds$", full.names = TRUE)
new_files = new_files[!basename(new_files) %in%
  paste0(apply(expand.grid(DJ_SIGNALS, WAVELETS, BOUNDARIES,
                           c("hard", "soft", "semisoft"), MODES),
               1, paste, collapse = "__"), ".rds")]
v2 = if (length(new_files) > 0) {
  do.call(rbind, lapply(new_files, readRDS))
} else {
  data.frame()
}

shared = intersect(names(v1), names(v2))
benchmark_rlifting = rbind(v1[, shared, drop = FALSE],
                           v2[, shared, drop = FALSE])
rownames(benchmark_rlifting) = NULL
message(sprintf("Final: %d rows x %d cols (v1: %d, v2: %d)",
                nrow(benchmark_rlifting), ncol(benchmark_rlifting),
                nrow(v1), nrow(v2)))
save(benchmark_rlifting, file = "data/benchmark_rlifting.rda",
     compress = "xz")
message("Saved data/benchmark_rlifting.rda")
