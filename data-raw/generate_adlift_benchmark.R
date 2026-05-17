
# Generates data/benchmark_adlift.rda
#
# Design:
#   Signals : 4 Donoho-Johnstone (doppler, heavisine, bumps, blocks)
#   Configs : pred x neigh x int x clo x rule = 4x3x2x2x2 = 96
#             pred  : LinearPred, QuadPred, CubicPred, AdaptPred
#             neigh : 1, 2, 3
#             int   : TRUE, FALSE
#             clo   : FALSE, TRUE
#             rule  : "median", "mean"
#             (*mp variants excluded — incompatible API)
#   Reps    : 1000 per config x signal
#   Output  : 1 row per (signal x config) with summary stats over reps
#             Columns: Signal, Pkg, Wavelet, Boundary, N,
#                      MSE_{min,q1,median,mean,q3,max,se},
#                      Time_{min,q1,median,mean,q3,max,se}
#   Total rows : 96 x 4 = 384
#
# Checkpoint/resume: each job writes a 1-row .rds to data/tmp_adlift/.
# Runtime estimate: ~13 h on 11 cores (adlift ~1.3 s/call).
#
# source("data-raw/generate_adlift_benchmark.R")

suppressPackageStartupMessages({
  library(parallel)
  library(adlift)
  library(microbenchmark)
  library(rLifting)
})

# ── Summary helper ────────────────────────────────────────────────────────────

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

# ── Configuration grid ────────────────────────────────────────────────────────

PRED_LIST = list(LinearPred=LinearPred, QuadPred=QuadPred,
                 CubicPred=CubicPred,   AdaptPred=AdaptPred)

CONFIGS = do.call(rbind, lapply(names(PRED_LIST), function(pname) {
  expand.grid(pred_name=pname, neigh=1:3,
              int=c(TRUE,FALSE), clo=c(FALSE,TRUE),
              rule=c("median","mean"), stringsAsFactors=FALSE)
}))
CONFIGS$label = with(CONFIGS, paste(
  tolower(sub("Pred","",pred_name)), paste0("n",neigh),
  ifelse(int,"int","noint"), ifelse(clo,"clo","noclo"), rule, sep="_"))

DJ_SIGNALS = c("doppler","heavisine","bumps","blocks")
N_SIM      = 1000
T_PTS      = 1024
SIGMA      = 0.3
x_grid     = seq_len(T_PTS)

TMPDIR = "data/tmp_adlift"
dir.create(TMPDIR, showWarnings=FALSE, recursive=TRUE)

jobs = list()
for (sig in DJ_SIGNALS) {
  for (ci in seq_len(nrow(CONFIGS))) {
    cfg  = CONFIGS[ci,]
    fout = file.path(TMPDIR, sprintf("%s__%s.rds", sig, cfg$label))
    if (!file.exists(fout))
      jobs[[length(jobs)+1]] = list(signal=sig, cfg=cfg, fout=fout)
  }
}

message(sprintf("Configs: %d | Signals: %d | Total rows: %d",
                nrow(CONFIGS), length(DJ_SIGNALS), nrow(CONFIGS)*length(DJ_SIGNALS)))
message(sprintf("Jobs to run: %d  (skipping: %d already done)",
                length(jobs), nrow(CONFIGS)*length(DJ_SIGNALS)-length(jobs)))
message(sprintf("Estimated runtime: %.1f h on %d cores",
                length(jobs)*N_SIM*1.3/max(1L,detectCores()-1L)/3600,
                max(1L,detectCores()-1L)))

# ── Worker ────────────────────────────────────────────────────────────────────

run_job = function(job) {
  sig  = job$signal; cfg = job$cfg
  pred = PRED_LIST[[cfg$pred_name]]
  mse_vec  = numeric(N_SIM)
  time_vec = numeric(N_SIM)

  for (i in seq_len(N_SIM)) {
    pure  = rLifting:::.generate_signal(sig, n=T_PTS)
    noisy = pure + rnorm(T_PTS, sd=SIGMA)
    res = tryCatch({
      mb = microbenchmark(
        out <- as.vector(adlift::denoise(x_grid, noisy,
          pred=pred, neigh=cfg$neigh, int=cfg$int, clo=cfg$clo, rule=cfg$rule)),
        times=1L)
      list(time=mb$time[1]/1e9, mse=mean((pure-out)^2))
    }, error=function(e) NULL)
    mse_vec[i]  = if (!is.null(res)) res$mse  else NA_real_
    time_vec[i] = if (!is.null(res)) res$time else NA_real_
  }

  row = c(
    list(Signal=sig, Pkg="adlift", Wavelet=cfg$pred_name,
         Boundary=cfg$label, N=sum(!is.na(mse_vec))),
    col_stats(mse_vec,  "MSE_"),
    col_stats(time_vec, "Time_")
  )
  result = as.data.frame(row, stringsAsFactors=FALSE)
  saveRDS(result, job$fout)
  result
}

# ── Parallel execution ────────────────────────────────────────────────────────

if (length(jobs) > 0) {
  RNGkind("L'Ecuyer-CMRG"); set.seed(2025)
  n_cores = max(1L, detectCores()-1L)
  t0 = proc.time()
  mclapply(jobs, run_job, mc.cores=n_cores,
           mc.preschedule=FALSE, mc.set.seed=TRUE)
  elapsed = (proc.time()-t0)[["elapsed"]]
  message(sprintf("Done in %.1f s  (%.1f h)", elapsed, elapsed/3600))
}

# ── Combine and save ──────────────────────────────────────────────────────────

files = list.files(TMPDIR, pattern="\\.rds$", full.names=TRUE)
message(sprintf("Combining %d checkpoint files...", length(files)))
benchmark_adlift = do.call(rbind, lapply(files, readRDS))
message(sprintf("Final: %d rows x %d cols", nrow(benchmark_adlift), ncol(benchmark_adlift)))
save(benchmark_adlift, file="data/benchmark_adlift.rda", compress="xz")
message("Saved data/benchmark_adlift.rda  — remove data/tmp_adlift/ when satisfied.")
