
# Generates data/benchmark_wavethresh.rda
#
# Design:
#   Signals  : 4 Donoho-Johnstone (doppler, heavisine, bumps, blocks)
#   Wavelets : DaubExPhase fn 1-10, DaubLeAsymm fn 4-10, Coiflets fn 1-5 = 22
#   Policies :
#     universal, cv, fdr  x  {hard, soft} = 132 configs
#     sure                x  soft only    =  22 configs  (hard unsupported by SURE)
#     BayesThresh (EbayesThresh) x soft   =  22 configs  (type internally ignored)
#   Total configs : 176
#   Reps          : 1000 per config x signal
#   Output        : 1 row per (signal x config) with summary stats over reps
#                   Columns: Signal, Pkg, Wavelet, Boundary, N,
#                            MSE_{min,q1,median,mean,q3,max,se},
#                            Time_{min,q1,median,mean,q3,max,se}
#   Total rows    : 176 x 4 = 704
#
# Checkpoint/resume: each job writes a 1-row .rds to data/tmp_wt/.
#   Re-running the script skips finished jobs automatically.
#
# source("data-raw/generate_wavethresh_benchmark.R")

suppressPackageStartupMessages({
  library(parallel)
  library(wavethresh)
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

# ── Wavelet grid ──────────────────────────────────────────────────────────────

make_label = function(fam, fn)
  paste0(switch(fam, DaubExPhase="db", DaubLeAsymm="la", Coiflets="co"), fn)

WAVELET_GRID = rbind(
  data.frame(family="DaubExPhase", fn=1:10,  stringsAsFactors=FALSE),
  data.frame(family="DaubLeAsymm", fn=4:10,  stringsAsFactors=FALSE),
  data.frame(family="Coiflets",    fn=1:5,   stringsAsFactors=FALSE)
)
WAVELET_GRID$wavelet = mapply(make_label, WAVELET_GRID$family, WAVELET_GRID$fn)

POLICY_TYPE = rbind(
  expand.grid(policy=c("universal","cv","fdr"), type=c("soft","hard"),
              stringsAsFactors=FALSE),
  data.frame(policy="sure",        type="soft", stringsAsFactors=FALSE),
  data.frame(policy="BayesThresh", type="soft", stringsAsFactors=FALSE)
)

CONFIGS    = merge(WAVELET_GRID, POLICY_TYPE)
DJ_SIGNALS = c("doppler", "heavisine", "bumps", "blocks")
N_SIM      = 1000
T_PTS      = 1024
SIGMA      = 0.3

TMPDIR = "data/tmp_wt"
dir.create(TMPDIR, showWarnings=FALSE, recursive=TRUE)

# ── Job list (skip already-done checkpoints) ──────────────────────────────────

jobs = list()
for (sig in DJ_SIGNALS) {
  for (ci in seq_len(nrow(CONFIGS))) {
    cfg  = CONFIGS[ci, ]
    key  = sprintf("%s__%s__%s__%s", sig, cfg$wavelet, cfg$policy, cfg$type)
    fout = file.path(TMPDIR, paste0(key, ".rds"))
    if (!file.exists(fout))
      jobs[[length(jobs)+1]] = list(signal=sig, cfg=cfg, fout=fout)
  }
}

message(sprintf("Configs: %d | Signals: %d | Total rows: %d",
                nrow(CONFIGS), length(DJ_SIGNALS), nrow(CONFIGS)*length(DJ_SIGNALS)))
message(sprintf("Jobs to run: %d  (skipping: %d already done)",
                length(jobs),
                nrow(CONFIGS)*length(DJ_SIGNALS) - length(jobs)))

# ── Worker ────────────────────────────────────────────────────────────────────

run_job = function(job) {
  sig = job$signal; cfg = job$cfg
  mse_vec  = numeric(N_SIM)
  time_vec = numeric(N_SIM)

  for (i in seq_len(N_SIM)) {
    pure  = rLifting:::.generate_signal(sig, n=T_PTS)
    noisy = pure + rnorm(T_PTS, sd=SIGMA)
    res = tryCatch({
      mb = microbenchmark({
        wd_obj = wavethresh::wd(noisy, filter.number=cfg$fn,
                                family=cfg$family, type="wavelet")
        th_obj = wavethresh::threshold(wd_obj, policy=cfg$policy, type=cfg$type)
        out    = wavethresh::wr(th_obj)
      }, times=3L)
      list(time=median(mb$time)/1e9, mse=mean((pure-out)^2))
    }, error=function(e) NULL)
    mse_vec[i]  = if (!is.null(res)) res$mse  else NA_real_
    time_vec[i] = if (!is.null(res)) res$time else NA_real_
  }

  row = c(
    list(Signal=sig, Pkg="wavethresh", Wavelet=cfg$wavelet,
         Boundary=paste0(cfg$policy,"_",cfg$type), N=sum(!is.na(mse_vec))),
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
  message(sprintf("Done in %.1f s  (%.1f min)", elapsed, elapsed/60))
}

# ── Combine and save ──────────────────────────────────────────────────────────

files = list.files(TMPDIR, pattern="\\.rds$", full.names=TRUE)
message(sprintf("Combining %d checkpoint files...", length(files)))
benchmark_wavethresh = do.call(rbind, lapply(files, readRDS))
message(sprintf("Final: %d rows x %d cols", nrow(benchmark_wavethresh),
                ncol(benchmark_wavethresh)))
save(benchmark_wavethresh, file="data/benchmark_wavethresh.rda", compress="xz")
message("Saved data/benchmark_wavethresh.rda  — remove data/tmp_wt/ when satisfied.")
