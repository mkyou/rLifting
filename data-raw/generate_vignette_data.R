
# Script to generate comparison data for vignettes
# Run this from the package root: source("data-raw/generate_vignette_data.R")

library(rLifting)
library(dplyr)
library(microbenchmark)

set.seed(2025)

# --- Check for Suggested Packages ---
pkgs = c("wavethresh", "adlift", "nlt")
installed_pkgs = pkgs[pkgs %in% rownames(installed.packages())]

message("Packages available for benchmark: ", paste(installed_pkgs, collapse = ", "))

# Load all installed packages globally to avoid scope issues
for(p in installed_pkgs) {
  library(p, character.only = TRUE)
}

# --- 1. Doppler Example (Base) ---
message("Generating Doppler example...")
n_samples = 2048
signal_pure = rLifting:::.generate_signal("doppler", n_samples)
noise = rnorm(n_samples, sd = 0.5)
signal_noisy = signal_pure + noise

doppler_example = data.frame(
  index = 1:n_samples,
  original = signal_pure,
  noisy = signal_noisy
)

# --- 2. Offline Benchmark (Speed & Accuracy) ---
message("Running offline benchmark (Haar, N=1024)...")

# Config
N_SIM = 50
T_PTS = 1024
SIG_TYPE = "doppler"
LEVELS_OFFLINE = floor(log2(T_PTS))

results_offline = list()
counter = 1

for(i in 1:N_SIM) {
  if(i %% 10 == 0) message("  Simulation ", i, "/", N_SIM)

  # Generate Data
  pure = rLifting:::.generate_signal(SIG_TYPE, n = T_PTS)
  noisy = pure + rnorm(T_PTS, sd = 0.3)
  x_grid = 1:T_PTS

  # -- rLifting --
  # Compute result for MSE
  res_rl = denoise_signal_offline(
    noisy, lifting_scheme("haar"), levels = LEVELS_OFFLINE, method = "semisoft"
  )
  # Time with microbenchmark (nanosecond precision)
  mb_rl = microbenchmark(
    denoise_signal_offline(noisy, lifting_scheme("haar"), levels = LEVELS_OFFLINE, method = "semisoft"),
    times = 5L
  )
  t_rl = median(mb_rl$time) / 1e9
  mse_rl = mean((pure - res_rl)^2)
  results_offline[[counter]] = data.frame(Pkg="rLifting", Time=t_rl, MSE=mse_rl)
  counter = counter + 1

  # -- wavethresh --
  if("wavethresh" %in% installed_pkgs) {
    tryCatch({
      out_wt = {
        wd_obj = wavethresh::wd(noisy, filter.number=1, family="DaubExPhase", type="wavelet")
        th_obj = wavethresh::threshold(wd_obj, policy="universal", type="soft")
        wavethresh::wr(th_obj)
      }
      mb_wt = microbenchmark({
        wd_obj = wavethresh::wd(noisy, filter.number=1, family="DaubExPhase", type="wavelet")
        th_obj = wavethresh::threshold(wd_obj, policy="universal", type="soft")
        wavethresh::wr(th_obj)
      }, times = 5L)
      t_wt = median(mb_wt$time) / 1e9
      mse_wt = mean((pure - out_wt)^2)
      results_offline[[counter]] = data.frame(Pkg="wavethresh", Time=t_wt, MSE=mse_wt)
    }, error = function(e) {
      message("wavethresh failed: ", e$message)
      results_offline[[counter]] <<- data.frame(Pkg="wavethresh", Time=NA, MSE=NA)
    })
  } else {
    results_offline[[counter]] = data.frame(Pkg="wavethresh", Time=NA, MSE=NA)
  }
  counter = counter + 1

  # -- adlift --
  if("adlift" %in% installed_pkgs) {
    tryCatch({
      out_al = as.vector(denoise(x_grid, noisy))
      mb_al = microbenchmark(denoise(x_grid, noisy), times = 1L)
      t_al = mb_al$time[1] / 1e9
      mse_al = mean((pure - out_al)^2)
      results_offline[[counter]] = data.frame(Pkg="adlift", Time=t_al, MSE=mse_al)
    }, error = function(e) {
      message("adlift failed: ", e$message)
      results_offline[[counter]] <<- data.frame(Pkg="adlift", Time=NA, MSE=NA)
    })
  } else {
    results_offline[[counter]] = data.frame(Pkg="adlift", Time=NA, MSE=NA)
  }
  counter = counter + 1

  # -- nlt --
  if("nlt" %in% installed_pkgs) {
    tryCatch({
      out_nlt = as.vector(denoiseperm(x_grid, noisy))
      mb_nlt = microbenchmark(denoiseperm(x_grid, noisy), times = 1L)
      t_nlt = mb_nlt$time[1] / 1e9
      mse_nlt = mean((pure - out_nlt)^2)
      results_offline[[counter]] = data.frame(Pkg="nlt", Time=t_nlt, MSE=mse_nlt)
    }, error = function(e) {
      message("nlt failed: ", e$message)
      results_offline[[counter]] <<- data.frame(Pkg="nlt", Time=NA, MSE=NA)
    })
  } else {
    results_offline[[counter]] = data.frame(Pkg="nlt", Time=NA, MSE=NA)
  }
  counter = counter + 1
}

benchmark_offline = do.call(rbind, results_offline)

# --- 3. Causal Benchmark (Speed + MSE) ---
message("Running causal benchmark...")
W_SIZE = 128
LEVELS_CAUSAL = floor(log2(W_SIZE))
N_TEST = 500
pure = rLifting:::.generate_signal("heavisine", n = N_TEST)
noisy = pure + rnorm(N_TEST, sd = 0.3)

# -- rLifting causal --
res_rl_causal = denoise_signal_causal(
  noisy, lifting_scheme("haar"), window_size = W_SIZE, levels = LEVELS_CAUSAL
)
mb_rl_causal = microbenchmark(
  denoise_signal_causal(noisy, lifting_scheme("haar"),
                        window_size = W_SIZE, levels = LEVELS_CAUSAL),
  times = 10L
)
mse_rl_causal = mean((pure - res_rl_causal)^2)

# -- Wavethresh Naive Loop --
t_wt_naive = NA
mse_wt_naive = NA
if("wavethresh" %in% installed_pkgs) {
  naive_causal_func = function(series, w) {
    out = numeric(length(series))
    out[1:(w-1)] = series[1:(w-1)]
    for(k in w:length(series)) {
      block = series[(k-w+1):k]
      suppressWarnings({
        wd_obj = wavethresh::wd(block, filter.number=1,
                                  family="DaubExPhase", type="wavelet")
        th_obj = wavethresh::threshold(wd_obj, policy="universal", type="soft")
        rec = wavethresh::wr(th_obj)
      })
      out[k] = tail(rec, 1)
    }
    return(out)
  }
  tryCatch({
    res_naive = naive_causal_func(noisy, W_SIZE)
    mb_naive = microbenchmark(naive_causal_func(noisy, W_SIZE), times = 3L)
    t_wt_naive = median(mb_naive$time) / 1e9
    mse_wt_naive = mean((pure - res_naive)^2)
  }, error = function(e) message("naive causal error: ", e$message))
}

benchmark_causal = list(
  rLifting_Time_Avg = median(mb_rl_causal$time) / 1e9,
  Wavethresh_Naive_Time = t_wt_naive,
  Speedup_Factor = if(!is.na(t_wt_naive)) t_wt_naive / (median(mb_rl_causal$time) / 1e9) else NA,
  rLifting_MSE = mse_rl_causal,
  Wavethresh_Naive_MSE = mse_wt_naive
)

# --- 4. Leakage test (counterfactual) ---
# Run each filter on two signals that share identical noise for t < t_change.
#   Signal A: noise only (no event)
#   Signal B: same noise, but with a step of amplitude 5 added at t >= t_change
# Leakage = sum((output_B[pre] - output_A[pre])^2)
# A causal filter → leakage = 0
# An offline filter with long support → leakage > 0

message("Running leakage test...")
set.seed(2025)
n_leak = 256
t_change = 128
LEVELS_LEAK = floor(log2(n_leak))

noise = rnorm(n_leak, sd = 0.5)
signal_A = noise
signal_B = noise + c(rep(0, t_change), rep(5, n_leak - t_change))

pre = 1:(t_change - 1)

# --- rLifting causal (CDF 9/7) ---
out_A_causal = denoise_signal_causal(signal_A, lifting_scheme("cdf97"),
  window_size = 64, levels = floor(log2(64)), method = "semisoft")
out_B_causal = denoise_signal_causal(signal_B, lifting_scheme("cdf97"),
  window_size = 64, levels = floor(log2(64)), method = "semisoft")

# --- rLifting offline (CDF 9/7) ---
out_A_offline = denoise_signal_offline(signal_A, lifting_scheme("cdf97"),
  levels = LEVELS_LEAK, method = "semisoft")
out_B_offline = denoise_signal_offline(signal_B, lifting_scheme("cdf97"),
  levels = LEVELS_LEAK, method = "semisoft")

# --- wavethresh offline (D8, comparable support to CDF 9/7) ---
out_A_wt = rep(NA_real_, n_leak)
out_B_wt = rep(NA_real_, n_leak)
if("wavethresh" %in% installed_pkgs) {
  tryCatch({
    wd_A = wavethresh::wd(signal_A, filter.number=4, family="DaubExPhase", type="wavelet")
    th_A = wavethresh::threshold(wd_A, policy="universal", type="soft")
    out_A_wt = wavethresh::wr(th_A)

    wd_B = wavethresh::wd(signal_B, filter.number=4, family="DaubExPhase", type="wavelet")
    th_B = wavethresh::threshold(wd_B, policy="universal", type="soft")
    out_B_wt = wavethresh::wr(th_B)
  }, error = function(e) message("wavethresh leakage error: ", e$message))
}

leakage_results = data.frame(
  Method = c("rLifting causal (CDF 9/7)", "rLifting offline (CDF 9/7)", "wavethresh offline (D8)"),
  Leakage = c(
    sum((out_B_causal[pre] - out_A_causal[pre])^2),
    sum((out_B_offline[pre] - out_A_offline[pre])^2),
    sum((out_B_wt[pre] - out_A_wt[pre])^2, na.rm = TRUE)
  )
)

message("Saving data...")
usethis::use_data(doppler_example, benchmark_offline, benchmark_causal,
                  leakage_results, overwrite = TRUE)
message("Done!")




