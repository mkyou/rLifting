library(rLifting)
library(dplyr)
library(microbenchmark)

set.seed(2025)

pkgs_ext = c("wavethresh", "adlift", "nlt")
avail = pkgs_ext[pkgs_ext %in% rownames(installed.packages())]
message("External packages available: ", paste(avail, collapse = ", "))
for (p in avail) library(p, character.only = TRUE)

DJ_SIGNALS = c("doppler", "heavisine", "bumps", "blocks")
RL_WAVELETS = c("haar", "db2", "cdf53", "cdf97", "dd4")
RL_BOUNDARIES = c("symmetric", "local_linear", "one_sided")
SIGMA = 0.3
T_PTS = 1024
LEVELS_OFFLINE = floor(log2(T_PTS))
N_SIM_RL = 50
N_SIM_EXT = 20

message("Generating Doppler example...")
n_samples = 2048
signal_pure = rLifting:::.generate_signal("doppler", n_samples)
doppler_example = data.frame(
  index = seq_len(n_samples),
  original = signal_pure,
  noisy = signal_pure + rnorm(n_samples, sd = 0.5)
)

message("Running offline benchmark...")

rows = list()

for (sig_type in DJ_SIGNALS) {
  message("  Signal: ", sig_type)

  for (wname in RL_WAVELETS) {
    for (bmode in RL_BOUNDARIES) {
      sch = lifting_scheme(wname)
      for (i in seq_len(N_SIM_RL)) {
        pure = rLifting:::.generate_signal(sig_type, n = T_PTS)
        noisy = pure + rnorm(T_PTS, sd = SIGMA)

        mb = microbenchmark(
          res <- denoise_signal_offline(
            noisy, sch, levels = LEVELS_OFFLINE,
            method = "semisoft", extension = bmode
          ),
          times = 3L
        )
        rows[[length(rows) + 1]] = data.frame(
          Signal = sig_type,
          Pkg = "rLifting",
          Wavelet = wname,
          Boundary = bmode,
          Sim = i,
          Time = median(mb$time) / 1e9,
          MSE = mean((pure - res)^2)
        )
      }
    }
  }

  if ("wavethresh" %in% avail) {
    for (i in seq_len(N_SIM_EXT)) {
      pure = rLifting:::.generate_signal(sig_type, n = T_PTS)
      noisy = pure + rnorm(T_PTS, sd = SIGMA)
      tryCatch({
        mb = microbenchmark({
          wd_obj = wavethresh::wd(noisy, filter.number = 1,
            family = "DaubExPhase", type = "wavelet")
          th_obj = wavethresh::threshold(wd_obj, policy = "universal",
            type = "soft")
          out_wt = wavethresh::wr(th_obj)
        }, times = 3L)
        rows[[length(rows) + 1]] = data.frame(
          Signal = sig_type, Pkg = "wavethresh", Wavelet = "haar",
          Boundary = NA_character_, Sim = i,
          Time = median(mb$time) / 1e9,
          MSE = mean((pure - out_wt)^2)
        )
      }, error = function(e) message("    wavethresh failed: ", e$message))
    }
  }

  if ("adlift" %in% avail) {
    x_grid = seq_len(T_PTS)
    for (i in seq_len(N_SIM_EXT)) {
      pure = rLifting:::.generate_signal(sig_type, n = T_PTS)
      noisy = pure + rnorm(T_PTS, sd = SIGMA)
      tryCatch({
        mb = microbenchmark(
          out_al <- as.vector(adlift::denoise(x_grid, noisy)),
          times = 1L)
        rows[[length(rows) + 1]] = data.frame(
          Signal = sig_type, Pkg = "adlift", Wavelet = "haar",
          Boundary = NA_character_, Sim = i,
          Time = mb$time[1] / 1e9,
          MSE = mean((pure - out_al)^2)
        )
      }, error = function(e) message("    adlift failed: ", e$message))
    }
  }

  if ("nlt" %in% avail) {
    x_grid = seq_len(T_PTS)
    for (i in seq_len(N_SIM_EXT)) {
      pure = rLifting:::.generate_signal(sig_type, n = T_PTS)
      noisy = pure + rnorm(T_PTS, sd = SIGMA)
      tryCatch({
        mb = microbenchmark(
          out_nlt <- as.vector(nlt::denoiseperm(x_grid, noisy)),
          times = 1L
        )
        rows[[length(rows) + 1]] = data.frame(
          Signal = sig_type, Pkg = "nlt", Wavelet = "haar",
          Boundary = NA_character_, Sim = i,
          Time = mb$time[1] / 1e9,
          MSE = mean((pure - out_nlt)^2)
        )
      }, error = function(e) message("    nlt failed: ", e$message))
    }
  }
}

benchmark_offline = do.call(rbind, rows)

message("Running causal benchmark...")
W_SIZE = 128
LEVELS_CAUSAL = floor(log2(W_SIZE))
N_TEST = 500
pure_caus = rLifting:::.generate_signal("heavisine", n = N_TEST)
noisy_caus = pure_caus + rnorm(N_TEST, sd = SIGMA)

res_rl_causal = denoise_signal_causal(
  noisy_caus, lifting_scheme("haar"),
  window_size = W_SIZE, levels = LEVELS_CAUSAL
)
mb_rl_causal = microbenchmark(
  denoise_signal_causal(noisy_caus, lifting_scheme("haar"),
    window_size = W_SIZE, levels = LEVELS_CAUSAL),
  times = 10L
)

t_wt_naive = NA_real_
mse_wt_naive = NA_real_
if ("wavethresh" %in% avail) {
  naive_causal = function(series, w) {
    out = numeric(length(series))
    out[seq_len(w - 1)] = series[seq_len(w - 1)]
    for (k in w:length(series)) {
      block = series[(k - w + 1):k]
      suppressWarnings({
        wd_obj = wavethresh::wd(block, filter.number = 1,
          family = "DaubExPhase", type = "wavelet")
        th_obj = wavethresh::threshold(wd_obj, policy = "universal",
          type = "soft")
        out[k] = tail(wavethresh::wr(th_obj), 1)
      })
    }
    out
  }
  tryCatch({
    res_naive = naive_causal(noisy_caus, W_SIZE)
    mb_naive = microbenchmark(naive_causal(noisy_caus, W_SIZE), times = 3L)
    t_wt_naive = median(mb_naive$time) / 1e9
    mse_wt_naive = mean((pure_caus - res_naive)^2)
  }, error = function(e) message("naive causal error: ", e$message))
}

benchmark_causal = list(
  rLifting_Time_Avg = median(mb_rl_causal$time) / 1e9,
  Wavethresh_Naive_Time = t_wt_naive,
  Speedup_Factor = if (!is.na(t_wt_naive))
    t_wt_naive / (median(mb_rl_causal$time) / 1e9) else NA,
  rLifting_MSE = mean((pure_caus - res_rl_causal)^2),
  Wavethresh_Naive_MSE = mse_wt_naive
)

message("Running leakage test...")
set.seed(2025)
n_leak = 256
t_change = 128
noise_lk = rnorm(n_leak, sd = 0.5)
signal_A = noise_lk
signal_B = noise_lk + c(rep(0, t_change), rep(5, n_leak - t_change))
pre = seq_len(t_change - 1)

out_A_causal = denoise_signal_causal(signal_A, lifting_scheme("cdf97"),
  window_size = 64, levels = floor(log2(64)), method = "semisoft")
out_B_causal = denoise_signal_causal(signal_B, lifting_scheme("cdf97"),
  window_size = 64, levels = floor(log2(64)), method = "semisoft")
out_A_offline = denoise_signal_offline(signal_A, lifting_scheme("cdf97"),
  levels = floor(log2(n_leak)), method = "semisoft")
out_B_offline = denoise_signal_offline(signal_B, lifting_scheme("cdf97"),
  levels = floor(log2(n_leak)), method = "semisoft")

out_A_wt = rep(NA_real_, n_leak)
out_B_wt = rep(NA_real_, n_leak)
if ("wavethresh" %in% avail) {
  tryCatch({
    wd_A = wavethresh::wd(signal_A, filter.number = 4,
      family = "DaubExPhase", type = "wavelet")
    out_A_wt = wavethresh::wr(wavethresh::threshold(wd_A,
      policy = "universal", type = "soft"))
    wd_B = wavethresh::wd(signal_B, filter.number = 4,
      family = "DaubExPhase", type = "wavelet")
    out_B_wt = wavethresh::wr(wavethresh::threshold(wd_B,
      policy = "universal", type = "soft"))
  }, error = function(e) message("wavethresh leakage error: ", e$message))
}

leakage_results = data.frame(
  Method = c("rLifting causal (CDF 9/7)", "rLifting offline (CDF 9/7)",
             "wavethresh offline (D8)"),
  Leakage = c(
    sum((out_B_causal[pre]  - out_A_causal[pre])^2),
    sum((out_B_offline[pre] - out_A_offline[pre])^2),
    sum((out_B_wt[pre]      - out_A_wt[pre])^2, na.rm = TRUE)
  )
)

message("Saving data...")
usethis::use_data(doppler_example, benchmark_offline, benchmark_causal,
                  leakage_results, overwrite = TRUE, compress = "xz")
message("Done!")
