
# Script to generate comparison data for vignettes
# Run this from the package root: source("data-raw/generate_vignette_data.R")

library(rLifting)
library(dplyr)
library(microbenchmark)

set.seed(2025)

# --- Check for Suggested Packages ---
pkgs <- c("wavethresh", "wavelets", "adlift", "waveslim", "nlt", "CNLTreg")
installed_pkgs <- pkgs[pkgs %in% rownames(installed.packages())]

message("Packages available for benchmark: ", paste(installed_pkgs, collapse = ", "))

# --- 1. Doppler Example (Base) ---
message("Generating Doppler Example...")
n_samples <- 2048
signal_pure <- rLifting:::.generate_signal("doppler", n = n_samples)
noise <- rnorm(n_samples, sd = 0.5)
signal_noisy <- signal_pure + noise

doppler_example <- data.frame(
  index = 1:n_samples,
  original = signal_pure,
  noisy = signal_noisy
)

# --- 2. Offline Benchmark (Speed & Accuracy) ---
message("Running Offline Benchmark (Haar, N=1024)...")

# Config
N_SIM <- 20 
T_PTS <- 1024
SIG_TYPE <- "doppler"
LEVELS_OFFLINE <- floor(log2(T_PTS))

results_offline <- list()
counter <- 1

for(i in 1:N_SIM) {
  # Generate Data
  pure <- rLifting:::.generate_signal(SIG_TYPE, n = T_PTS)
  noisy <- pure + rnorm(T_PTS, sd = 0.3)
  x_grid <- 1:T_PTS
  
  # -- rLifting --
  t_rl <- system.time({
    res_rl <- denoise_signal_offline(
      noisy, lifting_scheme("haar"), levels = LEVELS_OFFLINE, method = "soft"
    )
  })['elapsed']
  
  mse_rl <- mean((pure - res_rl)^2)
  results_offline[[counter]] <- data.frame(Pkg="rLifting", Time=as.numeric(t_rl), MSE=mse_rl)
  counter <- counter + 1
  
  # -- wavethresh --
  if("wavethresh" %in% installed_pkgs) {
    tryCatch({
      t_wt <- system.time({
        wd_obj <- wavethresh::wd(noisy, filter.number=1, family="DaubExPhase", type="wavelet")
        # Wavethresh handles levels differently (nlevels). usually nlevels = log2(N).
        # We can just threshold the whole object.
        th_obj <- wavethresh::threshold(wd_obj, policy="universal", type="soft")
        res_wt <- wavethresh::wr(th_obj)
      })['elapsed']
      mse_wt <- mean((pure - res_wt)^2)
      results_offline[[counter]] <- data.frame(Pkg="wavethresh", Time=as.numeric(t_wt), MSE=mse_wt)
      counter <- counter + 1
    }, error = function(e) message("wavethresh error: ", e$message))
  }
  
  # -- wavelets --
  if("wavelets" %in% installed_pkgs) {
    tryCatch({
      t_wv <- system.time({
        wt_obj <- wavelets::dwt(noisy, filter="haar", n.levels=LEVELS_OFFLINE)
        res_inv <- wavelets::idwt(wt_obj) 
      })['elapsed']
      results_offline[[counter]] <- data.frame(Pkg="wavelets", Time=as.numeric(t_wv), MSE=NA)
      counter <- counter + 1
    }, error = function(e) message("wavelets error: ", e$message))
  }

  # -- waveslim --
  if("waveslim" %in% installed_pkgs) {
    tryCatch({
      t_ws <- system.time({
        ws_obj <- waveslim::dwt(noisy, wf="haar", n.levels=LEVELS_OFFLINE)
        # Simple soft thresholding manually
        univ_thresh <- sqrt(2 * log(T_PTS))
        ws_thresh <- lapply(ws_obj, function(x) sign(x) * pmax(abs(x) - univ_thresh, 0))
        res_ws <- waveslim::idwt(ws_thresh)
      })['elapsed']
      mse_ws <- mean((pure - res_ws[1:T_PTS])^2)
      results_offline[[counter]] <- data.frame(Pkg="waveslim", Time=as.numeric(t_ws), MSE=mse_ws)
      counter <- counter + 1
    }, error = function(e) message("waveslim error: ", e$message))
  }

  # -- adlift --
  if("adlift" %in% installed_pkgs) {
    tryCatch({
      t_al <- system.time({
        # adlift is adaptive
        al_obj <- adlift::amt(noisy, x_grid, filternumber=1, family="DaubExPhase")
        al_den <- adlift::denoise(al_obj, f=1) 
        res_al <- al_den
      })['elapsed']
      mse_al <- mean((pure - res_al)^2)
      results_offline[[counter]] <- data.frame(Pkg="adlift", Time=as.numeric(t_al), MSE=mse_al)
      counter <- counter + 1
    }, error = function(e) message("adlift error: ", e$message))
  }

  # -- nlt --
  if("nlt" %in% installed_pkgs) {
    tryCatch({
      t_nlt <- system.time({
         nlt_res <- nlt::denoise.nlt(noisy, family="DaubExPhase", filter=1)
         res_nlt <- nlt_res$denoised
      })['elapsed']
      mse_nlt <- mean((pure - res_nlt)^2)
      results_offline[[counter]] <- data.frame(Pkg="nlt", Time=as.numeric(t_nlt), MSE=mse_nlt)
      counter <- counter + 1
    }, error = function(e) message("nlt error: ", e$message))
  }
}

benchmark_offline <- do.call(rbind, results_offline)

# --- 3. Causal Benchmark ---
message("Running Causal Benchmark...")
W_SIZE <- 128
N_TEST <- 500
LEVELS_CAUSAL <- floor(log2(W_SIZE))
pure <- rLifting:::.generate_signal("heavisine", n = N_TEST)
noisy <- pure + rnorm(N_TEST, sd = 0.3)

# -- rLifting --
t_rl_causal <- microbenchmark(
  Process = {
    res <- denoise_signal_causal(noisy, lifting_scheme("haar"), window_size=W_SIZE, levels=LEVELS_CAUSAL)
  },
  times = 5
)

# -- Wavethresh Naive --
t_wt_naive <- NA
if("wavethresh" %in% installed_pkgs) {
  naive_causal_func <- function(series, w) {
    out <- numeric(length(series))
    out[1:(w-1)] <- series[1:(w-1)]
    for(k in w:length(series)) {
      block <- series[(k-w+1):k]
      wd_obj <- wavethresh::wd(block, filter.number=1, family="DaubExPhase", type="wavelet")
      th_obj <- wavethresh::threshold(wd_obj, policy="universal", type="soft")
      rec <- wavethresh::wr(th_obj)
      out[k] <- tail(rec, 1)
    }
    return(out)
  }
  tryCatch({
    start_t <- Sys.time()
    res_naive <- naive_causal_func(noisy, W_SIZE)
    end_t <- Sys.time()
    t_wt_naive <- as.numeric(difftime(end_t, start_t, units="secs"))
  }, error = function(e) message("naive causal error: ", e$message))
}

benchmark_causal <- list(
  rLifting_Time_Avg = mean(t_rl_causal$time) / 1e9,
  Wavethresh_Naive_Time = t_wt_naive,
  Speedup_Factor = if(!is.na(t_wt_naive)) t_wt_naive / (mean(t_rl_causal$time) / 1e9) else NA
)

# --- 4. Leakage ---
message("Calculating Leakage...")
n_imp <- 256
impulse <- rep(0, n_imp)
trigger_idx <- 128
impulse[trigger_idx] <- 10

out_causal <- denoise_signal_causal(impulse, lifting_scheme("cdf97"), window_size=64, levels=3, beta=0)
out_offline <- denoise_signal_offline(impulse, lifting_scheme("cdf97"), levels=3, beta=0)

leakage_results <- data.frame(
  Mode = c("Causal", "Offline"),
  Pre_Trigger_Energy = c(sum(out_causal[1:(trigger_idx-1)]^2), sum(out_offline[1:(trigger_idx-1)]^2))
)

message("Saving data...")
usethis::use_data(doppler_example, benchmark_offline, benchmark_causal, leakage_results, overwrite = TRUE)
message("Done!")
