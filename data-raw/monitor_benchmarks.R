suppressPackageStartupMessages(library(dplyr))

check_pkg = function(pkg, tmpdir, total_jobs) {
  files = list.files(tmpdir, pattern = "\\.rds$", full.names = TRUE)
  n_done = length(files)
  if (n_done == 0) {
    cat(sprintf("  %-12s %4d / %4d  (0%%) -- not started\n",
      pkg, 0L, total_jobs))
    return(invisible(NULL))
  }

  pct = round(100 * n_done / total_jobs)
  df = do.call(rbind, lapply(files, readRDS))

  cat(sprintf("  %-12s %4d / %4d  (%d%%)\n", pkg, n_done, total_jobs, pct))
  cat(sprintf(
    "    MSE  median  -- min: %.5f  mean: %.5f  max: %.5f\n",
    min(df$MSE_median, na.rm = TRUE),
    mean(df$MSE_median, na.rm = TRUE),
    max(df$MSE_median, na.rm = TRUE)))
  cat(sprintf(
    "    Time median  -- min: %.3f ms  mean: %.3f ms  max: %.3f ms\n",
    min(df$Time_median, na.rm = TRUE) * 1e3,
    mean(df$Time_median, na.rm = TRUE) * 1e3,
    max(df$Time_median, na.rm = TRUE) * 1e3))

  top = df[order(df$MSE_median), ][seq_len(min(3L, nrow(df))), ]
  cat("    Best MSE so far:\n")
  for (i in seq_len(nrow(top))) {
    cat(sprintf("      [%d] %s | %s | %s  MSE=%.5f\n",
      i, top$Signal[i], top$Wavelet[i], top$Boundary[i],
      top$MSE_median[i]))
  }
  invisible(df)
}

cat("=== Benchmark progress ===\n")
cat(sprintf("    %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

check_pkg("wavethresh", "data/tmp_wt",       total_jobs = 704L)
check_pkg("adlift",     "data/tmp_adlift",   total_jobs = 384L)
check_pkg("nlt",        "data/tmp_nlt",      total_jobs = 384L)
check_pkg("rLifting",   "data/tmp_rlifting", total_jobs = 360L)
