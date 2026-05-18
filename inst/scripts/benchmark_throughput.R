library(rLifting)
library(microbenchmark)

cat("=== rLifting: HIGH FREQUENCY THROUGHPUT TEST ===\n")

W = 1024
LEVELS = 3
SCHEME = lifting_scheme("db2")
N_EVENTS = 100000

cat(sprintf("Config: Window=%d, Levels=%d, Wavelet=DB2\n", W, LEVELS))

proc = new_wavelet_stream(
  SCHEME, window_size = W, levels = LEVELS,
  method = "semisoft", update_freq = 1
)

input_stream = rnorm(N_EVENTS)

start_time = Sys.time()

dummy_out = numeric(N_EVENTS)
for (i in 1:N_EVENTS) {
  dummy_out[i] = proc(input_stream[i])
}

end_time = Sys.time()
duration = as.numeric(end_time - start_time)

eps = N_EVENTS / duration

cat(sprintf("\n--- RESULTS ---\n"))
cat(sprintf("Total Time:  %.4f s\n", duration))
cat(sprintf("Events (N):  %d\n", N_EVENTS))
cat(sprintf("Throughput:  %.2f events/s\n", eps))
cat(sprintf("Avg Latency: %.2f us/event\n", (duration/N_EVENTS)*1e6))

if (eps > 50000) {
  cat("\nSTATUS: EXCELLENT\n")
} else if (eps > 10000) {
  cat("\nSTATUS: GOOD\n")
} else {
  cat("\nSTATUS: PERFORMANCE WARNING\n")
}
