### 2. Script de Teste de Carga (`inst/scripts/benchmark_throughput.R`)

# Script de stress test: high frequency throughput
# Objetivo: medir a capacidade máxima de eventos por segundo (EPS)

library(rLifting)
library(microbenchmark)

cat("=== rLifting: HIGH FREQUENCY THROUGHPUT TEST ===\n")

# Configuração
W = 1024
LEVELS = 3
SCHEME = lifting_scheme("db2") # Um pouco mais pesado que Haar
N_EVENTS = 100000 # 100 mil eventos

cat(sprintf("Config: Window=%d, Levels=%d, Wavelet=DB2\n", W, LEVELS))
cat("Inicializando Engine C++...\n")

# Inicializa processador
proc = new_wavelet_stream(
  SCHEME, window_size = W, levels = LEVELS,
  method = "semisoft", update_freq = 1
  # Update threshold a cada passo (Pior caso)
)

# Gera dados dummy
input_stream = rnorm(N_EVENTS)

cat("Iniciando processamento de 100k eventos (loop R simples)...\n")

# Medição de Tempo
start_time = Sys.time()

# Loop explícito para simular chegada sequencial
# (sapply seria mais rapido, mas menos realista para streaming real)
dummy_out = numeric(N_EVENTS)
for(i in 1:N_EVENTS) {
  dummy_out[i] = proc(input_stream[i])
}

end_time = Sys.time()
duration = as.numeric(end_time - start_time) # segundos

eps = N_EVENTS / duration

cat(sprintf("\n--- RESULTADOS ---\n"))
cat(sprintf("Tempo Total:   %.4f segundos\n", duration))
cat(sprintf("Eventos (N):   %d\n", N_EVENTS))
cat(sprintf("Throughput:    %.2f eventos/segundo\n", eps))
cat(
  sprintf(
    "Latência Méd:  %.2f microsegundos/evento\n",
    (duration/N_EVENTS)*1e6
    )
  )

if (eps > 50000) {
  cat("\nSTATUS: EXCELENTE (Viável para Áudio/Sensores Alta Freq)\n")
} else if (eps > 10000) {
  cat("\nSTATUS: BOM (Viável para Sensores Industriais/Financeiro)\n")
} else {
  cat("\nSTATUS: ALERTA DE PERFORMANCE\n")
}
