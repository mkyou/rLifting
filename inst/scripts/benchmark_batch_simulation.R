# Objetivo: medir a capacidade de processar LOTES de curvas
# (Cenário de Backtest)
# Comparando: offline vs causal C++ vs Causal R (naive)

library(rLifting)

# CONFIGURAÇÃO
N_CURVES = 1000 # Quantidade de sinais (ex: 1000 dias de trading)
T_LEN = 1024 # Tamanho do sinal
SCHEME = lifting_scheme("db2")
LEVELS = 3

cat("=== BENCHMARK: SIMULAÇÃO MASSIVA (BATCH) ===\n")
cat(
  sprintf(
    "Sinais: %d | Tamanho: %d pontos | Total Pontos: %d\n",
    N_CURVES, T_LEN, N_CURVES * T_LEN
    )
  )
cat("Gerando matriz de dados...\n")

# Matriz (Linhas = Tempo, Colunas = Sinais Diferentes)
set.seed(123)
data_matrix = matrix(rnorm(N_CURVES * T_LEN), nrow = T_LEN, ncol = N_CURVES)

# MODO OFFLINE (referência de velocidade)
cat("\n1. Executando Modo OFFLINE (Referência)...\n")
time_off = system.time({
  # apply percorre as colunas
  res_off = apply(data_matrix, 2, function(x) {
    denoise_signal_offline(x, SCHEME, levels = LEVELS, method="semisoft")
  })
})
cat(sprintf("   Tempo: %.4f s\n", time_off["elapsed"]))


# MODO CAUSAL (C++ Batch)
# Aqui usamos denoise_signal_causal, que manda o vetor todo pro C++
# O C++ simula a chegada sequencial (Ring Buffer) mas roda compilado.
cat("\n2. Executando Modo CAUSAL TURBO (C++ Loop Interno)...\n")
time_caus_cpp = system.time({
  res_caus = apply(data_matrix, 2, function(x) {
    denoise_signal_causal(
      x, SCHEME, levels = LEVELS,
      window_size = 256, method="semisoft"
      )
  })
})
cat(sprintf("   Tempo: %.4f s\n", time_caus_cpp["elapsed"]))


# MODO CAUSAL NAIVE (R Loop - Apenas 10 curvas para não travar)
# Apenas para demonstrar a diferença de arquitetura
cat("\n3. Executando Modo CAUSAL NAIVE (R Loop - Amostra de 10 curvas)...")
cat("   (Estimando tempo para 1000 curvas...)\n")

# Cria processador
proc = new_wavelet_stream(SCHEME, window_size = 256, levels = LEVELS)

time_naive_sample = system.time({
  # Processa apenas as primeiras 10 colunas
  for(j in 1:10) {
    signal = data_matrix[, j]
    # Simula chegada ponto a ponto no R
    out = numeric(T_LEN)
    for(i in 1:T_LEN) {
      out[i] = proc(signal[i])
    }
  }
})

# Projeção linear para 1000 curvas
est_time_naive = time_naive_sample["elapsed"] * (N_CURVES / 10)
cat(
  sprintf(
    "   Tempo (Estimado para %d): %.4f s\n", N_CURVES,
    est_time_naive
    )
  )


# RELATÓRIO FINAL
cat("\n=== RESUMO DE PERFORMANCE (1.000 CURVAS) ===\n")
fps_off  = N_CURVES / time_off["elapsed"]
fps_caus = N_CURVES / time_caus_cpp["elapsed"]
fps_naive = N_CURVES / est_time_naive

cat(
  sprintf(
    "%-20s | %-10s | %-15s\n", "MÉTODO", "TEMPO TOTAL", "CURVAS/SEG"
    )
  )
cat("----------------------------------------------------------\n")
cat(
  sprintf(
    "%-20s | %7.4f s   | %9.1f cps\n", "Offline (Global)",
    time_off["elapsed"],
    fps_off
    )
  )

cat(
  sprintf(
    "%-20s | %7.4f s   | %9.1f cps\n", "Causal (C++ Turbo)",
    time_caus_cpp["elapsed"],
    fps_caus
    )
  )

cat(
  sprintf(
    "%-20s | %7.4f s   | %9.1f cps (Est)\n", "Causal (R Loop)",
    est_time_naive,
    fps_naive
    )
  )

ratio = fps_caus / fps_naive
cat(sprintf("\n>>> Otimização C++ vs R Loop: %.1fx mais rápido\n", ratio))
cat(
  sprintf(
    ">>> Custo da Causalidade (vs Offline): %.1f%% mais lento (Esperado)\n",
    (1 - fps_caus/fps_off)*100)
  )
