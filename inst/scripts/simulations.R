# Script que executa o benchmark estatístico
# por velocidade, rodamos execuções em paralelo

# Setup
if(!require("dplyr")) install.packages("dplyr")
if(!require("parallel")) install.packages("parallel")

library(rLifting)
library(parallel)
library(dplyr)

N_SIMS = 1000 # total de curvas por sinal
T_POINTS = 2048 # pontos por curva
SIGMA_NOISE = 0.2 # desvio padrão do ruído
OUTPUT_FILE = "vignettes/sim_results_N1000.rds"

# Frequência de atualização do thresholding causal (para ganhar velocidade)
UPDATE_FREQ_CAUSAL = 16

# Paralelismo
NUM_CORES = parallel::detectCores(logical = FALSE) - 1
if (NUM_CORES < 1) NUM_CORES = 1

# Tamanho do bloco para atualização da barra de progresso
CHUNK_SIZE = 100

set.seed(42)

# Funções auxiliares

calc_metrics = function(original, processado) {
  erro = original - processado
  mse = mean(erro^2)
  rmse = sqrt(mse)

  pot_sinal = sum(original^2)
  pot_ruido = sum(erro^2)
  if(pot_ruido < 1e-15) pot_ruido = 1e-15
  snr = 10 * log10(pot_sinal / pot_ruido)
  mbe = mean(erro)

  return(c(RMSE = rmse, SNR = snr, MBE = mbe))
}

run_single_simulation = function(
    sinal_ruidoso, sinal_puro,
    grid, update_freq_val
    ) {
  n_grid = nrow(grid)
  matriz_res = matrix(NA, nrow = n_grid + 1, ncol = 3)

  matriz_res[1, ] = calc_metrics(sinal_puro, sinal_ruidoso)

  for (i in 1:n_grid) {
    p = grid[i, ]
    sch = lifting_scheme(p$Wavelet)

    if (p$Mode == "Offline") {
      proc = denoise_signal_offline(
        sinal_ruidoso, sch,
        method = p$Method, extension = p$Extension
      )
    } else {
      proc = denoise_signal_causal(
        sinal_ruidoso, sch,
        window_size = 256,
        method = p$Method,
        extension = p$Extension,
        update_freq = update_freq_val
      )
    }
    matriz_res[i + 1, ] = calc_metrics(sinal_puro, proc)
  }
  return(matriz_res)
}

# Definição de grid

grid_params = expand.grid(
  Wavelet = c("haar", "db2", "cdf53", "cdf97", "dd4"),
  Mode = c("Offline", "Causal"),
  Method = c("hard", "soft", "semisoft"),
  Extension = c("symmetric", "periodic"),
  stringsAsFactors = FALSE
)

cl = makeCluster(NUM_CORES)
clusterEvalQ(cl, { library(rLifting) })
clusterExport(
  cl, c(
    "run_single_simulation", "calc_metrics",
    "grid_params", "UPDATE_FREQ_CAUSAL"
    )
  )

# Execução
cat(
  sprintf(
    "INICIANDO SIMULAÇÃO (N=%d, Cores=%d, Stride=%d)\n",
    N_SIMS, NUM_CORES, UPDATE_FREQ_CAUSAL
    )
  )

resultados_globais = list()
tipos_sinal = c("bumps", "doppler", "heavisine")

for (tipo in tipos_sinal) {
  cat(sprintf("\nProcessando: %s\n", toupper(tipo)))

  # Gerando dados em lote
  sinal_puro = rLifting:::.generate_signal(tipo, n = T_POINTS)
  matriz_ruido = matrix(
    rnorm(T_POINTS * N_SIMS, sd = SIGMA_NOISE),
    ncol = N_SIMS
    )
  matriz_sinais = matriz_ruido + sinal_puro

  # Exportar sinal puro para workers (específico deste loop)
  clusterExport(cl, "sinal_puro", envir = environment())

  # Processamento em chunks
  indices = 1:N_SIMS
  chunks = split(indices, ceiling(seq_along(indices)/CHUNK_SIZE))

  lista_res_parcial = list()

  pb = txtProgressBar(min = 0, max = N_SIMS, style = 3)

  for (chunk_idx in seq_along(chunks)) {
    idx_batch = chunks[[chunk_idx]]

    lista_inputs = as.list(
      as.data.frame(
        matriz_sinais[, idx_batch, drop=FALSE]
        )
      )

    res_chunk = parLapply(cl, lista_inputs, function(col) {
      run_single_simulation(col, sinal_puro, grid_params, UPDATE_FREQ_CAUSAL)
    })

    lista_res_parcial = c(lista_res_parcial, res_chunk)
    setTxtProgressBar(pb, max(idx_batch))
  }
  close(pb)

  # Agregação
  soma_metricas = Reduce(`+`, lista_res_parcial)
  media_metricas = soma_metricas / N_SIMS

  df_base = data.frame(
    Sinal = tipo, Wavelet = "None", Mode = "Baseline",
    Method = "None", Extension = "None",
    RMSE = media_metricas[1, 1],
    SNR = media_metricas[1, 2],
    MBE = media_metricas[1, 3]
  )

  df_grid = grid_params
  df_grid$Sinal = tipo
  df_grid$RMSE = media_metricas[2:nrow(media_metricas), 1]
  df_grid$SNR  = media_metricas[2:nrow(media_metricas), 2]
  df_grid$MBE  = media_metricas[2:nrow(media_metricas), 3]

  resultados_globais[[tipo]] = rbind(df_base, df_grid)
}

stopCluster(cl)

# Salvar
df_final = dplyr::bind_rows(resultados_globais)
saveRDS(df_final, file = OUTPUT_FILE)

cat(sprintf("\nConcluído! Resultados salvos em: %s\n", OUTPUT_FILE))
