# inst/scripts/denoising_performance.R
if(!require("ggplot2")) install.packages("ggplot2")
if(!require("tidyr"))   install.packages("tidyr")
if(!require("dplyr"))   install.packages("dplyr")
if(!require("gridExtra")) install.packages("gridExtra")

devtools::load_all(quiet = TRUE)

# FUNCOES AUXILIARES (Metricas)
calc_metrics = function(clean, denoised) {
  err = clean - denoised
  rmse = sqrt(mean(err^2))

  pow_sig = sum(clean^2)
  pow_err = sum(err^2)
  snr = if(pow_err < 1e-15) 99.0 else 10 * log10(pow_sig / pow_err)

  list(rmse = rmse, snr = snr)
}

# ENGINE DE VISUALIZACAO E TESTE
run_visual_benchmark = function(
    signal_name,
    wavelet_name,
    ext = "symmetric"
    ) {

  # 1. Setup
  clean = rLifting:::.generate_signal(signal_name, n = 512)

  # SNR de Entrada ~15dB
  sigma = sqrt(mean(clean^2) / (10^(15/10)))
  set.seed(42)
  noisy = clean + rnorm(length(clean), 0, sigma)

  metrics_in = calc_metrics(clean, noisy)

  # 2. Processamento
  sch = lifting_scheme(wavelet_name)
  # Parametros do artigo
  alpha = 0.3; beta = 1.2; W = 128
  L = 3 # Niveis de decomposicao padronizados

  # A. Offline
  out_off = denoise_signal_offline(
    noisy, sch,
    levels = L,
    alpha = alpha, beta = beta,
    method = "semisoft", extension = ext
  )
  metrics_off = calc_metrics(clean, out_off)

  # B. Causal
  out_causal = denoise_signal_causal(
    noisy, sch,
    levels = L,
    window_size = W,
    alpha = alpha, beta = beta,
    method = "semisoft", extension = ext
  )
  metrics_causal = calc_metrics(clean, out_causal)

  # 3. Consolidacao de Dados para Plot
  df = data.frame(
    Index = 1:length(clean),
    Clean = clean,
    Noisy = noisy,
    Offline = out_off,
    Causal = out_causal
  )

  df_long = df |>
    pivot_longer(
      cols = -Index, names_to = "Type", values_to = "Value"
    ) |>
    mutate(
      Type = factor(Type, levels = c("Clean", "Noisy", "Offline", "Causal"))
    )

  # Titulo com Metricas
  title_str = sprintf(
    "%s + %s (L=%d, %s)\nIn: %.1fdB | Off: %.1fdB | Causal: %.1fdB",
    toupper(signal_name), toupper(wavelet_name), L, toupper(ext),
    metrics_in$snr, metrics_off$snr, metrics_causal$snr
  )

  # 4. Plot
  p = ggplot(
    df_long, aes(x = Index, y = Value, color = Type, alpha = Type)
    ) +
    geom_line() +
    scale_color_manual(values = c(
      "Clean" = "black",
      "Noisy" = "grey70",
      "Offline" = "blue",
      "Causal" = "red"
    )) +
    scale_alpha_manual(values = c(1, 0.4, 0.8, 0.8)) +
    theme_bw() +
    labs(title = title_str, x = NULL, y = NULL) +
    theme(legend.position = "bottom")

  return(
    list(plot = p, metrics = list(off = metrics_off, causal = metrics_causal))
  )
}

# MATRIZ DE TESTES
signals  = c("doppler", "heavisine", "bumps")
wavelets = c("haar", "db2", "cdf97")
extensions = c("symmetric", "zero")

plot_list = list()
idx = 1

cat(sprintf("%-10s | %-8s | %-8s | %-9s | %-9s\n",
            "SINAL", "WAVELET", "EXT", "SNR(Off)", "SNR(Caus)"))
cat(paste(rep("-", 60), collapse=""), "\n")

for (sig in signals) {
  for (wav in wavelets) {
    for (ext in extensions) {

      res = run_visual_benchmark(sig, wav, ext)

      # Report Tabular
      cat(
        sprintf(
          "%-10s | %-8s | %-8s | %5.2f dB  | %5.2f dB\n",
          sig, toupper(wav), ext,
          res$metrics$off$snr, res$metrics$causal$snr
        )
      )

      if (ext == "symmetric") {
        plot_list[[idx]] = res$plot
        idx = idx + 1
      }
    }
  }
  cat(paste(rep("-", 60), collapse=""), "\n")
}

grid.arrange(grobs = plot_list, ncol = 3)
