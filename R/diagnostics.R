#' Gera sinais padrao para testes de ondaletas (Benchmark Donoho-Johnstone)
#'
#' Cria sinais sinteticos classicos para validacao de wavelets.
#' As formulas sao baseadas em Donoho & Johnstone e no artigo de Liu et al.
#'
#' @param type Tipo do sinal: "const", "ramp", "poly2", "poly3", "random",
#'             "impulse", "sine", "doppler", "heavisine", "bumps".
#' @param n Tamanho do sinal (padrao 512).
#'
#' @return Vetor numerico.
#' @keywords internal
.generate_signal = function(type, n = 512) {
  t = seq(0, 1, length.out = n)

  # Basicos
  if (type == "const")   return(rep(100, n))
  if (type == "ramp")    return(seq(1, 100, length.out = n))
  if (type == "poly2")   return(100 * t^2)
  if (type == "poly3")   return(100 * t^3)
  if (type == "random")  return(rnorm(n))
  if (type == "impulse") {
    x = rep(0, n)
    x[floor(n/2)] = 1
    return(x)
  }
  if (type == "sine")    return(sin(4 * pi * t))

  # Benchmarks Avancados (Ref: Liu et al., Appendix B)

  # Doppler: Frequencia variavel
  if (type == "doppler") {
    eps = 0.05
    return(sqrt(t * (1 - t)) * sin((2 * pi * 1.05) / (t + eps)))
  }

  # HeaviSine: Senoide com saltos (descontinuidades)
  if (type == "heavisine") {
    return(4 * sin(4 * pi * t) - sign(t - 0.3) - sign(0.72 - t))
  }

  # Bumps: Picos localizados
  if (type == "bumps") {
    pos = c(0.1, 0.13, 0.15, 0.23, 0.25, 0.40, 0.44, 0.65, 0.76, 0.78, 0.81)
    h   = c(4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 5.1, -4.2)
    w   = c(0.005, 0.005, 0.006, 0.01, 0.01, 0.03, 0.01, 0.01, 0.005,
            0.008, 0.005)

    x = numeric(n)
    for (j in seq_along(pos)) {
      x = x + h[j] * (1 + abs((t - pos[j]) / w[j])^4)^(-1)
    }
    return(x)
  }

  stop(sprintf("Tipo de sinal '%s' desconhecido.", type))
}

#' Valida a Reconstrucao Perfeita (Stress Test)
#'
#' Verifica a invertibilidade da ondaleta contra uma bateria de sinais.
#'
#' @param scheme Objeto da classe `lifting_scheme`.
#' @param tol Tolerancia de erro numerico (padrao 1e-10).
#'
#' @return Lista com status global e erro maximo encontrado.
#' @export
validate_perfect_reconstruction = function(scheme, tol = 1e-10) {
  signals = c("random", "ramp", "sine", "doppler", "heavisine", "bumps")
  max_err_global = 0
  failed_signals = c()

  for (sig_type in signals) {
    x = .generate_signal(sig_type, n = 512)

    # Usa periodic para focar na matematica e nao na borda
    res = lwt(x, scheme, extension = "periodic")
    rec = ilwt(res)

    err = max(abs(x - rec))
    max_err_global = max(max_err_global, err)

    if (err >= tol) {
      failed_signals = c(failed_signals, sig_type)
    }
  }

  passed = length(failed_signals) == 0
  msg = if (passed) {
    sprintf("Erro Max (Global): %.2e", max_err_global)
  } else {
    sprintf("FALHA em: %s. Erro Max: %.2e",
            paste(failed_signals, collapse=", "), max_err_global)
  }

  list(
    name = "Perfect Reconstruction (Stress Test)",
    passed = passed,
    metric = max_err_global,
    msg = msg
  )
}

#' Valida momentos nulos (Vanishing Moments)
#'
#' Verifica se a ondaleta anula polinomios de grau especifico.
#'
#' @param scheme Objeto da classe `lifting_scheme`.
#' @param degree Grau do polinomio (0=Constante, 1=Rampa, 2=Parab...).
#' @param tol Tolerancia de energia residual (padrao 1e-9).
#'
#' @return Lista com status e energia residual.
#' @export
validate_vanishing_moments = function(scheme, degree = 0, tol = 1e-9) {
  n = 128
  type_map = c("0"="const", "1"="ramp", "2"="poly2", "3"="poly3")
  sig_type = type_map[as.character(degree)]

  if (is.na(sig_type)) return(
    list(name = "VM", passed = FALSE, msg = "Grau invalido")
  )

  x = .generate_signal(sig_type, n)

  # Periodic para evitar vazamento de energia nas bordas
  res = lwt(x, scheme, extension = "periodic")
  d1 = res$coeffs$d1

  # Corte agressivo do miolo para garantir que estamos testando o filtro
  cut = floor(n * 0.15)
  d_core = d1[(cut+1):(length(d1)-cut)]

  energy = sum(d_core^2)
  passed = energy < tol

  list(
    name = sprintf("Vanishing Moments (Degree %d)", degree),
    passed = passed,
    metric = energy,
    msg = sprintf("Energia Residual: %.2e (Tol: %.2e)", energy, tol)
  )
}

#' Valida Ortogonalidade (Conservacao de Energia)
#'
#' Verifica o Teorema de Parseval. Apenas para ondaletas ortogonais.
#'
#' @param scheme Objeto da classe `lifting_scheme`.
#' @param expected Booleano. Se TRUE, exige ortogonalidade.
#' @param tol Tolerancia (padrao 1e-10).
#'
#' @return Lista com status e razao de energia (Out/In).
#' @export
validate_orthogonality = function(scheme, expected = TRUE, tol = 1e-10) {
  x = .generate_signal("random", n = 512)
  energy_in = sum(x^2)

  res = lwt(x, scheme, extension = "periodic")
  # Energia total nos coeficientes (Aprox + Detalhes)
  energy_out = sum(res$coeffs$a1^2) + sum(res$coeffs$d1^2)

  ratio = energy_out / energy_in
  is_ortho = abs(ratio - 1) < tol

  passed = (expected && is_ortho) || (!expected)
  if (expected && !is_ortho) passed = FALSE

  status_str = if (is_ortho) "Ortogonal" else "Nao-Ortogonal"

  list(
    name = "Orthogonality (Energy Conservation)",
    passed = passed,
    metric = ratio,
    msg = sprintf("Ratio: %.6f (%s) [Esp: %s]",
                  ratio, status_str, expected)
  )
}

#' Valida Suporte Compacto (FIR Compliance)
#'
#' Verifica se a resposta ao impulso e finita (Filtro FIR).
#'
#' @param scheme Objeto da classe `lifting_scheme`.
#' @param max_width Largura maxima esperada (numero de taps).
#'
#' @return Lista com status e numero de taps ativos.
#' @export
validate_compact_support = function(scheme, max_width) {
  n = 256
  coeffs_mock = list(a1 = rep(0, n/2), d1 = rep(0, n/2))
  # Impulso no detalhe, longe das bordas
  coeffs_mock$d1[n/4] = 1

  mock_lwt = structure(
    list(
      coeffs = coeffs_mock, scheme = scheme,
      levels = 1, original_len = n, extension = "zero"
    ),
    class = "lwt"
  )

  psi = ilwt(mock_lwt)

  # Conta coeficientes nao-nulos (acima de ruido numerico)
  nonzero = sum(abs(psi) > 1e-10)

  # Margem de tolerancia de +2 taps para variacoes de fase
  passed = nonzero <= (max_width + 2) && nonzero >= 1

  list(
    name = "Compact Support (FIR)",
    passed = passed,
    metric = nonzero,
    msg = sprintf("Active Taps: %d (Max Expected: %d)", nonzero, max_width)
  )
}

#' Valida Sensibilidade a Deslocamento (Shift Variance)
#'
#' Ondaletas decimadas nao sao invariantes a translacao.
#' Este teste quantifica a variacao de energia nos detalhes ao
#' deslocar o sinal de entrada em 1 amostra.
#'
#' @param scheme Objeto da classe `lifting_scheme`.
#'
#' @return Lista com status e variacao percentual.
#' @export
validate_shift_sensitivity = function(scheme) {
  # Sinal teste: Impulso centralizado
  n = 128
  x = rep(0, n)
  x[n/2] = 1

  # Energia original no nivel 1
  res1 = lwt(x, scheme, levels = 1, extension = "periodic")
  energy1 = sum(res1$coeffs$d1^2)

  # Sinal deslocado circularmente em 1 amostra
  x_shift = c(x[n], x[1:(n-1)])
  res2 = lwt(x_shift, scheme, levels = 1, extension = "periodic")
  energy2 = sum(res2$coeffs$d1^2)

  # Diferenca percentual
  # Se a energia for zero (sinal constante), diff e zero
  diff_pct = 0
  if (energy1 > 1e-15) {
    diff_pct = abs(energy1 - energy2) / energy1 * 100
  }

  # O teste "passa" sempre, pois shift-variance e esperado em DWT.
  # A metrica serve apenas de informacao.
  msg = sprintf("Variacao Energia (Shift 1): %.2f%%", diff_pct)

  list(
    name = "Shift Sensitivity",
    passed = TRUE,
    metric = diff_pct,
    msg = msg
  )
}

#' Visualiza as funcoes de base (Scaling e Wavelet)
#'
#' Plota a forma da onda iterando a reconstrucao por varios niveis.
#'
#' @param scheme Objeto `lifting_scheme`.
#' @param plot Booleano.
#' @param levels Numero de niveis de cascata.
#'
#' @export
visualize_wavelet_basis = function(scheme, plot = TRUE, levels = 8) {
  n = 2^levels

  # 1. Função Wavelet (Psi)
  coeffs_psi = list()
  for (j in 1:levels) coeffs_psi[[paste0("d", j)]] = numeric(n / 2^j)
  coeffs_psi[[paste0("a", levels)]] = numeric(n / 2^levels)
  coeffs_psi[[paste0("d", levels)]][1] = 1

  lwt_psi = structure(
    list(coeffs = coeffs_psi, scheme = scheme, levels = levels,
         original_len = n, extension = "zero"),
    class = "lwt"
  )
  psi = ilwt(lwt_psi)

  # 2. Função Scaling (Phi)
  coeffs_phi = list()
  for (j in 1:levels) coeffs_phi[[paste0("d", j)]] = numeric(n / 2^j)
  coeffs_phi[[paste0("a", levels)]] = numeric(n / 2^levels)
  coeffs_phi[[paste0("a", levels)]][1] = 1

  lwt_phi = structure(
    list(coeffs = coeffs_phi, scheme = scheme, levels = levels,
         original_len = n, extension = "zero"),
    class = "lwt"
  )
  phi = ilwt(lwt_phi)

  if (plot) {
    par(mfrow = c(1, 2))
    ts.plot(phi, main = paste("Scaling (Phi):", scheme$wavelet),
            ylab = "Amp", col = "blue", lwd = 1.5)
    grid()
    ts.plot(psi, main = paste("Wavelet (Psi):", scheme$wavelet),
            ylab = "Amp", col = "red", lwd = 1.5)
    grid()
    par(mfrow = c(1, 1))
  }
}

#' @return (Invisivel) Uma lista contendo os resultados de cada teste:
#' \itemize{
#'   \item Perfect Reconstruction
#'   \item Orthogonality
#'   \item Vanishing Moments (para cada grau)
#'   \item Compact Support
#'   \item Shift Sensitivity
#' }
#' Cada item possui: \code{$passed} (TRUE/FALSE), \code{$metric}
#' (valor numerico) e \code{$msg}.
diagnose_wavelet = function(wavelet_name, config, verbose = TRUE) {

  if (is.character(wavelet_name)) {
    sch = tryCatch(lifting_scheme(wavelet_name), error = function(e) NULL)
  } else {
    sch = wavelet_name
    wavelet_name = sch$wavelet
  }

  if (is.null(sch)) stop("Ondaleta invalida.")

  if (verbose) {
    cat(sprintf("\n=== DIAGNOSTICO: %s ===\n", toupper(wavelet_name)))
  }

  tests = list()
  tests[[1]] = validate_perfect_reconstruction(sch)
  tests[[2]] = validate_orthogonality(sch, expected = config$is_ortho)

  for (deg in config$vm_degrees) {
    idx = length(tests) + 1
    tests[[idx]] = validate_vanishing_moments(sch, degree = deg)
  }

  tests[[length(tests)+1]] = validate_compact_support(
    sch, max_width = config$max_taps
  )

  # Adicionado teste de Shift Sensitivity
  tests[[length(tests)+1]] = validate_shift_sensitivity(sch)

  if (verbose) {
    cat("Gerando visualizacao das bases (verifique a janela de plot)...\n")
    try(visualize_wavelet_basis(sch, plot = TRUE))
  }

  if (verbose) {
    for (res in tests) {
      status = if(res$passed) "[PASS]" else "[FAIL]"
      # Se for shift sensitivity, status e informativo
      if (res$name == "Shift Sensitivity") status = "[INFO]"

      cat(sprintf("%-6s %-35s | %s\n", status, res$name, res$msg))
    }
  }
  invisible(tests)
}
