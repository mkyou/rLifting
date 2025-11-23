#' Cria um processador de fluxo (stream) Wavelet Adaptativo
#'
#' Gera uma funcao stateful que implementa o metodo de Liu et al. completo:
#' Janela Movel + Decomposicao Lifting + Threshold Recursivo Adaptativo.
#'
#' @param scheme Objeto `lifting_scheme`.
#' @param window_size Tamanho da janela (W). Ex: 256.
#' @param alpha Parametro de decaimento do threshold (Eq 9).
#' @param beta Fator de ganho do threshold (Eq 9).
#' @param method Metodo de encolhimento: "hard", "soft", "semisoft".
#' @param extension Tratamento de borda ('symmetric', 'periodic').
#'
#' @return Uma funcao closure `processor(new_sample)` que aceita um unico
#' valor numerico e retorna o valor filtrado correspondente. O estado
#' (buffer) e mantido entre chamadas.
#' @export
new_wavelet_stream = function(
    scheme,
    window_size = 256,
    alpha = 0.3,
    beta = 1.2,
    method = "semisoft",
    extension = "symmetric"
    ) {

  buffer = numeric(0)

  processor = function(new_sample) {
    # 1. Atualiza buffer
    buffer <<- c(buffer, new_sample)
    n_curr = length(buffer)

    if (n_curr > window_size) {
      excess = n_curr - window_size
      buffer <<- buffer[(excess + 1):n_curr]
      n_curr = window_size
    }

    # Fase de aquecimento
    if (n_curr < window_size) {
      return(new_sample)
    }

    # 2. Forward LWT
    lwt_obj = lwt(buffer, scheme, levels = 1, extension = extension)

    # 3. Calculo dinamico do Threshold
    # O lambda muda a cada amostra dependendo do ruido atual da janela
    lambdas = compute_adaptive_threshold(lwt_obj, alpha, beta)

    # 4. Aplicar Threshold (suporta multi-nivel se necessario)
    for (lvl_name in names(lambdas)) {
      coeffs_old = lwt_obj$coeffs[[lvl_name]]
      lam = lambdas[[lvl_name]]

      # Aplica a funcao de encolhimento (soft/hard/semisoft)
      coeffs_new = threshold(coeffs_old, lam, method)

      lwt_obj$coeffs[[lvl_name]] = coeffs_new
    }

    # 5. Reconstrucao
    rec = ilwt(lwt_obj)

    # 6. Output Causal
    return(rec[length(rec)])
  }

  return(processor)
}

#' Denoising causal em lote (simulacao fiel)
#'
#' Processa um sinal completo simulando a chegada sequencial dos dados.
#' Garante que o resultado seja identico ao uso em tempo real.
#'
#' @param signal Vetor completo do sinal ruidoso.
#' @param scheme Objeto `lifting_scheme`.
#' @param window_size Tamanho da janela.
#' @param alpha Parametro de decaimento do threshold (Eq 9).
#' @param beta Fator de ganho do threshold (Eq 9).
#' @param method Metodo de thresholding.
#' @param extension Tratamento de borda ('symmetric', 'periodic').
#'
#' @return Vetor filtrado (mesmo comprimento).
#' @export
denoise_signal_causal = function(
    signal,
    scheme,
    window_size = 256,
    alpha = 0.3,
    beta = 1.2,
    method = "semisoft",
    extension = "symmetric"
    ) {

  n = length(signal)
  output = numeric(n)

  # Passa alpha/beta para o construtor
  stream_processor = new_wavelet_stream(
    scheme, window_size,
    alpha, beta, method, extension
    )

  pb = NULL
  if (n > 5000) {
    cat("Processando Batch Causal Adaptativo...\n")
    pb = txtProgressBar(min = 0, max = n, style = 3)
  }

  for (i in 1:n) {
    output[i] = stream_processor(signal[i])
    if (!is.null(pb) && i %% 100 == 0) setTxtProgressBar(pb, i)
  }

  if (!is.null(pb)) close(pb)
  return(output)
}
