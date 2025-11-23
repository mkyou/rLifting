#' Tratamento de padding/extensao de sinais
#'
#' @param x Sinal de entrada.
#' @param n_left Amostras a esquerda.
#' @param n_right Amostras a direita.
#' @param mode Modo de extensao: "symmetric", "periodic", "zero".
#'
#' @keywords internal
.pad_signal = function(x, n_left, n_right, mode = "symmetric") {
  n = length(x)
  if (n == 0) return(x)

  left_pad = numeric(0)
  right_pad = numeric(0)

  if (mode == "zero") {
    if (n_left > 0)  left_pad  = rep(0, n_left)
    if (n_right > 0) right_pad = rep(0, n_right)

  } else if (mode == "periodic") {
    # Logica Circular Robusta (funciona mesmo se n_pad > n)
    if (n_left > 0) {
      # Indices: n, n-1, ... (circular)
      # Formula: (i - 1) %% n + 1
      idx = (seq(n - n_left + 1, n) - 1) %% n + 1
      # Ajuste: queremos os ultimos n_left elementos, em ordem
      start_idx = n - (n_left %% n)
      # Sequencia simples repetida
      full_idx = rep(1:n, ceiling(n_left/n) + 1)
      # Pega o trecho final correto
      left_pad = tail(full_idx, n_left)
      left_pad = x[left_pad]
    }
    if (n_right > 0) {
      idx = (0:(n_right-1)) %% n + 1
      right_pad = x[idx]
    }

  } else if (mode == "symmetric") {
    # Half-Point Symmetric (Espelhamento incluindo a borda: 3 2 1 | 1 2 3)
    # Robusto para n_pad > n (padrao "Triangular Wave")

    if (n_left > 0) {
      # Gera indices decrescentes a partir de 1: 1, 2, 3...
      # Reflete: 1, 2, ... n, n, n-1...
      # A logica mais simples para "Lifting" costuma ser apenas
      # refletir a borda imediata
      # Sequencia: 1, 1, 2, 3... ou 1, 2, 3...?
      # Padrao Lifting (Sweldens): x[1-k] refere a x[k] (half-point)

      # Sequencia bruta de indices espelhados: 1, 2, 3... n, n, n-1...
      # Para n_left, queremos os indices voltando de 1.

      # Implementacao simples: pmin(pmax(...)) funciona para 1 reflexao.
      # Para multiplas, precisamos de logica de onda triangular.
      # Como filtros wavelet sao curtos, raramente n_pad > n.
      # Vamos manter a robustez para 1 reflexao completa (n_pad <= n)
      # e truncar se for maior, para evitar complexidade desnecessaria no R.

      # Se n_left > n, truncamos para n (limite fisico de espelhamento unico)
      nl = min(n_left, n)
      idx = seq(nl, 1, by = -1)
      left_pad = x[idx]

      # Se ainda faltar (caso extremo n_pad > n),
      # preenche com a ultima borda (clamp)
      if (n_left > n) {
        extra = rep(x[n], n_left - n) # Repete o extremo oposto? Ou o proximo?
        # Decisao de Design: Wavelet filters > sinal = Borda degenerada.
        # Repetir o valor da borda mais distante (Clamp) e seguro.
        left_pad = c(extra, left_pad)
      }
    }

    if (n_right > 0) {
      nr = min(n_right, n)
      idx = seq(n, n - nr + 1, by = -1)
      right_pad = x[idx]

      if (n_right > n) {
        extra = rep(x[1], n_right - n)
        right_pad = c(right_pad, extra)
      }
    }
  }

  c(left_pad, x, right_pad)
}

#' Aplica filtro lifting com suporte a bordas
#'
#' @param x Sinal input.
#' @param coeffs Coeficientes do filtro P ou U.
#' @param start_idx Offset relativo.
#' @param extension Modo de extensao.
#'
#' @keywords internal
.apply_filter_lifting = function(
    x,
    coeffs,
    start_idx = 0,
    extension = "symmetric"
) {
  n_x = length(x)
  n_c = length(coeffs)

  # 1. Definir janela de leitura
  min_idx = 1 + start_idx
  max_idx = n_x + start_idx + n_c - 1

  # 2. Calcular padding necessario
  pad_left = 0
  pad_right = 0

  if (min_idx < 1) pad_left = 1 - min_idx
  if (max_idx > n_x) pad_right = max_idx - n_x

  # 3. Aplicar Padding
  x_padded = .pad_signal(x, pad_left, pad_right, mode = extension)

  # 4. Convolucao via Reduce (Substituindo loop for)
  # O indice 1 original agora esta em (1 + pad_left) no x_padded
  base_idx = (1:n_x) + start_idx + pad_left

  # Gera lista de vetores ponderados e soma tudo
  components = lapply(1:n_c, function(k) {
    idx = base_idx + (k - 1)
    x_padded[idx] * coeffs[k]
  })

  y = Reduce(`+`, components)

  return(y)
}
