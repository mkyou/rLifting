#' Extensao Simetrica (Half-Point Symmetric Extension)
#'
#' Reflete os dados nas bordas para mitigar efeitos de contorno.
#'
#' @param x Vetor de entrada.
#' @param n_left Quantidade de amostras para estender a esquerda.
#' @param n_right Quantidade de amostras para estender a direita.
#'
#' @keywords internal
.sym_extension = function(x, n_left, n_right) {
  n = length(x)

  # Previne erros se a extensao for maior que o sinal
  # (Nesse caso, repete-se a logica de reflexao, mas simplificamos aqui)
  if (n_left > n) n_left = n
  if (n_right > n) n_right = n

  left_pad = numeric(0)
  right_pad = numeric(0)

  if (n_left > 0) left_pad = x[seq(n_left, 1, by = -1)]
  if (n_right > 0) right_pad = x[seq(n, n - n_right + 1, by = -1)]

  c(left_pad, x, right_pad)
}

#' Aplica filtro com tratamento de borda e shift
#'
#' Wrapper que gerencia padding e alinhamento temporal (shift).
#'
#' @param x Sinal de entrada.
#' @param coeffs Coeficientes do filtro.
#' @param start_idx Indice onde o filtro comeca relativo a posicao atual.
#'                  Ex: 0 = alinhado, -1 = olha 1 para tras.
#'
#' @return Sinal filtrado alinhado.
#' @keywords internal
.apply_filter_lifting = function(x, coeffs, start_idx = 0) {
  n_x = length(x)
  n_c = length(coeffs)

  # 1. Definir a janela de leitura necessária
  # Para calcular y[1], precisamos de x[1 + start_idx] até
  # x[1 + start_idx + n_c - 1]
  # O menor índice acessado será (1 + start_idx)
  # O maior índice acessado será (n_x + start_idx + n_c - 1)

  min_idx = 1 + start_idx
  max_idx = n_x + start_idx + n_c - 1

  # 2. Calcular quanto padding precisamos
  pad_left = 0
  pad_right = 0

  if (min_idx < 1) pad_left = 1 - min_idx
  if (max_idx > n_x) pad_right = max_idx - n_x

  # Adicionamos a extensão simétrica
  x_padded = .sym_extension(x, pad_left, pad_right)

  # 3. Convolução Manual Vetorizada
  # y inicializado com zeros
  y = numeric(n_x)

  # O indice '1' de x_padded corresponde ao índice '1 - pad_left' do original.
  # Queremos acessar x_original[i + start_idx + k].
  # No padded, isso equivale a x_padded[i + start_idx + k + pad_left].

  base_idx = (1:n_x) + start_idx + pad_left

  for (k in 1:n_c) {
    # Pega a fatia do vetor deslocada para o k-ésimo coeficiente
    # Indices para o coeficiente k (0-based no loop seria k-1, aqui 1-based é k)
    # A fórmula é: soma( x[pos + (k-1)] * coeff[k] )

    idx = base_idx + (k - 1)
    y = y + x_padded[idx] * coeffs[k]
  }

  return(y)
}
