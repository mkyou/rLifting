devtools::load_all(quiet = TRUE)

# Função auxiliar para visualização gráfica (texto)
visualize_pad = function(input, n_left, n_right, mode) {

  # Acessa a funcao interna (privada) do pacote
  output = rLifting:::.pad_signal(input, n_left, n_right, mode)

  n_in = length(input)
  n_out = length(output)

  cat(
    sprintf(
      "\n>>> MODO: %-10s | Pad Left: %d | Pad Right: %d <<<\n",
      toupper(mode), n_left, n_right
      )
    )

  # Cria marcadores visuais
  # [L] = Padding Esquerdo
  # [S] = Sinal Original
  # [R] = Padding Direito

  types = c(rep("[L]", n_left), rep("[S]", n_in), rep("[R]", n_right))

  # Formata para visualizacao alinhada
  df = data.frame(
    Index = 1:n_out,
    Region = types,
    Value = output
  )

  # Imprime em formato horizontal para facil leitura
  cat("Regiao: ", paste(format(df$Region, width=4), collapse=" "), "\n")
  cat("Valor:  ", paste(format(df$Value, width=4), collapse=" "), "\n")

  # Checagem de integridade (O miolo foi preservado?)
  miolo = output[(n_left + 1) : (n_left + n_in)]
  if (all(miolo == input)) {
    cat("Status:  INTEGRIDADE OK (Sinal original preservado no centro)\n")
  } else {
    cat("Status:  ERRO CRITICO (Sinal original corrompido!)\n")
  }
}

# BATERIA DE TESTES
#-------------------------------------------------------------
# Sinal simples para rastreio: 10, 20, 30
x = c(10, 20, 30)

cat("SINAL ORIGINAL: 10 20 30\n")

# TESTE 1: SYMMETRIC (Reflexao / Espelho)
# Logica esperada (Half-Point): Reflete repetindo a borda.
# Esquerda (2): 20, 10 ...
# Direita (2): ... 30, 20
visualize_pad(x, 2, 2, "symmetric")

# TESTE 2: PERIODIC (Circular / Modulo)
# Logica esperada: O sinal da volta.
# Esquerda (2): Pega o fim (20, 30)
# Direita (2): Pega o inicio (10, 20)
visualize_pad(x, 2, 2, "periodic")

# TESTE 3: ZERO (Padding com 0)
visualize_pad(x, 2, 2, "zero")

#-------------------------------------------------------------
# TESTE 4: STRESS TEST (Padding > Tamanho do Sinal)
cat("\nSTRESS TEST (Pad=5 > Len=3)\n")

# Symmetric Stress:
# Deve refletir e, se acabar o sinal, truncar ou repetir a ponta
# (depende da implementacao).
# Nossa implementacao faz clamp na ultima borda se exceder.
visualize_pad(x, 5, 0, "symmetric")

# Periodic Stress:
# Deve dar varias voltas no circulo.
# Pad 5 (Esq) num sinal de 3:
# Volta 1 (3 amostras): 10 20 30
# Volta 2 (2 amostras): 20 30
# Esperado na esq: 20 30 10 20 30
visualize_pad(x, 5, 0, "periodic")
