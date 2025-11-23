devtools::load_all(quiet = TRUE)

# Gera um vetor linear de -10 a 10 (simulando coeficientes)
x = seq(-10, 10, length.out = 100)
lambda = 3

y_hard = threshold(x, lambda, "hard")
y_soft = threshold(x, lambda, "soft")
y_semi = threshold(x, lambda, "semisoft")

# Visualizacao tabular compacta (focada na transicao positiva)
# Vamos olhar a regiao entre 0 e 5 (onde lambda=3 atua)
idx_focus = which(x >= 0 & x <= 6)
df = data.frame(
  Input = x[idx_focus],
  Hard  = y_hard[idx_focus],
  Soft  = y_soft[idx_focus],
  Semi  = y_semi[idx_focus]
)

# Filtra linhas para nao poluir o console
print(df[seq(1, nrow(df), length.out = 15), ], digits = 2)

# Check de Sanidade
cat("\n--- Verificacao de Comportamento (Lambda = 3) ---\n")

# 1. Abaixo do limiar (Input = 2)
val_2 = threshold(2, 3, "hard")
cat(
  sprintf(
    "Input 2 (< 3): Hard=%.1f (Esp: 0.0) -> %s\n",
    val_2, if(val_2==0) "OK" else "FAIL"
    )
  )

# 2. No limiar (Input = 3)
# Hard mantem 3? Soft vira 0? Semi vira 0?
cat(
  sprintf(
    "Input 3 (= 3): Hard=%.1f, Soft=%.1f, Semi=%.1f\n",
    threshold(3, 3, "hard"),
    threshold(3, 3, "soft"),
    threshold(3, 3, "semisoft")
    )
  )

# 3. Acima do limiar (Input = 5)
# Hard = 5. Soft = 5-3=2. Semi = sqrt(25-9)=4.
cat(
  sprintf(
    "Input 5 (> 3): Hard=%.1f, Soft=%.1f, Semi=%.1f (Esp: 4.0)\n",
    threshold(5, 3, "hard"),
    threshold(5, 3, "soft"),
    threshold(5, 3, "semisoft")
    )
  )
