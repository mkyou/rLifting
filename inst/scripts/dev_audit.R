devtools::load_all()

# HAAR
diagnose_wavelet("haar", list(
  is_ortho = TRUE,
  vm_degrees = c(0),
  max_taps = 2
))

# DB2
diagnose_wavelet("db2", list(
  is_ortho = TRUE,
  vm_degrees = c(0, 1),
  max_taps = 4
))

# CDF 9/7
diagnose_wavelet("cdf97", list(
  is_ortho = FALSE,
  vm_degrees = c(0, 1, 2, 3),
  max_taps = 11
))

# CDF 5/3
diagnose_wavelet("cdf53", list(
  is_ortho = FALSE,
  vm_degrees = c(0, 1, 2, 3),
  max_taps = 11
))

# DD4
diagnose_wavelet("dd4", list(
  is_ortho = FALSE,
  vm_degrees = c(0, 1, 2, 3),
  max_taps = 11
))

# Lazy
diagnose_wavelet("lazy", list(
  is_ortho = FALSE,
  vm_degrees = c(0, 1, 2, 3),
  max_taps = 11
))


# Predict: (0.5, 0.5). Position 'center' p/ length 2 assume idx 0. Correto.
p_53 = lift_step("predict", c(0.5, 0.5), position = "center")

# Update: (0.25, 0.25). Queremos acessar k-1 e k.
# Position 'left' calcula start_idx = -n + 1 = -2 + 1 = -1. Correto.
u_53 = lift_step("update", c(0.25, 0.25), position = "left")

minha_cdf53 = lifting_scheme(
  wavelet = "MyCDF5.3",
  custom_steps = list(p_53, u_53),
  custom_norm = c(sqrt(2), 1/sqrt(2))
)

# Diagnostico (espera-se nao-ortogonal, VM grau 0 e 1 ok, grau 2 falha)
diagnose_wavelet(minha_cdf53, list(
  is_ortho = FALSE,
  vm_degrees = c(0, 1, 2),
  max_taps = 5
))
