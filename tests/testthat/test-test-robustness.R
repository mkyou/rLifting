test_that("Stream processor handles NA and Inf gracefully", {
  sch = lifting_scheme("haar")
  # Janela de 16. As primeiras 15 chamadas serao BYPASS (retorna o input).
  # A 16a chamada processa.
  proc = new_wavelet_stream(sch, window_size = 16)

  # Envia valores normais (Bypass phase)
  v1 = proc(1)
  expect_equal(v1, 1) # Deve retornar 1 puro

  # Envia NA (deve retornar NA e emitir warning)
  expect_warning(val_na <- proc(NA))
  expect_true(is.na(val_na))

  # Envia Inf (deve retornar Inf e emitir warning)
  expect_warning(val_inf <- proc(Inf))
  expect_true(is.infinite(val_inf))
})

test_that("Stream processor rejects bad window sizes", {
  sch = lifting_scheme("haar")
  # O minimo absoluto agora e 8 (para garantir que quando rodar, rode bem)
  expect_error(new_wavelet_stream(sch, window_size = 7))
})
