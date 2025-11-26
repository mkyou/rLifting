#include "utils.h"
using namespace Rcpp;

// Helper: Calcula índice com tratamento de borda (Virtual Padding)
// Evita alocar memória para padding físico.
// i: indice desejado (pode ser negativo ou >= n)
// n: tamanho do sinal
// mode: 1=symmetric, 2=periodic, 3=zero
inline double get_val_at(const NumericVector& x, int i, int n, int mode) {

  // Dentro dos limites (Fast Path)
  if (i >= 0 && i < n) {
    return x[i];
  }

  // Zero Padding
  if (mode == 3) {
    return 0.0;
  }

  // Periodic (Circular)
  if (mode == 2) {
    // Ex: i=-1, n=10 -> (-1 % 10) = -1 + 10 = 9
    int idx = i % n;
    if (idx < 0) idx += n;
    return x[idx];
  }

  // Symmetric (Half-Point / Mirror)
  // Ex: [0 1 2]. i=-1 -> 0. i=-2 -> 1.
  if (mode == 1) {
    while (i < 0 || i >= n) {
      if (i < 0) {
        i = -1 - i; // Reflete no inicio (half-point)
      } else {
        i = 2 * n - 1 - i; // Reflete no fim
      }
    }
    // Clamp final caso a reflexao falhe (paranoia check)
    if (i < 0) i = 0;
    if (i >= n) i = n - 1;
    return x[i];
  }

  return 0.0;
}

//' Aplica filtro via convolucao (Core C++)
 //'
 //' @param x Sinal de entrada.
 //' @param coeffs Coeficientes do filtro.
 //' @param start_idx Indice inicial (offset R style, ajustado internamente).
 //' @param ext_mode Inteiro: 1=sym, 2=per, 3=zero.
 //' @keywords internal
 // [[Rcpp::export]]
 NumericVector apply_filter_cpp(
     NumericVector x,
     NumericVector coeffs,
     int start_idx,
     int ext_mode
     ) {
   int n = x.size();
   int k = coeffs.size();
   NumericVector y(n);

   // O start_idx vem do R (baseado em 1 ou lógica do Lifting).
   // No R utils.R: base_idx = (1:n) + start_idx.
   // Aqui: i (0..n-1). O índice base correspondente é i + start_idx.

   for (int i = 0; i < n; i++) {
     double sum = 0.0;

     // Loop do filtro (Convolução)
     for (int j = 0; j < k; j++) {
       // Indice virtual no sinal original
       // Em R era: base_idx + (j - 1). Aqui j começa em 0.
       // Ajuste fino necessario para bater com a logica anterior.
       // Se start_idx=0 (R) -> offset 0.
       // Vamos assumir que start_idx passado ja considera a logica do filtro.

       int read_idx = i + start_idx + j;

       double val = get_val_at(x, read_idx, n, ext_mode);
       sum += val * coeffs[j];
     }
     y[i] = sum;
   }

   return y;
 }
