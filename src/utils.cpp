#include "utils.h"
using namespace Rcpp;

// Nota: a logica de borda (get_val_safe) foi movida para utils.h
// para ser compartilhada entre Online, Offline e Utils.

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

   // Converte Rcpp::NumericVector para std::vector para usar a funcao
   // inline compartilhada
   // Isso tem um custo minimo de copia,
   // mas garante consistencia matematica absoluta.
   // Como esta funcao nao e usada no loop de
   // Denoising (so em LWT/ILWT isolados),
   // a seguranca vale mais que a micro-performance aqui.
   std::vector<double> x_std = as<std::vector<double>>(x);

   for (int i = 0; i < n; i++) {
     double sum = 0.0;

     for (int j = 0; j < k; j++) {
       // Aritmetica de ponteiro/indice
       int read_idx = i + start_idx + j;

       // Usa a funcao centralizada do utils.h
       sum += get_val_safe(x_std, read_idx, n, ext_mode) * coeffs[j];
     }
     y[i] = sum;
   }

   return y;
 }
