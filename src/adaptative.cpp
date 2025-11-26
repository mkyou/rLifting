#include <Rcpp.h>
#include <algorithm> // std::nth_element (para mediana rapida)
#include <cmath> // log, sqrt

using namespace Rcpp;

//' Helper: Mediana Absoluta (MAD) em C++
 //' @keywords internal
 double calc_mad(NumericVector x) {
   int n = x.size();
   if (n == 0) return 0.0;

   // Copia valores absolutos
   std::vector<double> abs_x(n);
   for(int i=0; i<n; i++) abs_x[i] = std::abs(x[i]);

   // Algoritmo de Selecao (mais rapido que sort total)
   // Coloca a mediana na posicao n/2
   int mid = n / 2;
   std::nth_element(abs_x.begin(), abs_x.begin() + mid, abs_x.end());

   double median = abs_x[mid];

   // Se n for par, a definicao do R e a media dos dois centrais.
   // Para velocidade em real-time, usar o elemento central e aceitavel,
   // mas para precisao exata, vamos fazer a media se par.
   if (n % 2 == 0) {
     std::nth_element(abs_x.begin(), abs_x.begin() + mid - 1, abs_x.end());
     median = (median + abs_x[mid - 1]) / 2.0;
   }

   return median;
 }

 //' Adaptive Threshold Calculation (C++)
 //'
 //' Calcula os limiares recursivos diretamente em C++.
 //'
 //' @param d1 Coeficientes de detalhe do nivel 1.
 //' @param max_level Numero maximo de niveis.
 //' @param alpha Parametro de decaimento.
 //' @param beta Parametro de escala.
 //' @keywords internal
 // [[Rcpp::export]]
 NumericVector compute_thresholds_cpp(
     NumericVector d1, int max_level,
     double alpha, double beta
 ) {
   NumericVector lambdas(max_level);

   // 1. Estimar Ruido (Sigma) via MAD
   double mad_val = calc_mad(d1);
   double sigma = mad_val / 0.6745;

   if (sigma < 1e-15) {
     return lambdas; // Retorna tudo zero
   }

   // 2. Lambda Nivel 1
   int n1 = d1.size();
   // Log natural em C++ e std::log
   double lambda_1 = beta * sigma * std::sqrt(2.0 * std::log((double)n1));
   lambdas[0] = lambda_1;

   // 3. Recursao para niveis superiores (i = 2 ate max)
   // Nota: No vetor C++ indices sao 0..max-1.
   // A formula usa 'i' como nivel (1, 2...).

   for (int k = 1; k < max_level; k++) {
     int level = k + 1; // Nivel atual (2, 3...)
     double prev = lambdas[k-1];

     // Formula Eq(9): L_i = L_{i-1} * (i-1) / (i + alpha - 1)
     double factor = (double)(level - 1) / (double)(level + alpha - 1);
     lambdas[k] = prev * factor;
   }

   return lambdas;
 }
