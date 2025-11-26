#include "utils.h"
#include <string>
#include <vector>

using namespace Rcpp;

//' Core LWT (C++)
 //'
 //' Executa a decomposicao multinivel inteiramente em C++.
 //'
 //' @param signal Sinal de entrada.
 //' @param steps Lista de passos (extraida do scheme).
 //' @param norm Vetor de normalizacao (c(K, 1/K)).
 //' @param levels Numero de niveis.
 //' @param ext_mode Modo de borda (inteiro).
 //' @keywords internal
 // [[Rcpp::export]]
 List lwt_cpp(
     NumericVector signal, List steps,
     NumericVector norm, int levels,
     int ext_mode
     ) {

   List coeffs_out;
   NumericVector current_app = clone(signal);
   // Copia para nao alterar o original

   // Loop de Niveis (j = 1 ate levels)
   for (int j = 1; j <= levels; j++) {
     int n = current_app.size();

     // 1. Lazy Wavelet (Split)
     // R (1-based): Impares=x[1,3..] (Pares orig),
     // Pares=x[2,4..] (Impares orig)
     // C++ (0-based):
     //   Even (R concept) -> indices 0, 2, 4...
     //   Odd  (R concept) -> indices 1, 3, 5...

     // Calcula tamanhos
     int n_even = (n + 1) / 2; // ceil(n/2)
     int n_odd  = n / 2;       // floor(n/2)

     NumericVector even(n_even);
     NumericVector odd(n_odd);

     for (int i = 0; i < n_even; i++) even[i] = current_app[2*i];
     for (int i = 0; i < n_odd; i++)  odd[i]  = current_app[2*i + 1];

     // 2. Lifting Steps
     int n_steps = steps.size();

     for (int k = 0; k < n_steps; k++) {
       List step = steps[k];
       std::string type = as<std::string>(step["type"]);
       NumericVector coeffs = step["coeffs"];
       int start_idx = step["start_idx"];

       if (type == "predict") {
         // d = odd - P(even)
         NumericVector pred = apply_filter_cpp(
           even, coeffs,
           start_idx, ext_mode
         );

         // Subtracao (com cuidado se tamanhos diferirem em 1 devido ao split)
         int len = std::min(odd.size(), pred.size());
         for(int m=0; m<len; m++) odd[m] -= pred[m];

       } else if (type == "update") {
         // a = even + U(d)
         NumericVector upd = apply_filter_cpp(
           odd, coeffs,
           start_idx, ext_mode
         );

         // Soma
         int len = std::min(even.size(), upd.size());
         for(int m=0; m<len; m++) even[m] += upd[m];
       }
     }

     // 3. Normalizacao
     // norm[0] = K (approx), norm[1] = 1/K (detail)
     for(int m=0; m<even.size(); m++) even[m] *= norm[0];
     for(int m=0; m<odd.size(); m++)  odd[m]  *= norm[1];

     // 4. Armazenar e Preparar proximo nivel
     // Nomes no R sao "d1", "d2"...
     std::string d_name = "d" + std::to_string(j);
     coeffs_out[d_name] = odd;

     current_app = even;
   }

   // Ultima aproximacao "an"
   std::string a_name = "a" + std::to_string(levels);
   coeffs_out[a_name] = current_app;

   return coeffs_out;
 }
