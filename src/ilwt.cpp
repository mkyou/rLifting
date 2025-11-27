#include "utils.h"
#include <string>
#include <vector>

using namespace Rcpp;

//' Core ILWT (C++)
 //'
 //' Executa a reconstrucao multinivel inteiramente em C++.
 //'
 //' @param coeffs_list Lista nomeada com coeficientes (d1..dn, an).
 //' @param steps Lista de passos (do scheme).
 //' @param norm Vetor de normalizacao.
 //' @param levels Numero de niveis.
 //' @param ext_mode Modo de borda.
 //' @param original_len Comprimento original do sinal (para corte final).
 //' @keywords internal
 // [[Rcpp::export]]
 NumericVector ilwt_cpp(List coeffs_list, List steps, NumericVector norm,
                        int levels, int ext_mode, int original_len) {

    // Recupera a aproximacao mais grosseira (an)
    std::string a_name = "a" + std::to_string(levels);
    NumericVector current_app = as<NumericVector>(coeffs_list[a_name]);

    // Loop reverso de niveis (j = levels ate 1)
    for (int j = levels; j >= 1; j--) {

       std::string d_name = "d" + std::to_string(j);
       NumericVector current_det = as<NumericVector>(coeffs_list[d_name]);

       // Clones para manipular
       // (Rcpp pointers sao perigosos se nao clonar aqui)
       NumericVector even = clone(current_app);
       NumericVector odd  = clone(current_det);

       // Desfazer normalizacao
       // norm[0]=K, norm[1]=1/K.
       // Na ida multiplicamos. Na volta dividimos.
       for(int i=0; i<even.size(); i++) even[i] /= norm[0];
       for(int i=0; i<odd.size(); i++)  odd[i]  /= norm[1];

       // Reverse lifting steps
       // Iterar de tras pra frente
       int n_steps = steps.size();
       for (int k = n_steps - 1; k >= 0; k--) {
          List step = steps[k];
          std::string type = as<std::string>(step["type"]);
          NumericVector coeffs = step["coeffs"];
          int start_idx = step["start_idx"];

          // Na ida: Predict faz d = d - P(e)
          // Na volta: d = d + P(e) (SOMA)
          if (type == "predict") {
             NumericVector pred = apply_filter_cpp(
                 even, coeffs,
                 start_idx, ext_mode
              );
             int len = std::min(odd.size(), pred.size());
             for(int m=0; m<len; m++) odd[m] += pred[m];
          }
          // Na ida: update faz a = a + U(d)
          // Na volta: a = a - U(d) (SUBTRACAO)
          else if (type == "update") {
             NumericVector upd = apply_filter_cpp(
                 odd, coeffs,
                 start_idx, ext_mode
              );
             int len = std::min(even.size(), upd.size());
             for(int m=0; m<len; m++) even[m] -= upd[m];
          }
       }

       // Merge (Inverse Lazy Wavelet)
       // Intercalar: even[0], odd[0], even[1], odd[1]...
       int n_total = even.size() + odd.size();
       NumericVector merged(n_total);

       for(int i=0; i<even.size(); i++) {
          if (2*i < n_total) merged[2*i] = even[i];
       }
       for(int i=0; i<odd.size(); i++) {
          if (2*i + 1 < n_total) merged[2*i + 1] = odd[i];
       }

       current_app = merged;
    }

    // Corte final para garantir tamanho exato do original
    if (current_app.size() > original_len) {
       // Cria vetor novo do tamanho certo
       NumericVector final_res(original_len);
       for(int i=0; i<original_len; i++) final_res[i] = current_app[i];
       return final_res;
    }

    return current_app;
 }
