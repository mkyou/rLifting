#include "utils.h"
#include <string>
#include <vector>

using namespace Rcpp;

// Helper para aplicar threshold no objeto lwt
void apply_threshold_in_place(
    List coeffs, NumericVector lambdas,
    std::string method
) {
  int max_level = lambdas.size();

  for(int j=0; j<max_level; j++) {
    std::string level_name = "d" + std::to_string(j+1);
    if (!coeffs.containsElementNamed(level_name.c_str())) continue;

    NumericVector details = coeffs[level_name];
    double lam = lambdas[j];
    NumericVector filtered;

    if (method == "hard") {
      filtered = threshold_hard_cpp(details, lam);
    } else if (method == "soft") {
      filtered = threshold_soft_cpp(details, lam);
    } else {
      filtered = threshold_semisoft_cpp(details, lam);
    }
    coeffs[level_name] = filtered;
  }
}

//' Causal Denoising Batch (C++ Optimized Ring Buffer)
 //'
 //' Otimizacoes:
 //' 1. Ring Buffer: Evita deslocamento de memoria (O(1) insercao).
 //' 2. Lazy Threshold: Recalcula estatisticas de ruido apenas a cada k passos.
 //'
 //' @param signal Sinal completo.
 //' @param steps Passos da ondaleta.
 //' @param norm Normalizacao.
 //' @param levels Niveis de decomposicao.
 //' @param window_size Tamanho da janela.
 //' @param alpha Parametro threshold.
 //' @param beta Parametro threshold.
 //' @param method Metodo de threshold.
 //' @param ext_mode Modo de extensao.
 //' @param update_freq Frequencia de atualizacao do threshold
 //' (ex: 1 = sempre, 10 = a cada 10).
 //' @keywords internal
 // [[Rcpp::export]]
 NumericVector run_causal_batch_cpp(
     NumericVector signal,
     List steps,
     NumericVector norm,
     int levels,
     int window_size,
     double alpha,
     double beta,
     std::string method,
     int ext_mode,
     int update_freq = 1 // Padrao: atualiza sempre (comportamento antigo)
 ) {

   int n = signal.size();
   NumericVector output(n);

   // --- OTIMIZACAO 1: RING BUFFER ---
   // Alocamos tamanho fixo e usamos ponteiro 'head'
   std::vector<double> ring_buffer(window_size, 0.0);
   int head = 0; // Aponta para onde sera escrita a proxima amostra
   int count = 0; // Quantos elementos validos temos

   // Cache para o Threshold (Lazy Update)
   NumericVector current_lambdas;

   for(int i = 0; i < n; i++) {

     // 1. Ingestao Circular
     ring_buffer[head] = signal[i];
     head = (head + 1) % window_size; // Avança circularmente
     if (count < window_size) count++;

     // 2. Warming Phase (Bypass)
     if (count < window_size) {
       output[i] = signal[i];
       continue;
     }

     // 3. Linearizar o Buffer (Custo necessario para LWT)
     // Precisamos passar um vetor contiguo para a transformada
     // Copia ordenada a partir do 'tail'
     NumericVector curr_sig(window_size);
     int tail = head;
     // O head atual e o mais antigo
     // (pois acabamos de sobrescrever o mais antigo anterior, nao, espera)
     // Correcao logica Ring Buffer:
     // Se head esta em 5, o ultimo escrito foi 4. O mais antigo e 5.
     // A ordem cronologica e: buffer[head],
     // buffer[head+1]... buffer[window-1], buffer[0]... buffer[head-1]

     for(int k=0; k<window_size; k++) {
       curr_sig[k] = ring_buffer[(head + k) % window_size];
     }

     // 4. Forward LWT
     List lwt_obj = lwt_cpp(curr_sig, steps, norm, levels, ext_mode);

     // 5. Thresholding (Otimizacao Lazy)
     bool need_update = (i % update_freq == 0) ||
       (current_lambdas.size() == 0);

     if (lwt_obj.containsElementNamed("d1")) {
       if (need_update) {
         // Recalcula estatistica (caro)
         NumericVector d1 = lwt_obj["d1"];
         int max_l = 0;
         while(
           lwt_obj.containsElementNamed(
             ("d" + std::to_string(max_l+1)).c_str()
           )
          ) max_l++;
         current_lambdas = compute_thresholds_cpp(d1, max_l, alpha, beta);
       }
       // Aplica (barato)
       apply_threshold_in_place(lwt_obj, current_lambdas, method);
     } else {
       output[i] = signal[i]; continue;
     }

     // 6. Inverse ILWT
     NumericVector rec = ilwt_cpp(
       lwt_obj, steps, norm,
       levels, ext_mode, window_size
     );
     output[i] = rec[rec.size() - 1];
   }

   return output;
 }
