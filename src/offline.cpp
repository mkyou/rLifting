#include "utils.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>

using namespace Rcpp;

// Helper para calcular MAD e Thresholds
std::vector<double> compute_thresholds_internal(
    const std::vector<double>& d1, int max_level, double alpha, double beta
) {
  std::vector<double> lambdas(max_level, 0.0);
  int n = d1.size();
  if (n == 0) return lambdas;

  // MAD
  std::vector<double> abs_x(n);
  for(int i=0; i<n; i++) abs_x[i] = std::abs(d1[i]);
  int mid = n / 2;
  std::nth_element(abs_x.begin(), abs_x.begin() + mid, abs_x.end());
  double mad = abs_x[mid];
  double sigma = mad / 0.6745;

  if (sigma < 1e-15) return lambdas;

  double lambda_1 = beta * sigma * std::sqrt(2.0 * std::log((double)n));
  lambdas[0] = lambda_1;

  for (int k = 1; k < max_level; k++) {
    int level = k + 1;
    double prev = lambdas[k-1];
    double factor = (double)(level - 1) / (double)(level + alpha - 1);
    lambdas[k] = prev * factor;
  }
  return lambdas;
}

//' Unified offline denoising
 //'
 //' Executa LWT -> threshold -> ILWT inteiramente em C++ sem overhead de R.
 //' @keywords internal
 // [[Rcpp::export]]
 NumericVector denoise_offline_cpp(
     NumericVector signal,
     List steps,
     NumericVector norm,
     int levels,
     double alpha,
     double beta,
     std::string method,
     int ext_mode
 ) {
   // Setup & Parsing (uma unica vez)
   std::vector<LiftingStep> cpp_steps;
   int n_steps = steps.size();
   for(int i=0; i<n_steps; i++) {
     List s = steps[i];
     LiftingStep step;
     step.type = as<std::string>(s["type"]);
     step.coeffs = as<std::vector<double>>(s["coeffs"]);
     step.start_idx = s["start_idx"];
     cpp_steps.push_back(step);
   }

   double norm_approx = norm[0];
   double norm_detail = norm[1];

   // Estruturas de dados
   // current_app comeca com o sinal original
   std::vector<double> current_app = as<std::vector<double>>(signal);
   // Armazena detalhes de cada nivel (d1 em details[0])
   std::vector<std::vector<double>> details(levels);

   // FORWARD LWT
   for (int j = 0; j < levels; j++) {
     int n = current_app.size();
     int n_even = (n + 1) / 2;
     int n_odd = n / 2;

     std::vector<double> even(n_even);
     std::vector<double> odd(n_odd);

     // Split
     for (int i = 0; i < n_even; i++) even[i] = current_app[2*i];
     for (int i = 0; i < n_odd; i++)  odd[i]  = current_app[2*i+1];

     // Lifting
     for (const auto& step : cpp_steps) {
       int k_filt = step.coeffs.size();
       if (step.type == "predict") {
         for(int i=0; i<n_odd; i++) {
           double sum = 0.0;
           for(int k=0; k<k_filt; k++) sum += get_val_safe(
             even, i + step.start_idx + k, n_even, ext_mode
           ) * step.coeffs[k];
           odd[i] -= sum;
         }
       } else {
         for(int i=0; i<n_even; i++) {
           double sum = 0.0;
           for(int k=0; k<k_filt; k++) sum += get_val_safe(
             odd, i + step.start_idx + k, n_odd, ext_mode
           ) * step.coeffs[k];
           even[i] += sum;
         }
       }
     }

     // Norm
     for(int i=0; i<n_even; i++) even[i] *= norm_approx;
     for(int i=0; i<n_odd; i++)  odd[i]  *= norm_detail;

     // Store
     details[j] = odd;
     current_app = even; // Proximo nivel opera na approx
   }

   // THRESHOLDING
   std::vector<double> lambdas = compute_thresholds_internal(
     details[0], levels,
     alpha, beta
   );

   for (int j = 0; j < levels; j++) {
     double lam = lambdas[j];
     double lam_sq = lam * lam;
     std::vector<double>& det = details[j];

     for (size_t i = 0; i < det.size(); i++) {
       double val = det[i];
       double abs_val = std::abs(val);

       if (abs_val < lam) {
         det[i] = 0.0;
       } else {
         if (method == "soft") {
           det[i] = (val > 0) ? (abs_val - lam) : -(abs_val - lam);
         } else if (method == "semisoft") {
           double s = std::sqrt(val*val - lam_sq);
           det[i] = (val > 0) ? s : -s;
         }
       }
     }
   }

   // INVERSE ILWT
   int original_len = signal.size();

   for (int j = levels - 1; j >= 0; j--) {
     std::vector<double>& even = current_app; // Apx do nivel anterior (j+1)
     std::vector<double>& odd = details[j];   // Detalhe deste nivel (j)

     // Un-Norm
     for(size_t i=0; i<even.size(); i++) even[i] /= norm_approx;
     for(size_t i=0; i<odd.size(); i++)  odd[i]  /= norm_detail;

     // Reverse Lifting
     for (int k = cpp_steps.size() - 1; k >= 0; k--) {
       const auto& step = cpp_steps[k];
       int k_filt = step.coeffs.size();
       int n_even = even.size();
       int n_odd = odd.size();

       if (step.type == "predict") { // Volta: Soma
         for(int i=0; i<n_odd; i++) {
           double sum = 0.0;
           for(int m=0; m<k_filt; m++) sum += get_val_safe(
             even, i + step.start_idx + m, n_even, ext_mode
           ) * step.coeffs[m];
           odd[i] += sum;
         }
       } else { // Volta: Subtrai
         for(int i=0; i<n_even; i++) {
           double sum = 0.0;
           for(int m=0; m<k_filt; m++) sum += get_val_safe(
             odd, i + step.start_idx + m, n_odd, ext_mode
           ) * step.coeffs[m];
           even[i] -= sum;
         }
       }
     }

     // Merge
     std::vector<double> merged(even.size() + odd.size());
     for(size_t i=0; i<even.size(); i++) merged[2*i] = even[i];
     for(size_t i=0; i<odd.size(); i++)  merged[2*i+1] = odd[i];

     current_app = merged;
   }

   // Ajuste de tamanho final (corte se necessario)
   if ((int)current_app.size() > original_len) {
     current_app.resize(original_len);
   }

   return wrap(current_app);
 }
