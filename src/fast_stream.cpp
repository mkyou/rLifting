#include "utils.h"
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>

using namespace Rcpp;

// CLASS WAVELET ENGINE (Stateful Processing)
class WaveletEngine {
public:
  // Setup Configuration
  std::vector<LiftingStep> steps;
  double norm_approx;
  double norm_detail;
  int levels;
  int window_size;
  int ext_mode;

  // Buffer State
  std::vector<double> ring_buffer;
  int head;  // Current write index
  int count; // Number of processed samples (for warming)

  // Workspaces (Zero-Allocation)
  // Pre-allocated vectors to avoid heap allocation during the loop
  std::vector<double> work_signal;
  std::vector<std::vector<double>> work_approx;
  std::vector<std::vector<double>> work_detail;

  // Adaptive Threshold Cache
  std::vector<double> current_lambdas;

  // Constructor
  WaveletEngine(
    List r_steps, NumericVector norm,
    int lvl, int w_size, int mode
  ) {
    window_size = w_size;
    levels = lvl;
    ext_mode = mode;
    norm_approx = norm[0];
    norm_detail = norm[1];

    // Parse steps
    int n_steps = r_steps.size();
    for(int i=0; i<n_steps; i++) {
      List s = r_steps[i];
      LiftingStep step;
      step.type = as<std::string>(s["type"]);
      step.coeffs = as<std::vector<double>>(s["coeffs"]);
      step.start_idx = s["start_idx"];
      steps.push_back(step);
    }

    // Memory Allocation (State)
    ring_buffer.resize(window_size, 0.0);
    head = 0;
    count = 0;

    // Workspaces Allocation
    work_signal.resize(window_size);

    work_approx.resize(levels + 1);
    work_detail.resize(levels + 1);
    current_lambdas.resize(levels, 0.0);

    int current_len = window_size;
    for(int j=0; j <= levels; j++) {
      work_approx[j].resize(current_len);
      work_detail[j].resize(current_len);
      current_len = (current_len + 1) / 2;
    }
  }

  // Helper Wrapper for inline utils function
  inline double get_val(const std::vector<double>& x, int i, int n) {
    return get_val_safe(x, i, n, ext_mode);
  }

  // Update Thresholds (Statistical Calculation)
  void update_thresholds(double alpha, double beta) {
    std::vector<double>& d1 = work_detail[0];
    int n1 = d1.size();
    if (n1 == 0) return;

    // Calculate MAD
    std::vector<double> abs_d1(n1);
    for(int i=0; i<n1; i++) abs_d1[i] = std::abs(d1[i]);

    int mid = n1 / 2;
    std::nth_element(abs_d1.begin(), abs_d1.begin() + mid, abs_d1.end());
    double mad = abs_d1[mid];
    double sigma = mad / 0.6745;

    if (sigma < 1e-15) {
      std::fill(current_lambdas.begin(), current_lambdas.end(), 0.0);
      return;
    }

    // Recursive Lambda Calculation
    double lambda_1 = beta * sigma * std::sqrt(2.0 * std::log((double)n1));
    current_lambdas[0] = lambda_1;

    for (int k = 1; k < levels; k++) {
      int lvl_idx = k + 1;
      double prev = current_lambdas[k-1];
      double factor = (double)(lvl_idx - 1) / (double)(lvl_idx + alpha - 1);
      current_lambdas[k] = prev * factor;
    }
  }

  // Core Processing Loop
  double push_and_process(
      double new_val, double alpha, double beta,
      std::string method, int update_freq, int step_iter
  ) {
    // Ingestion into Ring Buffer
    ring_buffer[head] = new_val;
    head = (head + 1) % window_size;
    if (count < window_size) count++;

    // Warming phase: Return raw data until buffer is full
    if (count < window_size) return new_val;

    // Linearization (Ring -> Flat)
    // Copy required to ensure contiguity for the Transform
    for(int i=0; i<window_size; i++) {
      work_signal[i] = ring_buffer[(head + i) % window_size];
    }
    work_approx[0] = work_signal;

    // Multi-Level Forward Decomposition
    for(int j=0; j < levels; j++) {
      const std::vector<double>& input = work_approx[j];
      std::vector<double>& even = work_approx[j+1];
      std::vector<double>& odd  = work_detail[j];

      int n = input.size();
      int n_even = (n + 1) / 2;
      int n_odd = n / 2;

      even.resize(n_even);
      odd.resize(n_odd);

      // Lazy Split
      for(int i=0; i<n_even; i++) even[i] = input[2*i];
      for(int i=0; i<n_odd; i++)  odd[i]  = input[2*i+1];

      // Lifting Steps
      for(const auto& step : steps) {
        int k_filt = step.coeffs.size();
        if (step.type == "predict") {
          for(int i=0; i<n_odd; i++) {
            double sum = 0.0;
            for(int k=0; k<k_filt; k++) sum += get_val(
              even, i + step.start_idx + k, n_even
            ) * step.coeffs[k];
            odd[i] -= sum;
          }
        } else {
          for(int i=0; i<n_even; i++) {
            double sum = 0.0;
            for(int k=0; k<k_filt; k++) sum += get_val(
              odd, i + step.start_idx + k, n_odd
            ) * step.coeffs[k];
            even[i] += sum;
          }
        }
      }

      // Normalization
      for(int i=0; i<n_even; i++) even[i] *= norm_approx;
      for(int i=0; i<n_odd; i++)  odd[i]  *= norm_detail;
    }

    // Thresholding
    // Lazy Update of thresholds
    if (step_iter % update_freq == 0) update_thresholds(alpha, beta);

    for(int j=0; j < levels; j++) {
      double lam = current_lambdas[j];
      double lam_sq = lam * lam;
      std::vector<double>& det = work_detail[j];
      int n_det = det.size();

      for(int i=0; i<n_det; i++) {
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
          // hard: keep unchanged
        }
      }
    }

    // Inverse Reconstruction (Backward)
    for(int j = levels - 1; j >= 0; j--) {
      std::vector<double>& even = work_approx[j+1];
      std::vector<double>& odd  = work_detail[j];

      int n_even = even.size();
      int n_odd = odd.size();

      // De-normalization
      for(int i=0; i<even.size(); i++) even[i] /= norm_approx;
      for(int i=0; i<odd.size(); i++)  odd[i]  /= norm_detail;

      // Reverse Lifting Steps
      for (int k = steps.size() - 1; k >= 0; k--) {
        const auto& step = steps[k];
        int k_filt = step.coeffs.size();
        if (step.type == "predict") {
          for(int i=0; i<n_odd; i++) {
            double sum = 0.0;
            for(int m=0; m<k_filt; m++) sum += get_val(
              even, i + step.start_idx + m, even.size()
            ) * step.coeffs[m];
            odd[i] += sum;
          }
        } else {
          for(int i=0; i<n_even; i++) {
            double sum = 0.0;
            for(int m=0; m<k_filt; m++) sum += get_val(
              odd, i + step.start_idx + m, odd.size()
            ) * step.coeffs[m];
            even[i] -= sum;
          }
        }
      }

      // Merge (Inverse Lazy)
      int target_size = work_approx[j].size();
      for(int i=0; i<even.size(); i++) if (2*i < target_size)
        work_approx[j][2*i] = even[i];
      for(int i=0; i<odd.size(); i++) if (2*i+1 < target_size)
        work_approx[j][2*i+1] = odd[i];
    }

    // Return Causal Result
    // The fully reconstructed signal is in work_approx[0].
    // We return the last valid sample (causal).
    return work_approx[0][window_size - 1];
  }
};



// R INTERFACE (EXPORTS)
// Creates a pointer to the processing engine (Online)
 // @export
 // [[Rcpp::export]]
 SEXP create_engine_cpp(
     List steps, NumericVector norm, int levels,
     int window_size, int ext_mode
 ) {
   // Creates object on Heap
   WaveletEngine* engine = new WaveletEngine(
     steps, norm, levels,
     window_size, ext_mode
   );
   // Returns XPtr (R-managed smart pointer)
   Rcpp::XPtr<WaveletEngine> ptr(engine, true);
   return ptr;
 }

 // Processes a single sample using the persistent engine
 // @export
 // [[Rcpp::export]]
 double process_sample_cpp(
     SEXP engine_ptr, double new_sample,
     double alpha, double beta,
     std::string method, int update_freq,
     int step_iter
 ) {
   Rcpp::XPtr<WaveletEngine> engine(engine_ptr);
   return engine->push_and_process(
       new_sample, alpha, beta, method,
       update_freq, step_iter
   );
 }

 // Batch Simulation (Causal / Turbo)
 //
 // Optimized version to simulate the sequential process on a full vector.
 // @keywords internal
 // [[Rcpp::export]]
 NumericVector run_causal_batch_cpp(
     NumericVector signal, List steps, NumericVector norm,
     int levels, int window_size, double alpha,
     double beta, std::string method, int ext_mode,
     int update_freq
 ) {
   int n = signal.size();
   NumericVector output(n);

   // Instantiate Engine locally
   WaveletEngine engine(steps, norm, levels, window_size, ext_mode);

   for(int i=0; i<n; i++) {
     output[i] = engine.push_and_process(
       signal[i], alpha, beta, method, i, i
     );
   }
   return output;
 }
