#include "../inst/include/rLifting/WaveletEngine.h"
#include "utils.h"
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

using namespace Rcpp;

// R INTERFACE (EXPORTS)
// Creates a pointer to the processing engine (Online)
// @export
// [[Rcpp::export]]
SEXP create_engine_cpp(List steps, NumericVector norm, int levels,
                       int window_size, int ext_mode, bool irregular,
                       int ll_k = 2) {
  WaveletEngine *engine =
      new WaveletEngine(steps, norm, levels, window_size, ext_mode, irregular, ll_k);
  Rcpp::XPtr<WaveletEngine> ptr(engine, true);
  return ptr;
}

// Processes a single sample using the persistent engine
// @export
// [[Rcpp::export]]
double process_sample_cpp(SEXP engine_ptr, double new_sample, double t_val,
                          double alpha, double beta, std::string method,
                          int update_freq, int step_iter,
                          std::string threshold_method = "universal") {
  Rcpp::XPtr<WaveletEngine> engine(engine_ptr);
  return engine->push_and_process(new_sample, t_val, alpha, beta, method,
                                  update_freq, step_iter, threshold_method);
}

// Batch Simulation (Causal / Turbo)
//
// Optimized version to simulate the sequential process on a full vector.
// @keywords internal
// [[Rcpp::export]]
NumericVector run_causal_batch_cpp(NumericVector signal, List steps,
                                   NumericVector norm, int levels,
                                   int window_size, double alpha, double beta,
                                   std::string method, int ext_mode,
                                   int update_freq, NumericVector t,
                                   int ll_k = 2,
                                   std::string threshold_method = "universal") {
  int n = signal.size();
  NumericVector output(n);

  bool irreg = ((int)t.size() == n);
  WaveletEngine engine(steps, norm, levels, window_size, ext_mode, irreg, ll_k);

  for (int i = 0; i < n; i++) {
    double t_val = irreg ? t[i] : 0.0;
    output[i] = engine.push_and_process(signal[i], t_val, alpha, beta, method,
                                        update_freq, i, threshold_method);
  }
  return output;
}
