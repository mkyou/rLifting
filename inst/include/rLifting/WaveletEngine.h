#ifndef RLIFTING_WAVELET_ENGINE_H
#define RLIFTING_WAVELET_ENGINE_H

#include "utils.h"
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

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
  int ll_k;
  bool irregular;
  double scad_a;
  bool lambdas_initialized;

  // Buffer State
  std::vector<double> ring_buffer;
  std::vector<double> ring_buffer_t;  // parallel t positions (irregular mode)
  int head;
  int count;

  // Workspaces (Zero-Allocation)
  std::vector<double> work_signal;
  std::vector<double> work_t;                      // linearized t window
  std::vector<std::vector<double>> work_approx;
  std::vector<std::vector<double>> work_detail;
  std::vector<std::vector<double>> work_t_approx;  // t at each level (t_even)
  std::vector<std::vector<double>> work_t_detail;  // t at each level (t_odd)

  // Adaptive Threshold Cache
  std::vector<double> current_lambdas;

  // Constructor
  WaveletEngine(List r_steps, NumericVector norm, int lvl, int w_size,
                int mode, bool irreg = false, int k = 2,
                double scad_a_val = 3.7) {
    window_size = w_size;
    levels = lvl;
    ext_mode = mode;
    ll_k = k;
    norm_approx = norm[0];
    norm_detail = norm[1];
    irregular = irreg;
    scad_a = scad_a_val;
    lambdas_initialized = false;

    int n_steps = r_steps.size();
    for (int i = 0; i < n_steps; i++) {
      List s = r_steps[i];
      LiftingStep step;
      step.type = as<std::string>(s["type"]);
      step.coeffs = as<std::vector<double>>(s["coeffs"]);
      step.start_idx = s["start_idx"];
      step.degree = s.containsElementNamed("degree") ? (int)s["degree"] : -1;
      steps.push_back(step);
    }

    ring_buffer.resize(window_size, 0.0);
    head = 0;
    count = 0;

    work_signal.resize(window_size);
    work_approx.resize(levels + 1);
    work_detail.resize(levels + 1);
    current_lambdas.resize(levels, 0.0);

    int current_len = window_size;
    for (int j = 0; j <= levels; j++) {
      work_approx[j].resize(current_len);
      work_detail[j].resize(current_len);
      current_len = (current_len + 1) / 2;
    }

    if (irregular) {
      ring_buffer_t.resize(window_size, 0.0);
      work_t.resize(window_size, 0.0);
      work_t_approx.resize(levels + 1);
      work_t_detail.resize(levels);
      current_len = window_size;
      for (int j = 0; j <= levels; j++) {
        work_t_approx[j].resize(current_len, 0.0);
        if (j < levels) work_t_detail[j].resize(current_len / 2, 0.0);
        current_len = (current_len + 1) / 2;
      }
    }
  }

  inline double get_val(const std::vector<double> &x, int i, int n) {
    return get_val_safe(x, i, n, ext_mode, ll_k);
  }

  void update_thresholds(double alpha, double beta,
                         const std::string &threshold_method) {
    if (threshold_method == "sure") {
      for (int j = 0; j < levels; j++) {
        current_lambdas[j] = compute_sure_lambda_level(work_detail[j]);
      }
      return;
    }

    std::vector<double> &d1 = work_detail[0];
    int n1 = d1.size();
    if (n1 == 0) return;

    double sigma = compute_mad(d1) / 0.6745;

    if (sigma < 1e-15) {
      std::fill(current_lambdas.begin(), current_lambdas.end(), 0.0);
      return;
    }

    double lambda_1 = beta * sigma * std::sqrt(2.0 * std::log((double)n1));
    current_lambdas[0] = lambda_1;

    for (int k = 1; k < levels; k++) {
      int lvl_idx = k + 1;
      double prev = current_lambdas[k - 1];
      double factor = (double)(lvl_idx - 1) / (double)(lvl_idx + alpha - 1);
      current_lambdas[k] = prev * factor;
    }
  }

  // Core Processing Loop
  double push_and_process(double new_val, double t_val, double alpha,
                          double beta, std::string method, int update_freq,
                          int step_iter,
                          const std::string &threshold_method = "universal") {
    ring_buffer[head] = new_val;
    if (irregular) ring_buffer_t[head] = t_val;
    head = (head + 1) % window_size;
    if (count < window_size) count++;

    if (count < window_size) return new_val;

    bool use_os = (ext_mode == 5);

    // Linearize ring buffer
    for (int i = 0; i < window_size; i++)
      work_signal[i] = ring_buffer[(head + i) % window_size];
    work_approx[0] = work_signal;

    if (irregular) {
      for (int i = 0; i < window_size; i++)
        work_t[i] = ring_buffer_t[(head + i) % window_size];
      work_t_approx[0] = work_t;
    }

    // Multi-Level Forward Decomposition
    for (int j = 0; j < levels; j++) {
      const std::vector<double> &input = work_approx[j];
      std::vector<double> &even = work_approx[j + 1];
      std::vector<double> &odd  = work_detail[j];

      int n = input.size();
      int n_even = (n + 1) / 2;
      int n_odd  = n / 2;

      even.resize(n_even);
      odd.resize(n_odd);

      for (int i = 0; i < n_even; i++) even[i] = input[2 * i];
      for (int i = 0; i < n_odd;  i++) odd[i]  = input[2 * i + 1];

      if (irregular) {
        const std::vector<double> &t_cur = work_t_approx[j];
        work_t_approx[j + 1].resize(n_even);
        work_t_detail[j].resize(n_odd);
        for (int i = 0; i < n_even; i++) work_t_approx[j + 1][i] = t_cur[2 * i];
        for (int i = 0; i < n_odd;  i++) work_t_detail[j][i]     = t_cur[2 * i + 1];
      }

      for (const auto &step : steps) {
        int k_filt = (int)step.coeffs.size();
        const double* c = step.coeffs.data();

        if (step.type == "predict") {
          if (use_os) {
            for (int i = 0; i < n_odd; i++)
              odd[i] -= onesided_conv(even, n_even, c, k_filt, step.start_idx, i);
          } else if (irregular && step.degree >= 0) {
            const std::vector<double> &t_even = work_t_approx[j + 1];
            const std::vector<double> &t_odd  = work_t_detail[j];
            std::vector<double> x_nbr(k_filt), t_nbr(k_filt);
            for (int i = 0; i < n_odd; i++) {
              for (int m = 0; m < k_filt; m++) {
                int idx = i + step.start_idx + m;
                x_nbr[m] = get_val_safe(even, idx, n_even, ext_mode, ll_k);
                t_nbr[m] = get_t_extrap(t_even, idx, n_even);
              }
              odd[i] -= interp_predict(x_nbr, t_nbr, k_filt, t_odd[i]);
            }
          } else {
            for (int i = 0; i < n_odd; i++) {
              double sum = 0.0;
              for (int k = 0; k < k_filt; k++)
                sum += get_val(even, i + step.start_idx + k, n_even) * c[k];
              odd[i] -= sum;
            }
          }
        } else {
          if (use_os) {
            for (int i = 0; i < n_even; i++)
              even[i] += onesided_conv(odd, n_odd, c, k_filt, step.start_idx, i);
          } else {
            for (int i = 0; i < n_even; i++) {
              double sum = 0.0;
              for (int k = 0; k < k_filt; k++)
                sum += get_val(odd, i + step.start_idx + k, n_odd) * c[k];
              even[i] += sum;
            }
          }
        }
      }

      for (int i = 0; i < n_even; i++) even[i] *= norm_approx;
      for (int i = 0; i < n_odd;  i++) odd[i]  *= norm_detail;
    }

    // Thresholding. update_freq <= 0 means "freeze after first update".
    bool should_update = (update_freq > 0)
                            ? (step_iter % update_freq == 0)
                            : !lambdas_initialized;
    if (should_update) {
      update_thresholds(alpha, beta, threshold_method);
      lambdas_initialized = true;
    }

    // SCAD shape parameter (Fan-Li 2001), configurable per engine.
    const double SCAD_A = scad_a;

    for (int j = 0; j < levels; j++) {
      double lam    = current_lambdas[j];
      double lam_sq = lam * lam;
      double two_lam = 2.0 * lam;
      double a_lam = SCAD_A * lam;
      double scad_denom = SCAD_A - 2.0;
      std::vector<double> &det = work_detail[j];
      int n_det = det.size();

      for (int i = 0; i < n_det; i++) {
        double val     = det[i];
        double abs_val = std::abs(val);
        if (abs_val < lam) {
          det[i] = 0.0;
        } else {
          if (method == "soft") {
            det[i] = (val > 0) ? (abs_val - lam) : -(abs_val - lam);
          } else if (method == "semisoft") {
            double s = std::sqrt(val * val - lam_sq);
            det[i] = (val > 0) ? s : -s;
          } else if (method == "scad") {
            double sgn = (val > 0) ? 1.0 : -1.0;
            if (abs_val <= two_lam) {
              det[i] = sgn * (abs_val - lam);
            } else if (abs_val <= a_lam) {
              det[i] = ((SCAD_A - 1.0) * val - sgn * a_lam) / scad_denom;
            }
            // else (|val| > a*lam): identity, det[i] unchanged
          }
        }
      }
    }

    // Inverse Reconstruction (Backward)
    for (int j = levels - 1; j >= 0; j--) {
      std::vector<double> &even = work_approx[j + 1];
      std::vector<double> &odd  = work_detail[j];

      int n_even = even.size();
      int n_odd  = odd.size();

      for (int i = 0; i < n_even; i++) even[i] /= norm_approx;
      for (int i = 0; i < n_odd;  i++) odd[i]  /= norm_detail;

      // t positions for this level are stored from the forward pass
      for (int k = (int)steps.size() - 1; k >= 0; k--) {
        const auto &step = steps[k];
        int k_filt = (int)step.coeffs.size();
        const double* c = step.coeffs.data();
        int sz_even = (int)even.size();
        int sz_odd  = (int)odd.size();

        if (step.type == "predict") {
          if (use_os) {
            for (int i = 0; i < n_odd; i++)
              odd[i] += onesided_conv(even, sz_even, c, k_filt, step.start_idx, i);
          } else if (irregular && step.degree >= 0) {
            const std::vector<double> &t_even = work_t_approx[j + 1];
            const std::vector<double> &t_odd  = work_t_detail[j];
            std::vector<double> x_nbr(k_filt), t_nbr(k_filt);
            for (int i = 0; i < sz_odd; i++) {
              for (int m = 0; m < k_filt; m++) {
                int idx = i + step.start_idx + m;
                x_nbr[m] = get_val_safe(even, idx, sz_even, ext_mode, ll_k);
                t_nbr[m] = get_t_extrap(t_even, idx, sz_even);
              }
              odd[i] += interp_predict(x_nbr, t_nbr, k_filt, t_odd[i]);
            }
          } else {
            for (int i = 0; i < sz_odd; i++) {
              double sum = 0.0;
              for (int m = 0; m < k_filt; m++)
                sum += get_val(even, i + step.start_idx + m, sz_even) * c[m];
              odd[i] += sum;
            }
          }
        } else {
          if (use_os) {
            for (int i = 0; i < n_even; i++)
              even[i] -= onesided_conv(odd, sz_odd, c, k_filt, step.start_idx, i);
          } else {
            for (int i = 0; i < n_even; i++) {
              double sum = 0.0;
              for (int m = 0; m < k_filt; m++)
                sum += get_val(odd, i + step.start_idx + m, sz_odd) * c[m];
              even[i] -= sum;
            }
          }
        }
      }

      int target_size = (int)work_approx[j].size();
      for (int i = 0; i < (int)even.size(); i++)
        if (2 * i < target_size) work_approx[j][2 * i] = even[i];
      for (int i = 0; i < (int)odd.size(); i++)
        if (2 * i + 1 < target_size) work_approx[j][2 * i + 1] = odd[i];
    }

    return work_approx[0][window_size - 1];
  }
};

#endif
