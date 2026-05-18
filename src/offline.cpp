#include "utils.h"
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

using namespace Rcpp;

// SureShrink: per-level lambdas computed independently via SURE.
// (compute_sure_lambda_level is now inline in utils.h for reuse by
// WaveletEngine.)
std::vector<double> compute_thresholds_sure_internal(
    const std::vector<std::vector<double>> &details, int max_level) {
  std::vector<double> lambdas(max_level, 0.0);
  for (int j = 0; j < max_level && j < (int)details.size(); j++) {
    lambdas[j] = compute_sure_lambda_level(details[j]);
  }
  return lambdas;
}

// Helper to compute thresholds internally
std::vector<double> compute_thresholds_internal(const std::vector<double> &d1,
                                                int max_level, double alpha,
                                                double beta) {
  std::vector<double> lambdas(max_level, 0.0);
  int n = d1.size();
  if (n == 0)
    return lambdas;

  std::vector<double> abs_x(n);
  for (int i = 0; i < n; i++)
    abs_x[i] = std::abs(d1[i]);
  int mid = n / 2;
  std::nth_element(abs_x.begin(), abs_x.begin() + mid, abs_x.end());
  double mad = abs_x[mid];
  double sigma = mad / 0.6745;

  if (sigma < 1e-15)
    return lambdas;

  double lambda_1 = beta * sigma * std::sqrt(2.0 * std::log((double)n));
  lambdas[0] = lambda_1;

  for (int k = 1; k < max_level; k++) {
    int level = k + 1;
    double prev = lambdas[k - 1];
    double factor = (double)(level - 1) / (double)(level + alpha - 1);
    lambdas[k] = prev * factor;
  }
  return lambdas;
}

// Unified Offline Denoising
//
// Executes LWT -> threshold -> ILWT entirely in C++ without R overhead.
// When t is non-empty (length == signal length), position-aware Lagrange
// interpolation is applied in predict steps with degree >= 0.
// @keywords internal
// [[Rcpp::export]]
NumericVector denoise_offline_cpp(NumericVector signal, List steps,
                                  NumericVector norm, int levels, double alpha,
                                  double beta, std::string method,
                                  int ext_mode, NumericVector t, int ll_k = 2,
                                  std::string threshold_method = "universal") {
  // Setup & Parsing
  std::vector<LiftingStep> cpp_steps;
  int n_steps = steps.size();
  for (int i = 0; i < n_steps; i++) {
    List s = steps[i];
    LiftingStep step;
    step.type = as<std::string>(s["type"]);
    step.coeffs = as<std::vector<double>>(s["coeffs"]);
    step.start_idx = s["start_idx"];
    step.degree = s.containsElementNamed("degree") ? (int)s["degree"] : -1;
    cpp_steps.push_back(step);
  }

  double norm_approx = norm[0];
  double norm_detail = norm[1];

  bool irregular = ((int)t.size() == signal.size());
  bool use_os = (ext_mode == 5);

  std::vector<double> current_app = as<std::vector<double>>(signal);
  std::vector<std::vector<double>> details(levels);

  // Per-level t positions: t_levels[j] holds positions before the j-th split.
  std::vector<std::vector<double>> t_levels(levels + 1);
  if (irregular) t_levels[0] = as<std::vector<double>>(t);

  // FORWARD LWT
  for (int j = 0; j < levels; j++) {
    int n = current_app.size();
    int n_even = (n + 1) / 2;
    int n_odd = n / 2;

    std::vector<double> even(n_even);
    std::vector<double> odd(n_odd);

    for (int i = 0; i < n_even; i++) even[i] = current_app[2 * i];
    for (int i = 0; i < n_odd; i++)  odd[i]  = current_app[2 * i + 1];

    std::vector<double> t_even, t_odd;
    if (irregular) {
      const std::vector<double>& t_cur = t_levels[j];
      t_even.resize(n_even); t_odd.resize(n_odd);
      for (int i = 0; i < n_even; i++) t_even[i] = t_cur[2 * i];
      for (int i = 0; i < n_odd;  i++) t_odd[i]  = t_cur[2 * i + 1];
      t_levels[j + 1] = t_even;
    }

    for (const auto &step : cpp_steps) {
      int k_filt = (int)step.coeffs.size();
      const double* c = step.coeffs.data();
      if (step.type == "predict") {
        if (use_os) {
          for (int i = 0; i < n_odd; i++)
            odd[i] -= onesided_conv(even, n_even, c, k_filt, step.start_idx, i);
        } else if (irregular && step.degree >= 0) {
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
              sum += get_val_safe(even, i + step.start_idx + k, n_even, ext_mode, ll_k) * c[k];
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
              sum += get_val_safe(odd, i + step.start_idx + k, n_odd, ext_mode, ll_k) * c[k];
            even[i] += sum;
          }
        }
      }
    }

    for (int i = 0; i < n_even; i++) even[i] *= norm_approx;
    for (int i = 0; i < n_odd; i++)  odd[i]  *= norm_detail;

    details[j] = odd;
    current_app = even;
  }

  // THRESHOLDING
  std::vector<double> lambdas;
  if (threshold_method == "sure") {
    lambdas = compute_thresholds_sure_internal(details, levels);
  } else {
    lambdas = compute_thresholds_internal(details[0], levels, alpha, beta);
  }

  // SCAD canonical shape parameter (Fan-Li 2001).
  const double SCAD_A = 3.7;

  for (int j = 0; j < levels; j++) {
    double lam = lambdas[j];
    double lam_sq = lam * lam;
    double two_lam = 2.0 * lam;
    double a_lam = SCAD_A * lam;
    double scad_denom = SCAD_A - 2.0;
    std::vector<double> &det = details[j];

    for (int i = 0; i < (int)det.size(); i++) {
      double val = det[i];
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

  // INVERSE ILWT
  int original_len = signal.size();

  for (int j = levels - 1; j >= 0; j--) {
    std::vector<double> &even = current_app;
    std::vector<double> &odd = details[j];

    for (int i = 0; i < (int)even.size(); i++) even[i] /= norm_approx;
    for (int i = 0; i < (int)odd.size();  i++) odd[i]  /= norm_detail;

    std::vector<double> t_even, t_odd;
    if (irregular) {
      const std::vector<double>& t_cur = t_levels[j];
      int n_cur = (int)t_cur.size();
      int n_e = (n_cur + 1) / 2, n_o = n_cur / 2;
      t_even.resize(n_e); t_odd.resize(n_o);
      for (int i = 0; i < n_e; i++) t_even[i] = t_cur[2 * i];
      for (int i = 0; i < n_o; i++) t_odd[i]  = t_cur[2 * i + 1];
    }

    for (int k = (int)cpp_steps.size() - 1; k >= 0; k--) {
      const auto &step = cpp_steps[k];
      int k_filt = (int)step.coeffs.size();
      int n_even = (int)even.size();
      int n_odd  = (int)odd.size();
      const double* c = step.coeffs.data();

      if (step.type == "predict") {
        if (use_os) {
          for (int i = 0; i < n_odd; i++)
            odd[i] += onesided_conv(even, n_even, c, k_filt, step.start_idx, i);
        } else if (irregular && step.degree >= 0) {
          std::vector<double> x_nbr(k_filt), t_nbr(k_filt);
          for (int i = 0; i < n_odd; i++) {
            for (int m = 0; m < k_filt; m++) {
              int idx = i + step.start_idx + m;
              x_nbr[m] = get_val_safe(even, idx, n_even, ext_mode, ll_k);
              t_nbr[m] = get_t_extrap(t_even, idx, n_even);
            }
            odd[i] += interp_predict(x_nbr, t_nbr, k_filt, t_odd[i]);
          }
        } else {
          for (int i = 0; i < n_odd; i++) {
            double sum = 0.0;
            for (int m = 0; m < k_filt; m++)
              sum += get_val_safe(even, i + step.start_idx + m, n_even, ext_mode, ll_k) * c[m];
            odd[i] += sum;
          }
        }
      } else {
        if (use_os) {
          for (int i = 0; i < n_even; i++)
            even[i] -= onesided_conv(odd, n_odd, c, k_filt, step.start_idx, i);
        } else {
          for (int i = 0; i < n_even; i++) {
            double sum = 0.0;
            for (int m = 0; m < k_filt; m++)
              sum += get_val_safe(odd, i + step.start_idx + m, n_odd, ext_mode, ll_k) * c[m];
            even[i] -= sum;
          }
        }
      }
    }

    std::vector<double> merged(even.size() + odd.size());
    for (int i = 0; i < (int)even.size(); i++) merged[2 * i]     = even[i];
    for (int i = 0; i < (int)odd.size();  i++) merged[2 * i + 1] = odd[i];

    current_app = merged;
  }

  if ((int)current_app.size() > original_len)
    current_app.resize(original_len);

  return wrap(current_app);
}
