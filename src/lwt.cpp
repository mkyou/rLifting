#include "utils.h"
#include <string>
#include <vector>

using namespace Rcpp;

// ── Irregular-grid predict helper ─────────────────────────────────────────────
// Applies one predict step using position-aware interpolation.
// Falls back to apply_filter_cpp when degree < 0 (fixed coefficients).
static NumericVector predict_irregular(
    const NumericVector& even,
    const std::vector<double>& t_even,
    const std::vector<double>& t_odd,
    const NumericVector& coeffs,
    int start_idx,
    int degree,
    int ext_mode,
    int ll_k
) {
  int n_odd  = (int)t_odd.size();
  int n_even = (int)even.size();
  int k      = (int)coeffs.size();
  NumericVector pred(n_odd);

  const std::vector<double> even_v = as<std::vector<double>>(even);
  std::vector<double> x_nbr(k), t_nbr(k);

  for (int i = 0; i < n_odd; i++) {
    for (int j = 0; j < k; j++) {
      int idx = i + start_idx + j;
      x_nbr[j] = get_val_safe(even_v, idx, n_even, ext_mode, ll_k);
      t_nbr[j] = get_t_extrap(t_even, idx, n_even);
    }
    pred[i] = interp_predict(x_nbr, t_nbr, k, t_odd[i]);
  }
  return pred;
}

// ── Core LWT ──────────────────────────────────────────────────────────────────
// @param signal Input signal.
// @param steps  List of lifting steps (each may contain 'degree' field).
// @param norm   Normalization vector.
// @param levels Number of levels.
// @param ext_mode Boundary mode.
// @param t Optional position vector (length 0 = regular grid).
// @keywords internal
// [[Rcpp::export]]
List lwt_cpp(NumericVector signal, List steps, NumericVector norm, int levels,
             int ext_mode, NumericVector t, int ll_k = 2) {

  bool irregular = (t.size() == signal.size());
  List coeffs_out;
  NumericVector current_app = clone(signal);
  std::vector<double> current_t;
  if (irregular) current_t = as<std::vector<double>>(t);

  for (int j = 1; j <= levels; j++) {
    int n = current_app.size();
    int n_even = (n + 1) / 2;
    int n_odd  = n / 2;

    NumericVector even(n_even), odd(n_odd);
    for (int i = 0; i < n_even; i++) even[i] = current_app[2 * i];
    for (int i = 0; i < n_odd;  i++) odd[i]  = current_app[2 * i + 1];

    std::vector<double> t_even, t_odd;
    if (irregular) {
      t_even.resize(n_even); t_odd.resize(n_odd);
      for (int i = 0; i < n_even; i++) t_even[i] = current_t[2 * i];
      for (int i = 0; i < n_odd;  i++) t_odd[i]  = current_t[2 * i + 1];
    }

    int n_steps = steps.size();
    for (int k = 0; k < n_steps; k++) {
      List step = steps[k];
      std::string type = as<std::string>(step["type"]);
      NumericVector coeffs = step["coeffs"];
      int start_idx = step["start_idx"];
      int degree = step.containsElementNamed("degree")
                     ? (int)step["degree"] : -1;

      if (type == "predict") {
        NumericVector pred;
        if (irregular && degree >= 0) {
          pred = predict_irregular(even, t_even, t_odd,
                                   coeffs, start_idx, degree, ext_mode, ll_k);
        } else {
          pred = apply_filter_cpp(even, coeffs, start_idx, ext_mode, ll_k);
        }
        int len = std::min((int)odd.size(), (int)pred.size());
        for (int m = 0; m < len; m++) odd[m] -= pred[m];

      } else if (type == "update") {
        NumericVector upd = apply_filter_cpp(odd, coeffs, start_idx, ext_mode, ll_k);
        int len = std::min((int)even.size(), (int)upd.size());
        for (int m = 0; m < len; m++) even[m] += upd[m];
      }
    }

    for (int m = 0; m < (int)even.size(); m++) even[m] *= norm[0];
    for (int m = 0; m < (int)odd.size();  m++) odd[m]  *= norm[1];

    coeffs_out["d" + std::to_string(j)] = odd;
    current_app = even;
    if (irregular) current_t = t_even;
  }

  coeffs_out["a" + std::to_string(levels)] = current_app;
  return coeffs_out;
}
