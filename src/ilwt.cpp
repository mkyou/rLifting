#include "utils.h"
#include <string>
#include <vector>

using namespace Rcpp;

// ── Core ILWT ─────────────────────────────────────────────────────────────────
// @param coeffs_list List of coefficients (d1..dn, an).
// @param steps List of lifting steps (each may contain 'degree' field).
// @param norm Normalization vector.
// @param levels Number of levels.
// @param ext_mode Boundary mode.
// @param original_len Original signal length (for final trim).
// @param t Optional position vector (length 0 = regular grid).
// @keywords internal
// [[Rcpp::export]]
NumericVector ilwt_cpp(List coeffs_list, List steps, NumericVector norm,
                       int levels, int ext_mode, int original_len,
                       NumericVector t, int ll_k = 2) {

  // Reconstruct per-level position vectors (even sub-sequences, finest first).
  // For irregular grids, t_levels[j] holds positions at decomposition level j.
  bool irregular = (t.size() == (std::size_t)original_len);
  std::vector<std::vector<double>> t_levels;

  if (irregular) {
    std::vector<double> cur = as<std::vector<double>>(t);
    for (int j = 0; j < levels; j++) {
      t_levels.push_back(cur);
      std::vector<double> next;
      for (int i = 0; i < (int)cur.size(); i += 2) next.push_back(cur[i]);
      cur = next;
    }
    // t_levels[j]: positions at level j (finest = 0, coarsest = levels-1)
  }

  std::string a_name = "a" + std::to_string(levels);
  NumericVector current_app = as<NumericVector>(coeffs_list[a_name]);

  for (int j = levels; j >= 1; j--) {
    NumericVector current_det =
        as<NumericVector>(coeffs_list["d" + std::to_string(j)]);

    NumericVector even = clone(current_app);
    NumericVector odd  = clone(current_det);

    for (int i = 0; i < (int)even.size(); i++) even[i] /= norm[0];
    for (int i = 0; i < (int)odd.size();  i++) odd[i]  /= norm[1];

    // Positions at this level (finest = j-1 in t_levels)
    std::vector<double> t_even, t_odd;
    if (irregular) {
      const std::vector<double>& t_cur = t_levels[j - 1];
      int n = (int)t_cur.size();
      int n_e = (n + 1) / 2, n_o = n / 2;
      t_even.resize(n_e); t_odd.resize(n_o);
      for (int i = 0; i < n_e; i++) t_even[i] = t_cur[2 * i];
      for (int i = 0; i < n_o; i++) t_odd[i]  = t_cur[2 * i + 1];
    }

    int n_steps = steps.size();
    for (int k = n_steps - 1; k >= 0; k--) {
      List step = steps[k];
      std::string type = as<std::string>(step["type"]);
      NumericVector coeffs = step["coeffs"];
      int start_idx = step["start_idx"];
      int degree = step.containsElementNamed("degree")
                     ? (int)step["degree"] : -1;

      if (type == "predict") {
        NumericVector pred;
        if (irregular && degree >= 0) {
          int n_odd_sz  = odd.size();
          int n_even_sz = even.size();
          int kf = (int)coeffs.size();
          pred = NumericVector(n_odd_sz);
          const std::vector<double> even_v = as<std::vector<double>>(even);
          std::vector<double> x_nbr(kf), t_nbr(kf);
          for (int i = 0; i < n_odd_sz; i++) {
            for (int m = 0; m < kf; m++) {
              int idx = i + start_idx + m;
              x_nbr[m] = get_val_safe(even_v, idx, n_even_sz, ext_mode, ll_k);
              t_nbr[m] = get_t_extrap(t_even, idx, n_even_sz);
            }
            pred[i] = interp_predict(x_nbr, t_nbr, kf, t_odd[i]);
          }
        } else {
          pred = apply_filter_cpp(even, coeffs, start_idx, ext_mode, ll_k);
        }
        int len = std::min((int)odd.size(), (int)pred.size());
        for (int m = 0; m < len; m++) odd[m] += pred[m];

      } else if (type == "update") {
        NumericVector upd = apply_filter_cpp(odd, coeffs, start_idx, ext_mode, ll_k);
        int len = std::min((int)even.size(), (int)upd.size());
        for (int m = 0; m < len; m++) even[m] -= upd[m];
      }
    }

    int n_total = even.size() + odd.size();
    NumericVector merged(n_total);
    for (int i = 0; i < (int)even.size(); i++)
      if (2*i < n_total) merged[2*i] = even[i];
    for (int i = 0; i < (int)odd.size(); i++)
      if (2*i+1 < n_total) merged[2*i+1] = odd[i];

    current_app = merged;
  }

  if ((int)current_app.size() > original_len) {
    NumericVector res(original_len);
    for (int i = 0; i < original_len; i++) res[i] = current_app[i];
    return res;
  }
  return current_app;
}
