#include "utils.h"
#include <string>
#include <vector>

using namespace Rcpp;

// Core LWT (C++)
 //
 // Executes multilevel decomposition entirely in C++.
 //
 // @param signal Input signal.
 // @param steps List of lifting steps.
 // @param norm Normalization vector.
 // @param levels Number of levels.
 // @param ext_mode Boundary mode.
 // @keywords internal
 // [[Rcpp::export]]
 List lwt_cpp(
       NumericVector signal, List steps,
       NumericVector norm, int levels,
       int ext_mode
 ) {

    List coeffs_out;
    NumericVector current_app = clone(signal);

    for (int j = 1; j <= levels; j++) {
       int n = current_app.size();

       int n_even = (n + 1) / 2;
       int n_odd  = n / 2;

       NumericVector even(n_even);
       NumericVector odd(n_odd);

       for (int i = 0; i < n_even; i++) even[i] = current_app[2*i];
       for (int i = 0; i < n_odd; i++)  odd[i]  = current_app[2*i + 1];

       int n_steps = steps.size();

       for (int k = 0; k < n_steps; k++) {
          List step = steps[k];
          std::string type = as<std::string>(step["type"]);
          NumericVector coeffs = step["coeffs"];
          int start_idx = step["start_idx"];

          if (type == "predict") {
             NumericVector pred = apply_filter_cpp(
                even, coeffs,
                start_idx, ext_mode
             );

             int len = std::min(odd.size(), pred.size());
             for(int m=0; m<len; m++) odd[m] -= pred[m];

          } else if (type == "update") {
             NumericVector upd = apply_filter_cpp(
                odd, coeffs,
                start_idx, ext_mode
             );

             int len = std::min(even.size(), upd.size());
             for(int m=0; m<len; m++) even[m] += upd[m];
          }
       }

       for(int m=0; m<even.size(); m++) even[m] *= norm[0];
       for(int m=0; m<odd.size(); m++)  odd[m]  *= norm[1];

       std::string d_name = "d" + std::to_string(j);
       coeffs_out[d_name] = odd;

       current_app = even;
    }

    std::string a_name = "a" + std::to_string(levels);
    coeffs_out[a_name] = current_app;

    return coeffs_out;
 }
