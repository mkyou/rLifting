#include "utils.h"
#include <string>
#include <vector>

using namespace Rcpp;

// Core ILWT (C++)
 //
 // Executes multilevel reconstruction entirely in C++.
 //
 // @param coeffs_list List of coefficients (d1..dn, an).
 // @param steps List of lifting steps.
 // @param norm Normalization vector.
 // @param levels Number of levels.
 // @param ext_mode Boundary mode.
 // @param original_len Original signal length (for final trim).
 // @keywords internal
 // [[Rcpp::export]]
 NumericVector ilwt_cpp(List coeffs_list, List steps, NumericVector norm,
                        int levels, int ext_mode, int original_len) {

     std::string a_name = "a" + std::to_string(levels);
     NumericVector current_app = as<NumericVector>(coeffs_list[a_name]);

     for (int j = levels; j >= 1; j--) {
         std::string d_name = "d" + std::to_string(j);
         NumericVector current_det = as<NumericVector>(coeffs_list[d_name]);

         NumericVector even = clone(current_app);
         NumericVector odd  = clone(current_det);

         for(int i=0; i<even.size(); i++) even[i] /= norm[0];
         for(int i=0; i<odd.size(); i++)  odd[i]  /= norm[1];

         int n_steps = steps.size();
         for (int k = n_steps - 1; k >= 0; k--) {
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
                 for(int m=0; m<len; m++) odd[m] += pred[m];
             }
             else if (type == "update") {
                 NumericVector upd = apply_filter_cpp(
                     odd, coeffs,
                     start_idx, ext_mode
                 );
                 int len = std::min(even.size(), upd.size());
                 for(int m=0; m<len; m++) even[m] -= upd[m];
             }
         }

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

     if (current_app.size() > original_len) {
         NumericVector final_res(original_len);
         for(int i=0; i<original_len; i++) final_res[i] = current_app[i];
         return final_res;
     }

     return current_app;
 }
