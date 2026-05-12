#include "utils.h"
using namespace Rcpp;

// Note: Boundary logic (get_val_safe) moved to utils.h

// Apply filter via convolution (Core C++)
 //
 // @param x Input signal.
 // @param coeffs Filter coefficients.
 // @param start_idx Start index offset.
 // @param ext_mode Integer: 1=sym, 2=per, 3=zero, 4=local_linear, 5=one_sided.
 // @keywords internal
 // [[Rcpp::export]]
 NumericVector apply_filter_cpp(
       NumericVector x,
       NumericVector coeffs,
       int start_idx,
       int ext_mode
 ) {
    int n = x.size();
    int k = coeffs.size();
    NumericVector y(n);

    std::vector<double> x_std = as<std::vector<double>>(x);

    if (ext_mode == 5) {  // one_sided: normalized convolution
        const double* c = &coeffs[0];
        for (int i = 0; i < n; i++)
            y[i] = onesided_conv(x_std, n, c, k, start_idx, i);
        return y;
    }

    for (int i = 0; i < n; i++) {
       double sum = 0.0;
       for (int j = 0; j < k; j++) {
          int read_idx = i + start_idx + j;
          sum += get_val_safe(x_std, read_idx, n, ext_mode) * coeffs[j];
       }
       y[i] = sum;
    }

    return y;
 }
