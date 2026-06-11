#include "utils.h"
#include <Rcpp.h>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

// Median Absolute Deviation (canonical) exposed to R for testing
// and downstream use. Delegates to compute_mad in inst/include/rLifting/utils.h.
// [[Rcpp::export]]
double compute_mad_cpp(NumericVector x) {
    std::vector<double> vals(x.begin(), x.end());
    return compute_mad(vals);
}

 // Adaptive Threshold Calculation (C++)
 //
 // Computes recursive thresholds directly in C++.
 //
 // @param d1 Detail coefficients from level 1.
 // @param max_level Maximum number of levels.
 // @param alpha Decay parameter.
 // @param beta Scale parameter.
 // @keywords internal
 // [[Rcpp::export]]
 NumericVector compute_thresholds_cpp(
       NumericVector d1, int max_level,
       double alpha, double beta
 ) {
    NumericVector lambdas(max_level);

    std::vector<double> d1_vals(d1.begin(), d1.end());
    double sigma = compute_mad(d1_vals) / 0.6745;

    if (sigma < 1e-15) {
       return lambdas;
    }

    int n1 = d1.size();
    double lambda_1 = beta * sigma * std::sqrt(2.0 * std::log((double)n1));
    lambdas[0] = lambda_1;

    for (int k = 1; k < max_level; k++) {
       int level = k + 1;
       double prev = lambdas[k-1];
       double factor = (double)(level - 1) / (double)(level + alpha - 1);
       lambdas[k] = prev * factor;
    }

    return lambdas;
 }
