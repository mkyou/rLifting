#include <Rcpp.h>
#include <algorithm> // std::nth_element
#include <cmath> // log, sqrt

using namespace Rcpp;

// Helper: Median Absolute Deviation (MAD) in C++
 // @keywords internal
 double calc_mad(NumericVector x) {
    int n = x.size();
    if (n == 0) return 0.0;

    // Copy absolute values
    std::vector<double> abs_x(n);
    for(int i=0; i<n; i++) abs_x[i] = std::abs(x[i]);

    // Selection algorithm (faster than full sort)
    int mid = n / 2;
    std::nth_element(abs_x.begin(), abs_x.begin() + mid, abs_x.end());

    double median = abs_x[mid];

    // For even n, average the two middle elements
    if (n % 2 == 0) {
       std::nth_element(abs_x.begin(), abs_x.begin() + mid - 1, abs_x.end());
       median = (median + abs_x[mid - 1]) / 2.0;
    }

    return median;
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

    double mad_val = calc_mad(d1);
    double sigma = mad_val / 0.6745;

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
