#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Hard Thresholding (C++)
 //
 // Sets coefficients below threshold to zero.
 //
 // @param x Coefficient vector.
 // @param lambda Threshold value.
 // @keywords internal
 // [[Rcpp::export]]
 NumericVector threshold_hard_cpp(NumericVector x, double lambda) {
    int n = x.size();
    NumericVector y(n);

    for(int i = 0; i < n; i++) {
       if (std::abs(x[i]) >= lambda) {
          y[i] = x[i];
       } else {
          y[i] = 0.0;
       }
    }
    return y;
 }

 // Soft Thresholding (C++)
 //
 // Shrinks coefficients towards zero.
 //
 // @param x Coefficient vector.
 // @param lambda Threshold value.
 // @keywords internal
 // [[Rcpp::export]]
 NumericVector threshold_soft_cpp(NumericVector x, double lambda) {
    int n = x.size();
    NumericVector y(n);

    for(int i = 0; i < n; i++) {
       double abs_x = std::abs(x[i]);
       if (abs_x < lambda) {
          y[i] = 0.0;
       } else {
          if (x[i] > 0) {
             y[i] = abs_x - lambda;
          } else {
             y[i] = -(abs_x - lambda);
          }
       }
    }
    return y;
 }

 // Semisoft Shrinkage (Hyperbolic C++)
 //
 // Implementation of Liu et al. (2014).
 //
 // @param x Coefficient vector.
 // @param lambda Threshold value.
 // @keywords internal
 // [[Rcpp::export]]
 NumericVector threshold_semisoft_cpp(NumericVector x, double lambda) {
    int n = x.size();
    NumericVector y(n);
    double lambda_sq = lambda * lambda;

    for(int i = 0; i < n; i++) {
       double val = x[i];
       double abs_val = std::abs(val);

       if (abs_val < lambda) {
          y[i] = 0.0;
       } else {
          double shrunk = std::sqrt(val * val - lambda_sq);
          if (val < 0) shrunk = -shrunk;
          y[i] = shrunk;
       }
    }
    return y;
 }

 // SCAD Shrinkage (Antoniadis & Fan 2001, Fan-Li 2001)
 //
 // Three-region shrinkage with continuity at lambda, 2*lambda and a*lambda.
 // Identity for |x| > a*lambda removes the soft-threshold bias on large
 // coefficients while preserving the sparsity-inducing zero region.
 //
 // @param x Coefficient vector.
 // @param lambda Positive threshold value.
 // @param a Shape parameter, a > 2. Fan-Li canonical value: 3.7.
 // @keywords internal
 // [[Rcpp::export]]
 NumericVector threshold_scad_cpp(NumericVector x, double lambda, double a) {
    if (a <= 2.0) {
       stop("SCAD parameter 'a' must be strictly greater than 2.");
    }
    int n = x.size();
    NumericVector y(n);
    double two_lambda = 2.0 * lambda;
    double a_lambda = a * lambda;
    double denom = a - 2.0;

    for (int i = 0; i < n; i++) {
       double val = x[i];
       double abs_val = std::abs(val);
       double sgn = (val > 0.0) ? 1.0 : (val < 0.0 ? -1.0 : 0.0);

       if (abs_val <= lambda) {
          y[i] = 0.0;
       } else if (abs_val <= two_lambda) {
          y[i] = sgn * (abs_val - lambda);
       } else if (abs_val <= a_lambda) {
          y[i] = ((a - 1.0) * val - sgn * a_lambda) / denom;
       } else {
          y[i] = val;
       }
    }
    return y;
 }
