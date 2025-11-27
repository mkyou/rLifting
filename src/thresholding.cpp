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
