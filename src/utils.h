#ifndef RLIFTING_UTILS_H
#define RLIFTING_UTILS_H

#include <Rcpp.h>
using namespace Rcpp;

// Funcoes utilitarias de convolucao
NumericVector apply_filter_cpp(NumericVector x, NumericVector coeffs, int start_idx, int ext_mode);

// Funcoes Core de Transformada (para uso interno em outros arquivos C++)
List lwt_cpp(
    NumericVector signal, List steps,
    NumericVector norm, int levels, int ext_mode
);

NumericVector ilwt_cpp(
    List coeffs_list, List steps,
    NumericVector norm, int levels, int ext_mode, int original_len
);

NumericVector compute_thresholds_cpp(
    NumericVector d1, int max_level, double alpha, double beta
);

NumericVector threshold_semisoft_cpp(NumericVector x, double lambda);
NumericVector threshold_soft_cpp(NumericVector x, double lambda);
NumericVector threshold_hard_cpp(NumericVector x, double lambda);

#endif
