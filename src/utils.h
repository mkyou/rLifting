#ifndef RLIFTING_UTILS_H
#define RLIFTING_UTILS_H

#include <Rcpp.h>
#include <vector>
#include <string>

// Shared structure for Lifting steps
struct LiftingStep {
        std::string type;
        std::vector<double> coeffs;
        int start_idx;
};

// Inline Helper Functions
// Centralized boundary logic to be shared between Online, Offline,
// and Utils engines.
inline double get_val_safe(
                const std::vector<double>& x,
                int i, int n, int mode
) {
        // Mode: 1=symmetric, 2=periodic, 3=zero

        if (i >= 0 && i < n) return x[i]; // Fast path

        if (mode == 3) return 0.0; // Zero

        if (mode == 2) { // Periodic
                int idx = i % n;
                if (idx < 0) idx += n;
                return x[idx];
        }

        // Symmetric (Reflection)
        while (i < 0 || i >= n) {
                if (i < 0) i = -1 - i;
                else i = 2 * n - 1 - i;
        }
        if (i < 0) i = 0;
        if (i >= n) i = n - 1;
        return x[i];
}

// Function Signatures

Rcpp::NumericVector apply_filter_cpp(
                Rcpp::NumericVector x, Rcpp::NumericVector coeffs,
                int start_idx, int ext_mode
);

Rcpp::List lwt_cpp(
                Rcpp::NumericVector signal, Rcpp::List steps,
                Rcpp::NumericVector norm,
                int levels, int ext_mode
);

Rcpp::NumericVector ilwt_cpp(
                Rcpp::List coeffs_list, Rcpp::List steps,
                Rcpp::NumericVector norm,
                int levels, int ext_mode, int original_len
);

Rcpp::NumericVector compute_thresholds_cpp(
                Rcpp::NumericVector d1, int max_level,
                double alpha, double beta
);

Rcpp::NumericVector threshold_semisoft_cpp(
                Rcpp::NumericVector x,
                double lambda
);
Rcpp::NumericVector threshold_soft_cpp(Rcpp::NumericVector x, double lambda);
Rcpp::NumericVector threshold_hard_cpp(Rcpp::NumericVector x, double lambda);

Rcpp::NumericVector denoise_offline_cpp(
                Rcpp::NumericVector signal,
                Rcpp::List steps,
                Rcpp::NumericVector norm,
                int levels,
                double alpha,
                double beta,
                std::string method,
                int ext_mode
);

#endif
