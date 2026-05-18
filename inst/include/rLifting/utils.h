#ifndef RLIFTING_UTILS_H
#define RLIFTING_UTILS_H

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <string>

// Shared structure for Lifting steps
struct LiftingStep {
    std::string type;
    std::vector<double> coeffs;
    int start_idx;
    int degree = -1;  // -1 = fixed coeffs (regular grid only)
                      //  0 = Haar (nearest-neighbour)
                      //  1 = linear interpolation
                      //  3 = cubic interpolation (Lagrange, 4-point)
};

// Position-aware interpolation for irregular-grid predict steps.
// Returns the predicted value at position t_target given k neighbours
// at positions t_nbr[] with signal values x_nbr[].
// degree 0: nearest left neighbour (Haar)
// degree 1: linear interpolation (2 points)
// degree 3: Lagrange cubic interpolation (4 points)
inline double interp_predict(
        const std::vector<double>& x_nbr,
        const std::vector<double>& t_nbr,
        int k, double t_target
) {
    if (k == 0) return 0.0;
    if (k == 1) return x_nbr[0];            // degree 0

    if (k == 2) {                            // degree 1: linear
        double span = t_nbr[1] - t_nbr[0];
        if (std::abs(span) < 1e-15) return 0.5 * (x_nbr[0] + x_nbr[1]);
        double w1 = (t_target - t_nbr[0]) / span;
        return x_nbr[0] * (1.0 - w1) + x_nbr[1] * w1;
    }

    // degree 3: 4-point Lagrange (k == 4)
    double result = 0.0;
    for (int i = 0; i < k; i++) {
        double w = 1.0;
        for (int j = 0; j < k; j++) {
            if (j != i) {
                double denom = t_nbr[i] - t_nbr[j];
                if (std::abs(denom) < 1e-15) { w = 0.0; break; }
                w *= (t_target - t_nbr[j]) / denom;
            }
        }
        result += w * x_nbr[i];
    }
    return result;
}

// Inline Helper Functions
// Centralized boundary logic to be shared between Online, Offline,
// and Utils engines.
inline double get_val_safe(
        const std::vector<double>& x,
        int i, int n, int mode, int ll_k = 2
) {
    // Mode: 1=symmetric, 2=periodic, 3=zero, 4=local_linear

    if (i >= 0 && i < n) return x[i]; // Fast path

    if (mode == 3) return 0.0; // Zero

    if (mode == 2) { // Periodic
        int idx = i % n;
        if (idx < 0) idx += n;
        return x[idx];
    }

    if (mode == 4) { // Local linear extrapolation (k-point OLS)
        if (n < 2) return x[0];
        int k = std::min(ll_k, n);
        if (k < 2) k = 2;
        if (i < 0) {
            double st = 0, sy = 0, stt = 0, sty = 0;
            for (int j = 0; j < k; j++) {
                st += j; sy += x[j];
                stt += (double)j * j; sty += (double)j * x[j];
            }
            double denom = k * stt - st * st;
            if (std::abs(denom) < 1e-15) return x[0];
            double slope = (k * sty - st * sy) / denom;
            double intercept = (sy - slope * st) / k;
            return intercept + slope * (double)i;
        } else {
            double st = 0, sy = 0, stt = 0, sty = 0;
            for (int j = 0; j < k; j++) {
                double t = (double)(n - k + j);
                double y = x[n - k + j];
                st += t; sy += y; stt += t * t; sty += t * y;
            }
            double denom = k * stt - st * st;
            if (std::abs(denom) < 1e-15) return x[n - 1];
            double slope = (k * sty - st * sy) / denom;
            double intercept = (sy - slope * st) / k;
            return intercept + slope * (double)i;
        }
    }

    // Symmetric (Reflection) — default
    while (i < 0 || i >= n) {
        if (i < 0) i = -1 - i;
        else i = 2 * n - 1 - i;
    }
    if (i < 0) i = 0;
    if (i >= n) i = n - 1;
    return x[i];
}

// Returns the extrapolated t position for an out-of-bounds index.
// Uses linear extrapolation from the nearest boundary spacing.
inline double get_t_extrap(const std::vector<double>& t, int idx, int n) {
    if (idx >= 0 && idx < n) return t[idx];
    if (n < 2) return t[0];
    if (idx < 0) {
        double dt = t[1] - t[0];
        return t[0] + (double)idx * dt;
    }
    double dt = t[n-1] - t[n-2];
    return t[n-1] + (double)(idx - (n-1)) * dt;
}

// Normalized (one-sided) convolution at position `pos`.
// Uses only in-bounds coefficients, renormalized by their sum.
// Fast path when the entire filter window is within bounds.
inline double onesided_conv(
        const std::vector<double>& x, int n,
        const double* c, int k, int start_idx, int pos
) {
    int first = pos + start_idx;
    int last  = first + k - 1;
    if (first >= 0 && last < n) {   // fast path: no boundary contact
        double s = 0.0;
        for (int j = 0; j < k; j++) s += x[first + j] * c[j];
        return s;
    }
    double vs = 0.0, ws = 0.0;     // boundary path: drop out-of-bounds taps
    for (int j = 0; j < k; j++) {
        int idx = first + j;
        if (idx >= 0 && idx < n) { vs += x[idx] * c[j]; ws += c[j]; }
    }
    return (ws > 1e-15) ? vs / ws : 0.0;
}

// Function Signatures

Rcpp::NumericVector apply_filter_cpp(
        Rcpp::NumericVector x, Rcpp::NumericVector coeffs,
        int start_idx, int ext_mode, int ll_k = 2
);

Rcpp::List lwt_cpp(
        Rcpp::NumericVector signal, Rcpp::List steps,
        Rcpp::NumericVector norm,
        int levels, int ext_mode, Rcpp::NumericVector t, int ll_k
);

Rcpp::NumericVector ilwt_cpp(
        Rcpp::List coeffs_list, Rcpp::List steps,
        Rcpp::NumericVector norm,
        int levels, int ext_mode, int original_len,
        Rcpp::NumericVector t, int ll_k
);

Rcpp::NumericVector compute_thresholds_cpp(
        Rcpp::NumericVector d1, int max_level,
        double alpha, double beta
);

// SURE-optimal threshold for a single decomposition level
// (Donoho & Johnstone 1995). Returns the lambda minimising SURE for
// soft thresholding, clamped to sigma * sqrt(2 log n). Sigma is
// estimated per-level via MAD.
inline double compute_sure_lambda_level(const std::vector<double> &d) {
    int n = d.size();
    if (n == 0) return 0.0;

    std::vector<double> abs_d(n);
    for (int i = 0; i < n; i++) abs_d[i] = std::abs(d[i]);

    std::vector<double> abs_for_mad = abs_d;
    int mid = n / 2;
    std::nth_element(abs_for_mad.begin(), abs_for_mad.begin() + mid,
                     abs_for_mad.end());
    double sigma = abs_for_mad[mid] / 0.6745;
    if (sigma < 1e-15) return 0.0;

    std::sort(abs_d.begin(), abs_d.end());
    double sigma_sq = sigma * sigma;
    double universal_cap = sigma * std::sqrt(2.0 * std::log((double)n));

    std::vector<double> cumsum_sq(n + 1, 0.0);
    for (int i = 0; i < n; i++)
        cumsum_sq[i + 1] = cumsum_sq[i] + abs_d[i] * abs_d[i];

    double n_sigma_sq = (double)n * sigma_sq;
    double best_sure = n_sigma_sq;
    double best_lambda = 0.0;

    for (int k = 0; k < n; k++) {
        double lam = abs_d[k];
        double sure = n_sigma_sq + cumsum_sq[k + 1] +
                      lam * lam * (double)(n - k - 1) -
                      2.0 * sigma_sq * (double)(k + 1);
        if (sure < best_sure) {
            best_sure = sure;
            best_lambda = lam;
        }
    }
    if (best_lambda > universal_cap) best_lambda = universal_cap;
    return best_lambda;
}

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
        int ext_mode,
        Rcpp::NumericVector t,
        int ll_k
);

#endif
