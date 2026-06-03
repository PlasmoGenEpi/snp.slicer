#ifndef SNP_SLICER_MCMC_KERNEL_COMMON_H
#define SNP_SLICER_MCMC_KERNEL_COMMON_H

#include <Rcpp.h>
#include <Rmath.h>
#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>

namespace snp_slicer {

enum ModelType { POISSON = 0, BINOMIAL = 1, NEGATIVE_BINOMIAL = 2, CATEGORICAL = 3 };

inline double mh_ratio(double logp1, double logp0) {
  if (!R_finite(logp1)) return 0.0;
  if (!R_finite(logp0)) return 1.0;
  const double maxlogp = std::max(logp1, logp0);
  const double e1 = std::exp(logp1 - maxlogp);
  const double e0 = std::exp(logp0 - maxlogp);
  return e1 / (e1 + e0);
}

inline int prop_bin(double prop) {
  if (prop <= 0.0) return 0;
  if (prop <= 0.99) return 1;
  return 2;
}

inline int y_bin(double y) {
  if (ISNA(y)) return -1;
  if (y == 0.0) return 0;
  if (y == 0.5) return 1;
  if (y == 1.0) return 2;
  return -1;
}

inline double dpois_log_native(double y, double lambda) {
  if (ISNA(y)) return 0.0;
  if (lambda <= 0.0) return (y == 0.0) ? 0.0 : R_NegInf;
  return y * std::log(lambda) - lambda - std::lgamma(y + 1.0);
}

inline double dbinom_log_native(double y, double n, double p) {
  if (ISNA(y)) return 0.0;
  if (p < 0.0 || p > 1.0 || n < 0.0) return R_NegInf;
  if (y < 0.0 || y > n) return R_NegInf;
  const double q = 1.0 - p;
  if (p == 0.0) return (y == 0.0) ? 0.0 : R_NegInf;
  if (p == 1.0) return (y == n) ? 0.0 : R_NegInf;
  return std::lgamma(n + 1.0) - std::lgamma(y + 1.0) - std::lgamma(n - y + 1.0)
         + y * std::log(p) + (n - y) * std::log(q);
}

inline double dnbinom_log_native(double y, double r, double prop) {
  if (ISNA(y)) return 0.0;
  const double p = 1.0 / (1.0 + prop);
  const double q = 1.0 - p;
  return std::lgamma(y + r) - std::lgamma(r) - std::lgamma(y + 1.0)
         + r * std::log(p) + y * std::log(q);
}

// loglik_const_* are terms depending only on (y, r), not prop (see prepare_obs_views).
inline double dpois_log_const(double y) {
  if (ISNA(y)) return 0.0;
  return -std::lgamma(y + 1.0);
}

inline double dbinom_log_const(double y, double n) {
  if (ISNA(y)) return 0.0;
  return std::lgamma(n + 1.0) - std::lgamma(y + 1.0) - std::lgamma(n - y + 1.0);
}

inline double dnbinom_log_const(double y, double r) {
  if (ISNA(y)) return 0.0;
  return std::lgamma(y + r) - std::lgamma(r) - std::lgamma(y + 1.0);
}

inline double dpois_log_fast(double prop, double r, double y, double loglik_const) {
  if (ISNA(y)) return 0.0;
  const double lambda = r * prop;
  if (lambda <= 0.0) return (y == 0.0) ? 0.0 : R_NegInf;
  return loglik_const + y * std::log(lambda) - lambda;
}

inline double dbinom_log_fast(double prop, double n, double y, double loglik_const) {
  if (ISNA(y)) return 0.0;
  if (prop < 0.0 || prop > 1.0 || n < 0.0) return R_NegInf;
  if (y < 0.0 || y > n) return R_NegInf;
  const double p = prop;
  const double q = 1.0 - p;
  if (p == 0.0) return (y == 0.0) ? 0.0 : R_NegInf;
  if (p == 1.0) return (y == n) ? 0.0 : R_NegInf;
  return loglik_const + y * std::log(p) + (n - y) * std::log(q);
}

inline double dnbinom_log_fast(double prop, double r, double y, double loglik_const) {
  if (ISNA(y)) return 0.0;
  const double p = 1.0 / (1.0 + prop);
  const double q = 1.0 - p;
  return loglik_const + r * std::log(p) + y * std::log(q);
}

inline double loglik_value_fast(double prop, double y, double r, double loglik_const,
                                int model_type) {
  switch (model_type) {
  case POISSON:
    return dpois_log_fast(prop, r, y, loglik_const);
  case BINOMIAL:
    return dbinom_log_fast(prop, r, y, loglik_const);
  case NEGATIVE_BINOMIAL:
    return dnbinom_log_fast(prop, r, y, loglik_const);
  default:
    return R_NaN;
  }
}

// Count models with y == 0 (obs_code == 1): no ISNA check; drops y*log(q) / y*log(lambda) terms.
inline double count_loglik_y0_fast(double prop, double r, double loglik_const,
                                   int model_type) {
  switch (model_type) {
  case POISSON: {
    const double lambda = r * prop;
    if (lambda <= 0.0) return 0.0;
    return loglik_const - lambda;
  }
  case BINOMIAL: {
    if (prop < 0.0 || prop > 1.0 || r < 0.0) return R_NegInf;
    const double q = 1.0 - prop;
    if (prop == 0.0) return 0.0;
    if (prop == 1.0) return R_NegInf;
    return loglik_const + r * std::log(q);
  }
  case NEGATIVE_BINOMIAL: {
    const double p = 1.0 / (1.0 + prop);
    return loglik_const + r * std::log(p);
  }
  default:
    return R_NaN;
  }
}

inline double loglik_value(double y, double r, double prop, int model_type) {
  switch (model_type) {
  case POISSON:
    return dpois_log_native(y, r * prop);
  case BINOMIAL:
    return dbinom_log_native(y, r, prop);
  case NEGATIVE_BINOMIAL:
    return dnbinom_log_native(y, r, prop);
  default:
    return R_NaN;
  }
}

inline double loglik_categorical_value(double y, double prop,
                                       const double* llik_tab, int llik_tab_nrow) {
  const int pb = prop_bin(prop);
  const int yb = y_bin(y);
  if (yb < 0) return 0.0;
  return llik_tab[pb + yb * llik_tab_nrow];
}

inline double loglik_categorical_value(double y, double prop,
                                       const Rcpp::NumericMatrix& llik_tab) {
  return loglik_categorical_value(y, prop, llik_tab.begin(), llik_tab.nrow());
}

inline double sample_grid_logweights(const std::vector<double>& lw,
                                     const std::vector<double>& xgrid) {
  const int n = static_cast<int>(lw.size());
  if (n == 0) return 0.5;
  double maxlw = R_NegInf;
  for (double w : lw) {
    if (w > maxlw) maxlw = w;
  }
  double sum = 0.0;
  std::vector<double> probs(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    probs[static_cast<std::size_t>(i)] = std::exp(lw[static_cast<std::size_t>(i)] - maxlw);
    sum += probs[static_cast<std::size_t>(i)];
  }
  const double u = unif_rand() * sum;
  double cum = 0.0;
  for (int i = 0; i < n; ++i) {
    cum += probs[static_cast<std::size_t>(i)];
    if (u <= cum) return xgrid[static_cast<std::size_t>(i)];
  }
  return xgrid[static_cast<std::size_t>(n - 1)];
}

inline double gridsample_bounds(double lb, double ub,
                                const std::function<double(double)>& logf) {
  if (!R_finite(lb) || !R_finite(ub) || lb >= ub) return 0.5;
  const int nsample = 100;
  std::vector<double> xgrid;
  std::vector<double> lw;
  xgrid.reserve(nsample - 2);
  lw.reserve(nsample - 2);
  for (int i = 1; i < nsample - 1; ++i) {
    const double x = lb + (ub - lb) * static_cast<double>(i) / static_cast<double>(nsample - 1);
    xgrid.push_back(x);
    lw.push_back(logf(x));
  }
  return sample_grid_logweights(lw, xgrid);
}

inline double logf_newfeature(double x, int N, double alpha) {
  double logf = 0.0;
  for (int i = 1; i <= N; ++i) {
    logf += std::pow(1.0 - x, i) / static_cast<double>(i);
  }
  return alpha * logf + (alpha - 1.0) * std::log(x) + N * std::log(1.0 - x);
}

inline double logf_oldfeature(double x, int m, int N) {
  return (m - 1.0) * std::log(x) + (N - m) * std::log(1.0 - x);
}

}  // namespace snp_slicer

#endif
