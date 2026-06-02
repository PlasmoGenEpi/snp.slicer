#include "mcmc_kernel_common.h"
#include "mcmc_kernel_internal.h"
#include <Rcpp.h>
#include <algorithm>
#include <vector>

using snp_slicer::gridsample_bounds;
using snp_slicer::logf_newfeature;
using snp_slicer::logf_oldfeature;

namespace snp_slicer {
namespace kernel {
namespace {

int get_kstar(const Rcpp::NumericMatrix& A) {
  int result = 0;
  for (int k = 0; k < A.ncol(); ++k) {
    double col_total = 0.0;
    for (int i = 0; i < A.nrow(); ++i) {
      col_total += A(i, k);
    }
    if (col_total > 0.0) result = k + 1;
  }
  return result;
}

double get_mustar(const Rcpp::NumericMatrix& A, const Rcpp::NumericVector& mu) {
  const int kstar = get_kstar(A);
  if (kstar == 0) return 0.0;
  return mu[kstar - 1];
}

void refresh_feature(Rcpp::NumericMatrix& A,
                     Rcpp::NumericMatrix& D,
                     int k,
                     double rho,
                     int P) {
  const int k_idx = k - 1;
  for (int i = 0; i < A.nrow(); ++i) {
    A(i, k_idx) = 0.0;
  }
  for (int p = 0; p < P; ++p) {
    D(k_idx, p) = unif_rand() < rho ? 1.0 : 0.0;
  }
}

Rcpp::NumericVector column_sums(const Rcpp::NumericMatrix& A) {
  Rcpp::NumericVector m(A.ncol());
  for (int k = 0; k < A.ncol(); ++k) {
    double total = 0.0;
    for (int i = 0; i < A.nrow(); ++i) {
      total += A(i, k);
    }
    m[k] = total;
  }
  return m;
}

}  // namespace

void update_s(SliceState& state, const ModelData& model) {
  Rcpp::NumericMatrix& A = state.A;
  Rcpp::NumericMatrix& D = state.D;
  Rcpp::NumericVector& mu = state.mu;
  int& ktrunc = state.ktrunc;
  int& kplus = state.kplus;
  const int N = A.nrow();
  const int P = model.P;
  const double alpha = model.alpha;
  const double rho = model.rho;

  const double mustar = get_mustar(A, mu);
  const double s = mustar * unif_rand();

  // Match R slice_update_s_r: loop while s < mu[k]; when k exceeds length(mu), stop
  // (do not also require k <= mu.size() — k and length grow together and never exit).
  int k = ktrunc;
  const int max_stick_steps = 10000;
  int stick_steps = 0;
  while (true) {
    if (k < 1 || k > static_cast<int>(mu.size())) break;
    if (s >= mu[k - 1]) break;
    if (++stick_steps > max_stick_steps) {
      Rcpp::warning("update_s: stick-breaking expansion exceeded %d steps; truncating",
                    max_stick_steps);
      break;
    }
    const double munext = gridsample_bounds(
      0.0,
      mu[k - 1],
      [&](double x) { return logf_newfeature(x, N, alpha); }
    );
    mu.push_back(munext);
    k++;
  }

  if (static_cast<int>(mu.size()) > ktrunc) {
    const int new_k = static_cast<int>(mu.size());
    const int old_a_cols = A.ncol();
    const int old_d_rows = D.nrow();
    Rcpp::NumericMatrix A_new(N, new_k);
    Rcpp::NumericMatrix D_new(new_k, P);

    const int copy_a = std::min(old_a_cols, new_k);
    const int copy_d = std::min(old_d_rows, new_k);
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < copy_a; ++j) {
        A_new(i, j) = A(i, j);
      }
    }
    for (int kk = 0; kk < copy_d; ++kk) {
      for (int p = 0; p < P; ++p) {
        D_new(kk, p) = D(kk, p);
      }
    }

    A = A_new;
    D = D_new;

    for (int kk = ktrunc + 1; kk <= new_k; ++kk) {
      refresh_feature(A, D, kk, rho, P);
    }
  }

  ktrunc = mu.size();
  kplus = ktrunc;
  for (int i = 0; i < mu.size(); ++i) {
    if (mu[i] < s) {
      kplus = i + 1;
      break;
    }
  }
}

void update_mu(SliceState& state, const ModelData& model) {
  Rcpp::NumericVector& mu = state.mu;
  const int kplus = state.kplus;
  const int ktrunc = state.ktrunc;
  const int N = model.N;
  const double alpha = model.alpha;
  const Rcpp::NumericVector m = column_sums(state.A);

  if (kplus > 1 && mu.size() >= 2 && R_finite(mu[1])) {
    mu[0] = gridsample_bounds(
      mu[1],
      1.0,
      [&](double x) { return logf_oldfeature(x, static_cast<int>(m[0]), N); }
    );
  }

  for (int k = 2; k <= kplus - 1; ++k) {
    const int idx = k - 1;
    if (k + 1 <= static_cast<int>(mu.size()) && k - 1 >= 1 &&
        R_finite(mu[k]) && R_finite(mu[k - 2])) {
      mu[idx] = gridsample_bounds(
        mu[k],
        mu[k - 2],
        [&](double x) { return logf_oldfeature(x, static_cast<int>(m[idx]), N); }
      );
    } else {
      mu[idx] = 0.5;
    }
  }

  for (int k = kplus; k <= ktrunc; ++k) {
    const int idx = k - 1;
    if (idx < 0 || idx >= static_cast<int>(mu.size())) continue;
    if (k - 1 >= 1 && R_finite(mu[k - 2])) {
      mu[idx] = gridsample_bounds(
        0.0,
        mu[k - 2],
        [&](double x) { return logf_newfeature(x, N, alpha); }
      );
    } else {
      mu[idx] = 0.5;
    }
  }
}

}  // namespace kernel
}  // namespace snp_slicer

// [[Rcpp::export]]
Rcpp::List cpp_update_s(Rcpp::NumericMatrix A,
                        Rcpp::NumericMatrix D,
                        Rcpp::NumericVector mu,
                        int ktrunc,
                        double alpha,
                        double rho,
                        int P) {
  snp_slicer::kernel::SliceState state;
  state.A = A;
  state.D = D;
  state.mu = mu;
  state.ktrunc = ktrunc;

  snp_slicer::kernel::ModelData model;
  model.alpha = alpha;
  model.rho = rho;
  model.N = A.nrow();
  model.P = P;

  snp_slicer::kernel::update_s(state, model);

  return Rcpp::List::create(
    Rcpp::_["A"] = state.A,
    Rcpp::_["D"] = state.D,
    Rcpp::_["mu"] = state.mu,
    Rcpp::_["ktrunc"] = state.ktrunc,
    Rcpp::_["kplus"] = state.kplus
  );
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_update_mu(Rcpp::NumericMatrix A,
                                  Rcpp::NumericVector mu,
                                  int kplus,
                                  int ktrunc,
                                  int N,
                                  double alpha) {
  snp_slicer::kernel::SliceState state;
  state.A = A;
  state.mu = mu;
  state.kplus = kplus;
  state.ktrunc = ktrunc;

  snp_slicer::kernel::ModelData model;
  model.N = N;
  model.alpha = alpha;

  snp_slicer::kernel::update_mu(state, model);
  return state.mu;
}
