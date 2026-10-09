// Compiled kernel for the multinomial observation model.
//
// The dictionary is an integer matrix of allele codes (strains x targets).
// Observations live in the expanded layout: `counts` is specimens x allele
// slots, `col_offset[p]` is the first slot of target p, and `n_alleles[p]` its
// slot count. Missing genotypes are NA across the whole block of a target.
//
// Categorical draws (new dictionary rows, the Gibbs step over alleles) use
// inversion in index order: one uniform, walk the cumulative probabilities.
// The R reference implementation draws the same way, so the two paths are
// comparable from the same seed.
// [[Rcpp::depends(RcppEigen)]]
#include "mcmc_kernel_common.h"
#include "mcmc_kernel_internal.h"
#include <RcppEigen.h>
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace snp_slicer {
namespace kernel {
namespace {

int sample_categorical_inversion(const double* prob, int M) {
  const double u = unif_rand();
  double cum = 0.0;
  for (int m = 0; m < M; ++m) {
    cum += prob[m];
    if (u <= cum) return m;
  }
  return M - 1;
}

// Draw one dictionary row from the per-target categorical priors.
// prior_pad is max_M x P, column p holding n_alleles[p] probabilities.
void draw_dictionary_row(Rcpp::IntegerMatrix& D, int k_idx,
                         const Rcpp::NumericMatrix& prior_pad,
                         const Rcpp::IntegerVector& n_alleles) {
  const int max_m = prior_pad.nrow();
  const double* prior = prior_pad.begin();
  for (int p = 0; p < D.ncol(); ++p) {
    D(k_idx, p) = sample_categorical_inversion(prior + p * max_m, n_alleles[p]);
  }
}

// Carrier count per specimen and allele slot: (A %*% one-hot(D)).
void carrier_counts(const Rcpp::NumericMatrix& A,
                    const Rcpp::IntegerMatrix& D,
                    const Rcpp::IntegerVector& col_offset,
                    int L,
                    std::vector<double>& M) {
  const int N = A.nrow();
  const int K = A.ncol();
  const int P = D.ncol();
  M.assign(static_cast<std::size_t>(N) * static_cast<std::size_t>(L), 0.0);
  for (int k = 0; k < K; ++k) {
    for (int i = 0; i < N; ++i) {
      const double a = A(i, k);
      if (a <= 0.0) continue;
      for (int p = 0; p < P; ++p) {
        const int slot = col_offset[p] + D(k, p);
        M[static_cast<std::size_t>(i) + static_cast<std::size_t>(slot) * N] += a;
      }
    }
  }
}

}  // namespace

// Slice-variable update for the multinomial model. Same stick-breaking logic
// as update_s(); differs only in how a fresh dictionary row is drawn.
void update_s_multinomial(SliceState& state,
                          Rcpp::IntegerMatrix& D,
                          const Rcpp::NumericMatrix& prior_pad,
                          const Rcpp::IntegerVector& n_alleles,
                          double alpha) {
  Rcpp::NumericMatrix& A = state.A;
  Rcpp::NumericVector& mu = state.mu;
  int& ktrunc = state.ktrunc;
  int& kplus = state.kplus;
  const int N = A.nrow();
  const int P = D.ncol();

  const double mustar = get_mustar(A, mu);
  const double s = mustar * unif_rand();

  int k = ktrunc;
  const int max_stick_steps = 10000;
  int stick_steps = 0;
  while (true) {
    if (k < 1 || k > static_cast<int>(mu.size())) break;
    if (s >= mu[k - 1]) break;
    if (++stick_steps > max_stick_steps) {
      Rcpp::warning("update_s_multinomial: stick-breaking expansion exceeded %d steps; truncating",
                    max_stick_steps);
      break;
    }
    const double munext = gridsample_newfeature(0.0, mu[k - 1], N, alpha);
    mu.push_back(munext);
    k++;
  }

  if (static_cast<int>(mu.size()) > ktrunc) {
    const int new_k = static_cast<int>(mu.size());
    const int copy_a = std::min(A.ncol(), new_k);
    const int copy_d = std::min(D.nrow(), new_k);
    Rcpp::NumericMatrix A_new(N, new_k);
    Rcpp::IntegerMatrix D_new(new_k, P);
    for (int j = 0; j < copy_a; ++j) {
      for (int i = 0; i < N; ++i) A_new(i, j) = A(i, j);
    }
    for (int p = 0; p < P; ++p) {
      for (int kk = 0; kk < copy_d; ++kk) D_new(kk, p) = D(kk, p);
    }
    A = A_new;
    D = D_new;
    for (int kk = ktrunc + 1; kk <= new_k; ++kk) {
      for (int i = 0; i < N; ++i) A(i, kk - 1) = 0.0;
      draw_dictionary_row(D, kk - 1, prior_pad, n_alleles);
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

// Gibbs update of every updatable dictionary cell.
//
// For strain k at target p with leave-k-out carrier counts c[h, m] over the
// specimens h carrying k, option m has log weight
//
//   log prior[m] + sum_h y[h, m] * (log(c[h, m] + 1) - log(c[h, m]))
//
// up to a constant shared by all options (the terms at every other slot and
// the normalisation by each specimen's strain count cancel). If some slot has
// reads in a carrying specimen but no other carrier, only that allele keeps
// the state possible and the cell is set to it without a draw. The carrier
// count matrix is kept incrementally as cells change.
void update_d_multinomial(const Rcpp::NumericMatrix& A,
                          Rcpp::IntegerMatrix& D,
                          const Rcpp::NumericMatrix& counts,
                          const Rcpp::IntegerVector& col_offset,
                          const Rcpp::IntegerVector& n_alleles,
                          const Rcpp::NumericMatrix& log_prior_pad,
                          const Rcpp::NumericMatrix& prior_pad,
                          int kmin,
                          int kstar) {
  const int N = A.nrow();
  const int K = A.ncol();
  const int P = D.ncol();
  const int L = counts.ncol();
  const int max_m = log_prior_pad.nrow();
  if (D.nrow() != K) {
    Rcpp::stop("update_d_multinomial: ncol(A)=%d must equal nrow(D)=%d", K, D.nrow());
  }
  if (kstar < kmin) return;

  std::vector<double> M;
  carrier_counts(A, D, col_offset, L, M);
  const double* y_col = counts.begin();

  std::vector<double> log_ratio(static_cast<std::size_t>(K) + 2, 0.0);
  for (int c = 1; c <= K + 1; ++c) {
    log_ratio[static_cast<std::size_t>(c)] =
      std::log(static_cast<double>(c) + 1.0) - std::log(static_cast<double>(c));
  }

  std::vector<int> hosts;
  std::vector<double> lw(static_cast<std::size_t>(max_m));
  std::vector<double> w(static_cast<std::size_t>(max_m));

  for (int k = kmin; k <= kstar; ++k) {
    const int k_idx = k - 1;
    hosts.clear();
    for (int i = 0; i < N; ++i) {
      if (A(i, k_idx) > 0.0) hosts.push_back(i);
    }
    if (hosts.empty()) {
      draw_dictionary_row(D, k_idx, prior_pad, n_alleles);
      continue;
    }

    for (int p = 0; p < P; ++p) {
      const int Mp = n_alleles[p];
      const int off = col_offset[p];
      const int cur = D(k_idx, p);

      // Alleles that some carrying specimen shows but no other strain of that
      // specimen carries: in a valid state that can only be `cur`.
      int required = -1;
      int n_required = 0;
      for (int m = 0; m < Mp; ++m) {
        const int slot = off + m;
        bool needed = false;
        for (int h : hosts) {
          const double y = y_col[h + slot * N];
          if (ISNA(y) || y == 0.0) continue;
          const double a = A(h, k_idx);
          const double c = M[static_cast<std::size_t>(h) + static_cast<std::size_t>(slot) * N] -
                           (m == cur ? a : 0.0);
          if (c <= 0.0) { needed = true; break; }
        }
        if (needed) { required = m; ++n_required; }
      }
      if (n_required > 1) continue;      // invalid state; leave the cell alone
      int chosen;
      if (n_required == 1) {
        chosen = required;
      } else {
        double max_lw = R_NegInf;
        for (int m = 0; m < Mp; ++m) {
          const int slot = off + m;
          double acc = log_prior_pad(m, p);
          for (int h : hosts) {
            const double y = y_col[h + slot * N];
            if (ISNA(y) || y == 0.0) continue;
            const double a = A(h, k_idx);
            const double c = M[static_cast<std::size_t>(h) + static_cast<std::size_t>(slot) * N] -
                             (m == cur ? a : 0.0);
            acc += y * log_ratio[static_cast<std::size_t>(static_cast<int>(c + 0.5))];
          }
          lw[static_cast<std::size_t>(m)] = acc;
          if (acc > max_lw) max_lw = acc;
        }
        if (!R_finite(max_lw)) continue;
        double total = 0.0;
        for (int m = 0; m < Mp; ++m) {
          w[static_cast<std::size_t>(m)] = std::exp(lw[static_cast<std::size_t>(m)] - max_lw);
          total += w[static_cast<std::size_t>(m)];
        }
        const double u = unif_rand() * total;
        double cum = 0.0;
        chosen = Mp - 1;
        for (int m = 0; m < Mp; ++m) {
          cum += w[static_cast<std::size_t>(m)];
          if (u <= cum) { chosen = m; break; }
        }
      }
      if (chosen != cur) {
        D(k_idx, p) = chosen;
        for (int h : hosts) {
          const double a = A(h, k_idx);
          M[static_cast<std::size_t>(h) + static_cast<std::size_t>(off + cur) * N] -= a;
          M[static_cast<std::size_t>(h) + static_cast<std::size_t>(off + chosen) * N] += a;
        }
      }
    }
  }
}

// Kernel term of the total log-likelihood: sum over specimens and slots with
// reads of y * log(carriers / strains). The multinomial coefficients are a
// constant added on the R side.
double loglik_multinomial(const Rcpp::NumericMatrix& A,
                          const Rcpp::IntegerMatrix& D,
                          const Rcpp::NumericMatrix& counts,
                          const Rcpp::IntegerVector& col_offset) {
  const int N = A.nrow();
  const int K = A.ncol();
  const int L = counts.ncol();
  std::vector<double> M;
  carrier_counts(A, D, col_offset, L, M);
  std::vector<double> an(static_cast<std::size_t>(N), 0.0);
  for (int k = 0; k < K; ++k) {
    for (int i = 0; i < N; ++i) an[static_cast<std::size_t>(i)] += A(i, k);
  }
  const double* y_col = counts.begin();
  double ll = 0.0;
  for (int slot = 0; slot < L; ++slot) {
    for (int i = 0; i < N; ++i) {
      const double y = y_col[i + slot * N];
      if (ISNA(y) || y == 0.0) continue;
      const double prop = M[static_cast<std::size_t>(i) + static_cast<std::size_t>(slot) * N] /
                          an[static_cast<std::size_t>(i)];
      if (prop <= 0.0) return R_NegInf;
      ll += y * std::log(prop);
    }
  }
  return ll;
}

}  // namespace kernel
}  // namespace snp_slicer

// [[Rcpp::export]]
Rcpp::IntegerMatrix cpp_update_d_multinomial(Rcpp::NumericMatrix A,
                                             Rcpp::IntegerMatrix D_codes,
                                             int kmin,
                                             int kstar,
                                             Rcpp::NumericMatrix counts,
                                             Rcpp::IntegerVector col_offset,
                                             Rcpp::IntegerVector n_alleles,
                                             Rcpp::NumericMatrix log_prior_pad,
                                             Rcpp::NumericMatrix prior_pad) {
  Rcpp::RNGScope rng_scope;
  Rcpp::IntegerMatrix D = Rcpp::clone(D_codes);
  snp_slicer::kernel::update_d_multinomial(A, D, counts, col_offset, n_alleles,
                                           log_prior_pad, prior_pad, kmin, kstar);
  return D;
}

// [[Rcpp::export]]
double cpp_loglik_multinomial(Rcpp::NumericMatrix A,
                              Rcpp::IntegerMatrix D_codes,
                              Rcpp::NumericMatrix counts,
                              Rcpp::IntegerVector col_offset) {
  return snp_slicer::kernel::loglik_multinomial(A, D_codes, counts, col_offset);
}

// [[Rcpp::export]]
Rcpp::List cpp_slice_iter_multinomial(Rcpp::NumericMatrix A,
                                      Rcpp::IntegerMatrix D_codes,
                                      Rcpp::NumericVector mu,
                                      Rcpp::IntegerVector mixed,
                                      int kplus,
                                      int kstar,
                                      int kmin,
                                      int ktrunc,
                                      Rcpp::NumericMatrix counts,
                                      Rcpp::IntegerVector col_offset,
                                      Rcpp::IntegerVector n_alleles,
                                      Rcpp::NumericMatrix log_prior_pad,
                                      Rcpp::NumericMatrix prior_pad,
                                      Rcpp::NumericMatrix r_totals,
                                      double alpha,
                                      int N) {
  Rcpp::RNGScope rng_scope;
  snp_slicer::kernel::SliceState state;
  state.A = Rcpp::clone(A);
  state.mu = Rcpp::clone(mu);
  state.kplus = kplus;
  state.kstar = kstar;
  state.kmin = kmin;
  state.ktrunc = ktrunc;
  Rcpp::IntegerMatrix D = Rcpp::clone(D_codes);

  snp_slicer::kernel::ModelData model;
  model.N = N;
  model.alpha = alpha;

  snp_slicer::kernel::update_s_multinomial(state, D, prior_pad, n_alleles, alpha);
  snp_slicer::kernel::update_a_multinomial(state, D, mixed, counts, col_offset, r_totals);
  snp_slicer::kernel::update_d_multinomial(state.A, D, counts, col_offset, n_alleles,
                                           log_prior_pad, prior_pad, state.kmin, state.kstar);
  snp_slicer::kernel::update_mu(state, model);

  return Rcpp::List::create(
    Rcpp::_["A"] = state.A,
    Rcpp::_["D"] = D,
    Rcpp::_["mu"] = state.mu,
    Rcpp::_["kplus"] = state.kplus,
    Rcpp::_["kstar"] = state.kstar,
    Rcpp::_["ktrunc"] = state.ktrunc
  );
}
