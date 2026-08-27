// [[Rcpp::depends(RcppEigen)]]
#include "mcmc_kernel_common.h"
#include "mcmc_kernel_internal.h"
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif
#include <RcppEigen.h>
#include <Rcpp.h>
#include <algorithm>
#include <vector>

using snp_slicer::CATEGORICAL;
using snp_slicer::NEGATIVE_BINOMIAL;
using snp_slicer::mh_ratio;
using snp_slicer::loglik_value;
using snp_slicer::loglik_value_fast;
using snp_slicer::count_loglik_y0_fast;
using snp_slicer::loglik_categorical_value;

namespace snp_slicer {
namespace kernel {
namespace {

double loglik_target_col_active(const std::vector<int>& active_hosts,
                                const double* ad_col,
                                const double* A_col,
                                double d_kp,
                                const int* an,
                                const ObsViews& obs,
                                int target_idx,
                                int model_type,
                                bool use_prop1) {
  const int p = target_idx;
  double ll = 0.0;
  if (model_type == CATEGORICAL) {
    for (int i : active_hosts) {
      const double a_val = A_col[i];
      const double ad0_i = ad_col[i] - a_val * d_kp;
      const double prop = use_prop1
        ? (ad0_i + a_val) / static_cast<double>(an[i])
        : ad0_i / static_cast<double>(an[i]);
      ll += loglik_categorical_value(
        obs.y_at(i, p), prop, obs.llik_tab, obs.llik_tab_nrow);
    }
    return ll;
  }
  for (int i : active_hosts) {
    const double a_val = A_col[i];
    const double ad0_i = ad_col[i] - a_val * d_kp;
    const double an_i = static_cast<double>(an[i]);
    const double prop = use_prop1
      ? (ad0_i + a_val) / an_i
      : ad0_i / an_i;
    if (obs.has_obs_code()) {
      const int code = obs.obs_code_at(i, p);
      if (code == OBS_SKIP) continue;
      if (obs.has_loglik_const()) {
        const double c = obs.loglik_const_at(i, p);
        const double r = obs.r_at(i, p);
        if (code == OBS_Y_ZERO) {
          ll += count_loglik_y0_fast(prop, r, c, model_type);
        } else {
          ll += loglik_value_fast(prop, obs.y_at(i, p), r, c, model_type);
        }
      } else {
        ll += loglik_value(obs.y_at(i, p), obs.r_at(i, p), prop, model_type);
      }
    } else {
      const double y = obs.y_at(i, p);
      const double r = obs.r_at(i, p);
      if (obs.has_loglik_const()) {
        const double c = obs.loglik_const_at(i, p);
        ll += loglik_value_fast(prop, y, r, c, model_type);
      } else {
        ll += loglik_value(y, r, prop, model_type);
      }
    }
  }
  return ll;
}

void build_active_hosts_by_col(const Eigen::MatrixXd& A_e,
                               int kmin,
                               int kstar,
                               std::vector<std::vector<int>>& active_hosts) {
  const int N = static_cast<int>(A_e.rows());
  const int K = static_cast<int>(A_e.cols());
  if (static_cast<int>(active_hosts.size()) != K) {
    active_hosts.assign(static_cast<std::size_t>(K), {});
  }
  for (int k = kmin; k <= kstar; ++k) {
    const int k_idx = k - 1;
    std::vector<int>& hosts = active_hosts[static_cast<std::size_t>(k_idx)];
    hosts.clear();
    const double* A_col = A_e.col(k_idx).data();
    for (int i = 0; i < N; ++i) {
      if (A_col[i] > 0.0) hosts.push_back(i);
    }
  }
}

struct ColSumTracker {
  std::vector<double> sums;
  std::vector<int> active_cols;

  void init(const Eigen::MatrixXd& A_e) {
    const int K = static_cast<int>(A_e.cols());
    sums.assign(static_cast<std::size_t>(K), 0.0);
    active_cols.clear();
    active_cols.reserve(static_cast<std::size_t>(K));
    for (int kk = 0; kk < K; ++kk) {
      sums[static_cast<std::size_t>(kk)] = A_e.col(kk).sum();
      if (sums[static_cast<std::size_t>(kk)] > 0.0) {
        active_cols.push_back(kk);
      }
    }
  }

  void apply_delta(int k_idx, double delta) {
    const double old = sums[static_cast<std::size_t>(k_idx)];
    sums[static_cast<std::size_t>(k_idx)] = old + delta;
    const double updated = sums[static_cast<std::size_t>(k_idx)];
    if (old == 0.0 && updated > 0.0) {
      auto it = std::lower_bound(active_cols.begin(), active_cols.end(), k_idx);
      active_cols.insert(it, k_idx);
    } else if (old > 0.0 && updated == 0.0) {
      auto it = std::find(active_cols.begin(), active_cols.end(), k_idx);
      if (it != active_cols.end()) active_cols.erase(it);
    }
  }

  int second_last_active_col() const {
    if (active_cols.empty()) return -1;
    if (active_cols.size() == 1) return active_cols[0];
    return active_cols[active_cols.size() - 2];
  }

  int kstar_from_active() const {
    return active_cols.empty() ? 0 : active_cols.back() + 1;
  }
};

}  // namespace

void update_a(SliceState& state, ModelData& model) {
  Rcpp::NumericMatrix& A = state.A;
  const Rcpp::NumericMatrix& D = state.D;
  const Rcpp::NumericVector& mu = state.mu;
  const Rcpp::IntegerVector& mixed = model.mixed;
  int& kstar = state.kstar;
  const int N = A.nrow();
  const int K = A.ncol();
  const int P = D.ncol();
  const int mu_len = mu.size();
  int kplus = state.kplus;
  if (kplus > K) kplus = K;
  if (kplus > mu_len) kplus = mu_len;
  if (kplus < 0) kplus = 0;
  const ObsViews& obs = model.obs;

  Eigen::Map<Eigen::MatrixXd> A_e(A.begin(), N, K);
  Eigen::Map<const Eigen::MatrixXd> D_e(D.begin(), K, P);
  const double* D_ptr = D.begin();

  std::vector<double> y_row(static_cast<std::size_t>(P));
  std::vector<double> r_row(static_cast<std::size_t>(P));
  std::vector<double> c_row(static_cast<std::size_t>(P));
  std::vector<int> code_row(static_cast<std::size_t>(P));
  Eigen::RowVectorXd full_row(P);

  ColSumTracker col_tracker;
  col_tracker.init(A_e);

  std::vector<double> log_mu(static_cast<std::size_t>(mu_len));
  std::vector<double> log_1m_mu(static_cast<std::size_t>(mu_len));
  for (int k = 0; k < mu_len; ++k) {
    log_mu[static_cast<std::size_t>(k)] = std::log(mu[k]);
    log_1m_mu[static_cast<std::size_t>(k)] = std::log(1.0 - mu[k]);
  }

  for (int m = 0; m < mixed.size(); ++m) {
    const int host = mixed[m] - 1;

    full_row = A_e.row(host) * D_e;
    double* full_ptr = full_row.data();
    double row_sum = A_e.row(host).sum();

    const bool use_obs_layout =
      model.model_type != CATEGORICAL &&
      obs.has_loglik_const() &&
      obs.has_obs_code();

    if (use_obs_layout) {
      for (int p = 0; p < P; ++p) {
        y_row[static_cast<std::size_t>(p)] = obs.y_at(host, p);
        r_row[static_cast<std::size_t>(p)] = obs.r_at(host, p);
        c_row[static_cast<std::size_t>(p)] = obs.loglik_const_at(host, p);
        code_row[static_cast<std::size_t>(p)] = obs.obs_code_at(host, p);
      }
    } else {
      for (int p = 0; p < P; ++p) {
        y_row[static_cast<std::size_t>(p)] = obs.y_at(host, p);
        r_row[static_cast<std::size_t>(p)] = obs.r_at(host, p);
      }
    }

    for (int k = 1; k <= kplus; ++k) {
      const int k_idx = k - 1;
      const double old_a = A(host, k_idx);
      const double a0 = row_sum - old_a;

      if (a0 == 0.0) {
        A(host, k_idx) = 1.0;
        if (A(host, k_idx) != old_a) {
          const double delta = A(host, k_idx) - old_a;
          full_row += delta * D_e.row(k_idx);
          row_sum += delta;
          col_tracker.apply_delta(k_idx, delta);
          if (A(host, k_idx) == 1.0) {
            kstar = std::max(kstar, k);
          } else if (k == kstar && col_tracker.sums[static_cast<std::size_t>(k_idx)] == 0.0) {
            kstar = col_tracker.kstar_from_active();
          }
        }
        continue;
      }

      double logp0 = log_1m_mu[static_cast<std::size_t>(k_idx)];
      double logp1 = log_mu[static_cast<std::size_t>(k_idx)];

      if (model.model_type == CATEGORICAL) {
        for (int p = 0; p < P; ++p) {
          if (obs.has_obs_code() && obs.obs_code_at(host, p) == OBS_SKIP) {
            continue;
          }
          const double d_kp = D_ptr[k_idx + p * K];
          const double ad0 = full_ptr[p] - old_a * d_kp;
          const double prop0 = ad0 / a0;
          const double prop1 = (ad0 + d_kp) / (a0 + 1.0);
          const double yp = obs.has_obs_code()
            ? obs.y_at(host, p)
            : y_row[static_cast<std::size_t>(p)];
          logp0 += loglik_categorical_value(yp, prop0, obs.llik_tab, obs.llik_tab_nrow);
          logp1 += loglik_categorical_value(yp, prop1, obs.llik_tab, obs.llik_tab_nrow);
        }
      } else if (use_obs_layout && model.model_type == NEGATIVE_BINOMIAL) {
        const double inv_a0 = 1.0 / a0;
        const double inv_a1 = 1.0 / (a0 + 1.0);
        for (int p = 0; p < P; ++p) {
          const int code = code_row[static_cast<std::size_t>(p)];
          if (code == OBS_SKIP) continue;
          const double d_kp = D_ptr[k_idx + p * K];
          const double ad0 = full_ptr[p] - old_a * d_kp;
          const double prop0 = ad0 * inv_a0;
          const double prop1 = (ad0 + d_kp) * inv_a1;
          const double c = c_row[static_cast<std::size_t>(p)];
          const double rp = r_row[static_cast<std::size_t>(p)];
          if (code == OBS_Y_ZERO) {
            const double p0 = 1.0 / (1.0 + prop0);
            const double p1 = 1.0 / (1.0 + prop1);
            logp0 += c + rp * std::log(p0);
            logp1 += c + rp * std::log(p1);
          } else {
            const double yp = y_row[static_cast<std::size_t>(p)];
            const double p0 = 1.0 / (1.0 + prop0);
            const double p1 = 1.0 / (1.0 + prop1);
            logp0 += c + rp * std::log(p0) + yp * std::log(1.0 - p0);
            logp1 += c + rp * std::log(p1) + yp * std::log(1.0 - p1);
          }
        }
      } else if (use_obs_layout) {
        const double inv_a0 = 1.0 / a0;
        const double inv_a1 = 1.0 / (a0 + 1.0);
        const int model_type = model.model_type;
        for (int p = 0; p < P; ++p) {
          const int code = code_row[static_cast<std::size_t>(p)];
          if (code == OBS_SKIP) continue;
          const double d_kp = D_ptr[k_idx + p * K];
          const double ad0 = full_ptr[p] - old_a * d_kp;
          const double prop0 = ad0 * inv_a0;
          const double prop1 = (ad0 + d_kp) * inv_a1;
          const double c = c_row[static_cast<std::size_t>(p)];
          const double rp = r_row[static_cast<std::size_t>(p)];
          if (code == OBS_Y_ZERO) {
            logp0 += count_loglik_y0_fast(prop0, rp, c, model_type);
            logp1 += count_loglik_y0_fast(prop1, rp, c, model_type);
          } else {
            const double yp = y_row[static_cast<std::size_t>(p)];
            logp0 += loglik_value_fast(prop0, yp, rp, c, model_type);
            logp1 += loglik_value_fast(prop1, yp, rp, c, model_type);
          }
        }
      } else {
        const double inv_a0 = 1.0 / a0;
        const double inv_a1 = 1.0 / (a0 + 1.0);
        for (int p = 0; p < P; ++p) {
          const double d_kp = D_ptr[k_idx + p * K];
          const double ad0 = full_ptr[p] - old_a * d_kp;
          const double prop0 = ad0 * inv_a0;
          const double prop1 = (ad0 + d_kp) * inv_a1;
          const double yp = y_row[static_cast<std::size_t>(p)];
          const double rp = r_row[static_cast<std::size_t>(p)];
          if (obs.has_loglik_const()) {
            const double c = obs.loglik_const_at(host, p);
            logp0 += loglik_value_fast(prop0, yp, rp, c, model.model_type);
            logp1 += loglik_value_fast(prop1, yp, rp, c, model.model_type);
          } else {
            logp0 += loglik_value(yp, rp, prop0, model.model_type);
            logp1 += loglik_value(yp, rp, prop1, model.model_type);
          }
        }
      }

      if (k == kstar && old_a == 1.0 &&
          col_tracker.sums[static_cast<std::size_t>(k_idx)] - old_a == 0.0) {
        logp1 -= log_mu[static_cast<std::size_t>(k_idx)];
        const int next_kstar = col_tracker.second_last_active_col();
        if (next_kstar >= 0) {
          logp0 -= log_mu[static_cast<std::size_t>(next_kstar)];
        }
      } else if (k > kstar) {
        logp1 -= log_mu[static_cast<std::size_t>(k_idx)];
        if (kstar >= 1) {
          logp0 -= log_mu[static_cast<std::size_t>(kstar - 1)];
        }
      }

      const double accept = mh_ratio(logp1, logp0);
      const double draw = unif_rand();
      A(host, k_idx) = draw < accept ? 1.0 : 0.0;

      if (A(host, k_idx) != old_a) {
        const double delta = A(host, k_idx) - old_a;
        full_row += delta * D_e.row(k_idx);
        row_sum += delta;
        col_tracker.apply_delta(k_idx, delta);
        if (A(host, k_idx) == 1.0) {
          kstar = std::max(kstar, k);
        } else if (k == kstar && col_tracker.sums[static_cast<std::size_t>(k_idx)] == 0.0) {
          kstar = col_tracker.kstar_from_active();
        }
      }
    }
  }
}

void update_d(SliceState& state, ModelData& model) {
  Rcpp::NumericMatrix& A = state.A;
  Rcpp::NumericMatrix& D = state.D;
  const int kmin = state.kmin;
  const int kstar = state.kstar;
  const int N = A.nrow();
  const int K = A.ncol();
  const int P = D.ncol();
  const ObsViews& obs = model.obs;

  std::vector<int> an(static_cast<std::size_t>(N));
  for (int i = 0; i < N; ++i) {
    double row_total = 0.0;
    for (int k = 0; k < K; ++k) {
      row_total += A(i, k);
    }
    an[static_cast<std::size_t>(i)] = static_cast<int>(row_total);
  }

  Eigen::Map<Eigen::MatrixXd> A_e(A.begin(), N, K);
  const int K_d = D.nrow();
  if (K != K_d) {
    Rcpp::stop("update_d: ncol(A)=%d must equal nrow(D)=%d", K, K_d);
  }
  Eigen::Map<Eigen::MatrixXd> D_e(D.begin(), K, P);
  Eigen::MatrixXd AD = A_e * D_e;

  std::vector<std::vector<int>> active_hosts;
  build_active_hosts_by_col(A_e, kmin, kstar, active_hosts);

  const double log_rho = std::log(model.rho);
  const double log_1m_rho = std::log(1.0 - model.rho);

  const int k_active_last = std::min(kstar, K);
  for (int k = kmin; k <= k_active_last; ++k) {
    const int k_idx = k - 1;
    if (k_idx < 0 || k_idx >= K) continue;
    const double* A_col = A_e.col(k_idx).data();
    const std::vector<int>& hosts = active_hosts[static_cast<std::size_t>(k_idx)];

    if (hosts.empty()) {
      for (int p = 0; p < P; ++p) {
        D_e(k_idx, p) = unif_rand() < model.rho ? 1.0 : 0.0;
      }
      continue;
    }

    for (int p = 0; p < P; ++p) {
      const double* ad_col = AD.col(p).data();
      const double d_kp = D_e(k_idx, p);

      const double logp0 = log_1m_rho + loglik_target_col_active(
        hosts, ad_col, A_col, d_kp, an.data(), obs, p, model.model_type, false);
      const double logp1 = log_rho + loglik_target_col_active(
        hosts, ad_col, A_col, d_kp, an.data(), obs, p, model.model_type, true);

      const double old_d = d_kp;
      const double new_d = unif_rand() < mh_ratio(logp1, logp0) ? 1.0 : 0.0;
      D_e(k_idx, p) = new_d;
      if (new_d != old_d) {
        AD.col(p) += A_e.col(k_idx) * (new_d - old_d);
      }
    }
  }
}

double cell_loglik_at_prop(const ObsViews& obs,
                           int host,
                           int target,
                           double prop,
                           int model_type) {
  const int p = target;
  if (model_type == CATEGORICAL) {
    return loglik_categorical_value(
      obs.y_at(host, p), prop, obs.llik_tab, obs.llik_tab_nrow);
  }
  if (obs.has_obs_code()) {
    const int code = obs.obs_code_at(host, p);
    if (code == OBS_SKIP) return 0.0;
    if (obs.has_loglik_const()) {
      const double c = obs.loglik_const_at(host, p);
      const double r = obs.r_at(host, p);
      if (code == OBS_Y_ZERO) {
        return count_loglik_y0_fast(prop, r, c, model_type);
      }
      return loglik_value_fast(prop, obs.y_at(host, p), r, c, model_type);
    }
  }
  const double y = obs.y_at(host, p);
  const double r = obs.r_at(host, p);
  if (obs.has_loglik_const()) {
    return loglik_value_fast(prop, y, r, obs.loglik_const_at(host, p), model_type);
  }
  return loglik_value(y, r, prop, model_type);
}

double total_loglik_matrix(const Eigen::MatrixXd& A_e,
                           const Eigen::MatrixXd& D_e,
                           const ObsViews& obs,
                           int model_type) {
  const int N = obs.N;
  const int P = obs.P;
  const Eigen::MatrixXd ad = A_e * D_e;
  const Eigen::VectorXd an = A_e.rowwise().sum();
  double ll = 0.0;
  for (int p = 0; p < P; ++p) {
    for (int i = 0; i < N; ++i) {
      const double an_i = an(i);
      const double prop = (an_i > 0.0) ? ad(i, p) / an_i : R_NaN;
      ll += cell_loglik_at_prop(obs, i, p, prop, model_type);
    }
  }
  return ll;
}

}  // namespace kernel
}  // namespace snp_slicer

// [[Rcpp::export]]
double cpp_loglik_total(Rcpp::NumericMatrix A,
                        Rcpp::NumericMatrix D,
                        Rcpp::NumericMatrix y,
                        Rcpp::NumericMatrix r,
                        int model_type,
                        Rcpp::NumericMatrix llik_tab,
                        Rcpp::Nullable<Rcpp::NumericVector> loglik_const = R_NilValue,
                        Rcpp::Nullable<Rcpp::IntegerVector> obs_code = R_NilValue) {
  snp_slicer::kernel::ModelData model;
  model.y = y;
  model.r = r;
  model.llik_tab = llik_tab;
  model.model_type = model_type;
  snp_slicer::kernel::prepare_obs_views(model, loglik_const, obs_code);

  const int N = A.nrow();
  const int K = A.ncol();
  const int P = D.ncol();
  Eigen::Map<const Eigen::MatrixXd> A_e(A.begin(), N, K);
  Eigen::Map<const Eigen::MatrixXd> D_e(D.begin(), K, P);
  return snp_slicer::kernel::total_loglik_matrix(A_e, D_e, model.obs, model_type);
}

// [[Rcpp::export]]
Rcpp::List cpp_build_kernel_obs_cache(Rcpp::NumericMatrix y,
                                      Rcpp::NumericMatrix r,
                                      int model_type) {
  snp_slicer::kernel::ModelData model;
  model.y = y;
  model.r = r;
  model.llik_tab = Rcpp::NumericMatrix(3, 3);
  std::fill(model.llik_tab.begin(), model.llik_tab.end(), 0.0);
  model.model_type = model_type;
  model.obs.bind(model.y, model.r, model.llik_tab);

  if (model_type == snp_slicer::CATEGORICAL) {
    snp_slicer::kernel::fill_obs_code_buffer(model);
    return Rcpp::List::create(
      Rcpp::_["loglik_const"] = Rcpp::NumericVector(0),
      Rcpp::_["obs_code"] = model.obs_code_vec
    );
  }

  snp_slicer::kernel::fill_loglik_const_buffer(model);
  snp_slicer::kernel::fill_obs_code_buffer(model);

  return Rcpp::List::create(
    Rcpp::_["loglik_const"] = Rcpp::NumericVector(
      model.loglik_const_storage.begin(), model.loglik_const_storage.end()),
    Rcpp::_["obs_code"] = model.obs_code_vec
  );
}

// [[Rcpp::export]]
Rcpp::List cpp_update_a(Rcpp::NumericMatrix A,
                        Rcpp::NumericMatrix D,
                        Rcpp::NumericVector mu,
                        Rcpp::IntegerVector mixed,
                        int kplus,
                        int kstar,
                        Rcpp::NumericMatrix y,
                        Rcpp::NumericMatrix r,
                        int model_type,
                        Rcpp::NumericMatrix llik_tab,
                        Rcpp::Nullable<Rcpp::NumericVector> loglik_const = R_NilValue,
                        Rcpp::Nullable<Rcpp::IntegerVector> obs_code = R_NilValue) {
  snp_slicer::kernel::SliceState state;
  state.A = A;
  state.D = D;
  state.mu = mu;
  state.kplus = kplus;
  state.kstar = kstar;

  snp_slicer::kernel::ModelData model;
  model.y = y;
  model.r = r;
  model.llik_tab = llik_tab;
  model.mixed = mixed;
  model.model_type = model_type;
  snp_slicer::kernel::prepare_obs_views(model, loglik_const, obs_code);

  snp_slicer::kernel::update_a(state, model);

  return Rcpp::List::create(
    Rcpp::_["A"] = state.A,
    Rcpp::_["kstar"] = state.kstar
  );
}

// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_update_d(Rcpp::NumericMatrix A,
                                 Rcpp::NumericMatrix D,
                                 Rcpp::IntegerVector an,
                                 int kmin,
                                 int kstar,
                                 Rcpp::NumericMatrix y,
                                 Rcpp::NumericMatrix r,
                                 double rho,
                                 int model_type,
                                 Rcpp::NumericMatrix llik_tab,
                                 Rcpp::Nullable<Rcpp::NumericVector> loglik_const = R_NilValue,
                                 Rcpp::Nullable<Rcpp::IntegerVector> obs_code = R_NilValue) {
  snp_slicer::kernel::SliceState state;
  state.A = A;
  state.D = D;
  state.kmin = kmin;
  state.kstar = kstar;

  snp_slicer::kernel::ModelData model;
  model.y = y;
  model.r = r;
  model.llik_tab = llik_tab;
  model.rho = rho;
  model.model_type = model_type;
  snp_slicer::kernel::prepare_obs_views(model, loglik_const, obs_code);

  snp_slicer::kernel::update_d(state, model);
  return state.D;
}

// [[Rcpp::export]]
Rcpp::List cpp_slice_iter(Rcpp::NumericMatrix A,
                          Rcpp::NumericMatrix D,
                          Rcpp::NumericVector mu,
                          Rcpp::IntegerVector mixed,
                          int kplus,
                          int kstar,
                          int kmin,
                          int ktrunc,
                          Rcpp::NumericMatrix y,
                          Rcpp::NumericMatrix r,
                          double rho,
                          double alpha,
                          int N,
                          int P,
                          int model_type,
                          Rcpp::NumericMatrix llik_tab,
                          Rcpp::Nullable<Rcpp::NumericVector> loglik_const = R_NilValue,
                          Rcpp::Nullable<Rcpp::IntegerVector> obs_code = R_NilValue) {
  Rcpp::RNGScope rng_scope;
  snp_slicer::kernel::SliceState state;
  state.A = A;
  state.D = D;
  state.mu = mu;
  state.kplus = kplus;
  state.kstar = kstar;
  state.kmin = kmin;
  state.ktrunc = ktrunc;

  snp_slicer::kernel::ModelData model;
  model.y = y;
  model.r = r;
  model.llik_tab = llik_tab;
  model.mixed = mixed;
  model.rho = rho;
  model.alpha = alpha;
  model.N = N;
  model.P = P;
  model.model_type = model_type;
  snp_slicer::kernel::prepare_obs_views(model, loglik_const, obs_code);

  snp_slicer::kernel::update_s(state, model);
  snp_slicer::kernel::update_a(state, model);
  snp_slicer::kernel::update_d(state, model);
  snp_slicer::kernel::update_mu(state, model);

  return Rcpp::List::create(
    Rcpp::_["A"] = state.A,
    Rcpp::_["D"] = state.D,
    Rcpp::_["mu"] = state.mu,
    Rcpp::_["kplus"] = state.kplus,
    Rcpp::_["kstar"] = state.kstar,
    Rcpp::_["ktrunc"] = state.ktrunc
  );
}
