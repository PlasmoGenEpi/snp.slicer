#ifndef SNP_SLICER_MCMC_KERNEL_INTERNAL_H
#define SNP_SLICER_MCMC_KERNEL_INTERNAL_H

#include "mcmc_kernel_common.h"
#include <Rcpp.h>
#include <vector>

namespace snp_slicer {
namespace kernel {

// obs_code layout (length N*P, column-major): 0 = NA/skip, 1 = y==0, 2 = y>0
constexpr int OBS_SKIP = 0;
constexpr int OBS_Y_ZERO = 1;
constexpr int OBS_Y_POS = 2;

struct ObsViews {
  int N = 0;
  int P = 0;
  const double* y_col = nullptr;
  const double* r_col = nullptr;
  const double* llik_tab = nullptr;
  const double* loglik_const = nullptr;
  const int* obs_code = nullptr;
  int llik_tab_nrow = 0;

  void bind(const Rcpp::NumericMatrix& y,
            const Rcpp::NumericMatrix& r,
            const Rcpp::NumericMatrix& tab) {
    N = y.nrow();
    P = y.ncol();
    y_col = y.begin();
    r_col = r.begin();
    llik_tab = tab.begin();
    llik_tab_nrow = tab.nrow();
    loglik_const = nullptr;
    obs_code = nullptr;
  }

  double y_at(int i, int p) const { return y_col[i + p * N]; }
  double r_at(int i, int p) const { return r_col[i + p * N]; }
  int obs_code_at(int i, int p) const {
    if (obs_code == nullptr) return OBS_Y_POS;
    return obs_code[i + p * N];
  }
  double loglik_const_at(int i, int p) const {
    if (loglik_const == nullptr) return 0.0;
    return loglik_const[i + p * N];
  }
  bool has_loglik_const() const { return loglik_const != nullptr; }
  bool has_obs_code() const { return obs_code != nullptr; }
};

struct SliceState {
  Rcpp::NumericMatrix A;
  Rcpp::NumericMatrix D;
  Rcpp::NumericVector mu;
  int kplus;
  int kstar;
  int ktrunc;
  int kmin;
};

struct ModelData {
  Rcpp::NumericMatrix y;
  Rcpp::NumericMatrix r;
  Rcpp::NumericMatrix llik_tab;
  Rcpp::IntegerVector mixed;
  ObsViews obs;
  Rcpp::NumericVector loglik_const_vec;
  Rcpp::IntegerVector obs_code_vec;
  std::vector<double> loglik_const_storage;
  double rho;
  double alpha;
  int N;
  int P;
  int model_type;
};

inline void fill_loglik_const_buffer(ModelData& model) {
  const int N = model.obs.N;
  const int P = model.obs.P;
  const std::size_t n_cells =
    static_cast<std::size_t>(N) * static_cast<std::size_t>(P);
  model.loglik_const_storage.resize(n_cells);
  double* const out = model.loglik_const_storage.data();
  for (int p = 0; p < P; ++p) {
    for (int i = 0; i < N; ++i) {
      const double y = model.obs.y_at(i, p);
      const double r = model.obs.r_at(i, p);
      const std::size_t idx = static_cast<std::size_t>(i) +
                              static_cast<std::size_t>(p) * static_cast<std::size_t>(N);
      switch (model.model_type) {
      case POISSON:
        out[idx] = dpois_log_const(y);
        break;
      case BINOMIAL:
        out[idx] = dbinom_log_const(y, r);
        break;
      case NEGATIVE_BINOMIAL:
        out[idx] = dnbinom_log_const(y, r);
        break;
      default:
        out[idx] = 0.0;
        break;
      }
    }
  }
  model.obs.loglik_const = model.loglik_const_storage.data();
}

inline void fill_obs_code_buffer(ModelData& model) {
  const int N = model.obs.N;
  const int P = model.obs.P;
  model.obs_code_vec = Rcpp::IntegerVector(N * P);
  int* const out = model.obs_code_vec.begin();
  for (int p = 0; p < P; ++p) {
    for (int i = 0; i < N; ++i) {
      const double y = model.obs.y_at(i, p);
      const std::size_t idx = static_cast<std::size_t>(i) +
                              static_cast<std::size_t>(p) * static_cast<std::size_t>(N);
      if (ISNA(y)) {
        out[idx] = OBS_SKIP;
      } else if (y == 0.0) {
        out[idx] = OBS_Y_ZERO;
      } else {
        out[idx] = OBS_Y_POS;
      }
    }
  }
  model.obs.obs_code = model.obs_code_vec.begin();
}

inline void prepare_obs_views(
    ModelData& model,
    Rcpp::Nullable<Rcpp::NumericVector> loglik_const_in = R_NilValue,
    Rcpp::Nullable<Rcpp::IntegerVector> obs_code_in = R_NilValue) {
  model.obs.bind(model.y, model.r, model.llik_tab);
  const int N = model.obs.N;
  const int P = model.obs.P;
  model.loglik_const_storage.clear();
  model.obs.loglik_const = nullptr;
  model.obs.obs_code = nullptr;
  model.loglik_const_vec = Rcpp::NumericVector();
  model.obs_code_vec = Rcpp::IntegerVector();

  const std::size_t n_cells =
    static_cast<std::size_t>(N) * static_cast<std::size_t>(P);

  if (obs_code_in.isNotNull()) {
    model.obs_code_vec = obs_code_in.get();
    if (static_cast<std::size_t>(model.obs_code_vec.size()) == n_cells) {
      model.obs.obs_code = model.obs_code_vec.begin();
    }
  }

  if (model.model_type == CATEGORICAL || N <= 0 || P <= 0) {
    if (model.obs.obs_code == nullptr && N > 0 && P > 0) {
      fill_obs_code_buffer(model);
    }
    return;
  }

  if (loglik_const_in.isNotNull()) {
    model.loglik_const_vec = loglik_const_in.get();
    if (static_cast<std::size_t>(model.loglik_const_vec.size()) == n_cells) {
      model.obs.loglik_const = model.loglik_const_vec.begin();
    }
  }
  if (model.obs.loglik_const == nullptr) {
    fill_loglik_const_buffer(model);
  }
  if (model.obs.obs_code == nullptr) {
    fill_obs_code_buffer(model);
  }
}

void update_s(SliceState& state, const ModelData& model);
void update_a(SliceState& state, ModelData& model);
void update_d(SliceState& state, ModelData& model);
void update_mu(SliceState& state, const ModelData& model);

}  // namespace kernel
}  // namespace snp_slicer

#endif
