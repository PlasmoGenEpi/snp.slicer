#' MCMC kernel seam
#'
#' The kernel module owns the hot slice-sampler loops. Callers cross this seam
#' via \code{model_obj$kernel}; tests can swap the R reference adapter for the
#' compiled adapter without changing \code{slice_iter()}.
#'
#' For count/categorical models with the compiled backend, log-likelihood layout
#' caches (\code{kernel_loglik_const}, \code{kernel_obs_code}) are built
#' automatically on first use (\code{resolve_mcmc_kernel}, \code{slice_iter}).
#'
#' @useDynLib snp.slicer, .registration = TRUE
#' @import Rcpp
#' @name mcmc_kernel
#' @keywords internal
NULL

NULL_LLIK_TAB <- matrix(0, nrow = 3, ncol = 3)

#' @return \code{TRUE} when compiled kernel entry points are loadable.
#' @keywords internal
cpp_kernel_available <- function() {
  exists("cpp_slice_iter", where = asNamespace("snp.slicer"), inherits = FALSE) &&
    exists("cpp_update_a", where = asNamespace("snp.slicer"), inherits = FALSE)
}

#' Map observation model name to compiled kernel type id
#' @keywords internal
model_type_id <- function(model_name) {
  switch(model_name,
    poisson = 0L,
    binomial = 1L,
    negative_binomial = 2L,
    categorical = 3L,
    stop("Model not supported by compiled MCMC kernel: ", model_name)
  )
}

#' Lookup table passed to compiled categorical updates
#' @keywords internal
kernel_llik_tab <- function(model_obj) {
  if (model_obj$name == "categorical") {
    model_obj$llik_tab
  } else {
    NULL_LLIK_TAB
  }
}

#' Whether a model can use the compiled kernel obs cache
#' @keywords internal
kernel_model_uses_cpp_cache <- function(model_obj) {
  model_obj$name %in% c(
    "poisson", "binomial", "negative_binomial", "categorical"
  ) && cpp_kernel_available()
}

#' Build cached loglik constants and per-cell obs layout for the C++ kernel
#' @keywords internal
build_kernel_obs_cache <- function(model_obj) {
  if (!kernel_model_uses_cpp_cache(model_obj)) {
    return(list(loglik_const = NULL, obs_code = NULL))
  }
  cpp_build_kernel_obs_cache(
    y = model_obj$y,
    r = model_obj$r,
    model_type = model_type_id(model_obj$name)
  )
}

#' Whether precomputed kernel obs cache is attached to a model
#' @keywords internal
kernel_obs_cache_ready <- function(model_obj) {
  n <- model_obj$N * model_obj$P
  if (n <= 0L) {
    return(TRUE)
  }
  obs_ok <- !is.null(model_obj$kernel_obs_code) &&
    length(model_obj$kernel_obs_code) == n
  if (!obs_ok) {
    return(FALSE)
  }
  if (model_obj$name == "categorical") {
    return(TRUE)
  }
  !is.null(model_obj$kernel_loglik_const) &&
    length(model_obj$kernel_loglik_const) == n
}

#' Ensure obs cache vectors exist on \code{model_obj} (idempotent, in-place safe)
#'
#' Called automatically before compiled kernel updates; callers do not need to
#' invoke this explicitly.
#' @keywords internal
ensure_kernel_obs_cache <- function(model_obj) {
  if (!kernel_model_uses_cpp_cache(model_obj) || kernel_obs_cache_ready(model_obj)) {
    return(model_obj)
  }
  cache <- build_kernel_obs_cache(model_obj)
  model_obj$kernel_loglik_const <- cache$loglik_const
  model_obj$kernel_obs_code <- cache$obs_code
  model_obj
}

#' @rdname ensure_kernel_obs_cache
#' @keywords internal
attach_kernel_obs_cache <- ensure_kernel_obs_cache

#' Arguments for cached obs layout passed into cpp kernel entry points
#' @keywords internal
kernel_obs_args <- function(model_obj) {
  model_obj <- ensure_kernel_obs_cache(model_obj)
  list(
    loglik_const = model_obj$kernel_loglik_const,
    obs_code = model_obj$kernel_obs_code
  )
}

#' Obs cache for a compiled kernel call (kernel object or model)
#' @keywords internal
kernel_cpp_obs_args <- function(kernel, model_obj) {
  n <- model_obj$N * model_obj$P
  if (!is.null(kernel$obs_code) && length(kernel$obs_code) == n) {
    return(list(
      loglik_const = kernel$obs_loglik_const,
      obs_code = kernel$obs_code
    ))
  }
  kernel_obs_args(model_obj)
}

#' Attach obs cache vectors to a compiled kernel adapter list
#' @keywords internal
kernel_with_obs_cache <- function(kernel, model_obj) {
  model_obj <- ensure_kernel_obs_cache(model_obj)
  kernel$obs_loglik_const <- model_obj$kernel_loglik_const
  kernel$obs_code <- model_obj$kernel_obs_code
  kernel
}

#' Preserve dimnames when compiled code returns a fresh matrix
#' @keywords internal
restore_matrix_dimnames <- function(updated, template) {
  dn <- dimnames(template)
  if (is.null(dn)) {
    return(updated)
  }
  out <- dn
  if (!is.null(dn[[1]]) && length(dn[[1]]) != nrow(updated)) {
    out[[1]] <- NULL
  }
  if (!is.null(dn[[2]]) && length(dn[[2]]) != ncol(updated)) {
    out[[2]] <- NULL
  }
  dimnames(updated) <- out
  updated
}

#' Resolve kernel adapter and attach cpp obs cache on \code{model_obj}
#'
#' Call before a manual \code{slice_iter()} loop when not using \code{run_chain()}.
#' \code{run_chain()} and \code{snp_slice()} invoke this automatically.
#'
#' @param model_obj Model object from \code{create_model()}.
#' @param backend One of \code{"auto"}, \code{"r"}, or \code{"cpp"}.
#' @return Updated model object (same list, with \code{kernel} and cache fields set).
#' @keywords internal
setup_mcmc_kernel <- function(model_obj, backend = getOption("snp.slicer.mcmc_kernel", "auto")) {
  model_obj$kernel <- resolve_mcmc_kernel(model_obj, backend = backend)
  if (identical(model_obj$kernel$name, "cpp")) {
    model_obj <- ensure_kernel_obs_cache(model_obj)
    model_obj$kernel <- kernel_with_obs_cache(model_obj$kernel, model_obj)
  }
  model_obj
}

#' Resolve MCMC kernel adapter for a model
#'
#' @param model_obj Model object from \code{create_model()}.
#' @param backend One of \code{"auto"}, \code{"r"}, or \code{"cpp"}.
#' @return Kernel adapter list.
#' @keywords internal
resolve_mcmc_kernel <- function(model_obj, backend = getOption("snp.slicer.mcmc_kernel", "auto")) {
  supported <- model_obj$name %in% c(
    "poisson", "binomial", "negative_binomial", "categorical"
  )

  use_cpp <- switch(backend,
    cpp = supported,
    r = FALSE,
    auto = supported && cpp_kernel_available(),
    stop("Unknown MCMC kernel backend: ", backend)
  )

  if (use_cpp) {
    kernel_with_obs_cache(mcmc_kernel_cpp(), model_obj)
  } else {
    mcmc_kernel_r()
  }
}

#' R reference kernel adapter
#' @keywords internal
mcmc_kernel_r <- function() {
  list(
    name = "r",
    update_s = slice_update_s_r,
    update_mu = slice_update_mu_r,
    update_a = slice_update_a_r,
    update_d = slice_update_d_r
  )
}

#' Compiled kernel adapter
#' @keywords internal
mcmc_kernel_cpp <- function() {
  if (!cpp_kernel_available()) {
    stop("Compiled MCMC kernel is not available; rebuild the package with Rcpp.")
  }
  list(
    name = "cpp",
    update_iter = kernel_update_iter_cpp,
    update_s = kernel_update_s_cpp,
    update_mu = kernel_update_mu_cpp,
    update_a = kernel_update_a_cpp,
    update_d = kernel_update_d_cpp
  )
}

#' Fused compiled slice iteration (s, a, d, mu in one .Call)
#' @keywords internal
kernel_update_iter_cpp <- function(state, model_obj) {
  obs <- kernel_cpp_obs_args(model_obj$kernel, model_obj)
  result <- cpp_slice_iter(
    A = state$A,
    D = state$D,
    mu = state$mu,
    mixed = as.integer(state$mixed),
    kplus = as.integer(state$kplus),
    kstar = as.integer(state$kstar),
    kmin = as.integer(state$kmin),
    ktrunc = as.integer(state$ktrunc),
    y = model_obj$y,
    r = model_obj$r,
    rho = model_obj$rho,
    alpha = model_obj$alpha,
    N = as.integer(model_obj$N),
    P = as.integer(model_obj$P),
    model_type = model_type_id(model_obj$name),
    llik_tab = kernel_llik_tab(model_obj),
    loglik_const = obs$loglik_const,
    obs_code = obs$obs_code
  )
  state$A <- result$A
  state$D <- restore_matrix_dimnames(result$D, state$D)
  state$mu <- result$mu
  state$kplus <- result$kplus
  state$kstar <- result$kstar
  state$ktrunc <- result$ktrunc
  state
}

#' @keywords internal
kernel_update_s_cpp <- function(state, model_obj) {
  result <- cpp_update_s(
    A = state$A,
    D = state$D,
    mu = state$mu,
    ktrunc = as.integer(state$ktrunc),
    alpha = model_obj$alpha,
    rho = model_obj$rho,
    P = as.integer(model_obj$P)
  )
  state$A <- result$A
  state$D <- restore_matrix_dimnames(result$D, state$D)
  state$mu <- result$mu
  state$ktrunc <- result$ktrunc
  state$kplus <- result$kplus
  state
}

#' @keywords internal
kernel_update_mu_cpp <- function(state, model_obj) {
  state$mu <- cpp_update_mu(
    A = state$A,
    mu = state$mu,
    kplus = as.integer(state$kplus),
    ktrunc = as.integer(state$ktrunc),
    N = as.integer(model_obj$N),
    alpha = model_obj$alpha
  )
  state
}

#' @keywords internal
kernel_update_a_cpp <- function(state, model_obj) {
  obs <- kernel_cpp_obs_args(model_obj$kernel, model_obj)
  result <- cpp_update_a(
    A = state$A,
    D = state$D,
    mu = state$mu,
    mixed = as.integer(state$mixed),
    kplus = as.integer(state$kplus),
    kstar = as.integer(state$kstar),
    y = model_obj$y,
    r = model_obj$r,
    model_type = model_type_id(model_obj$name),
    llik_tab = kernel_llik_tab(model_obj),
    loglik_const = obs$loglik_const,
    obs_code = obs$obs_code
  )
  state$A <- result$A
  state$kstar <- result$kstar
  state
}

#' @keywords internal
kernel_update_d_cpp <- function(state, model_obj) {
  obs <- kernel_cpp_obs_args(model_obj$kernel, model_obj)
  state$D <- cpp_update_d(
    A = state$A,
    D = state$D,
    an = as.integer(rowSums(state$A)),
    kmin = as.integer(state$kmin),
    kstar = as.integer(state$kstar),
    y = model_obj$y,
    r = model_obj$r,
    rho = model_obj$rho,
    model_type = model_type_id(model_obj$name),
    llik_tab = kernel_llik_tab(model_obj),
    loglik_const = obs$loglik_const,
    obs_code = obs$obs_code
  )
  state
}
