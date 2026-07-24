# Invariant scalar summaries used for convergence diagnostics.
#
# Strain labels are not identifiable across iterations or chains (column k of A
# is not a stable quantity), so R-hat / ESS are only meaningful on
# permutation-invariant summaries of each stored sample. Each extractor maps one
# stored sample to a scalar, except `coi`, which returns the per-host complexity
# of infection (a length-N vector that expands to coi[1]..coi[N]).
invariant_extractors <- list(
  logpost   = function(s) s$logpost,
  n_strains = function(s) sum(colSums(s$A) > 0),
  kstar     = function(s) s$kstar,
  ktrunc    = function(s) s$ktrunc,
  coi       = function(s) rowSums(s$A)
)

#' Assemble a posterior draws array from a SNP-Slice run
#'
#' @description
#' Converts the retained MCMC samples of every chain into a
#' \code{posterior::draws_array} (dimensions iteration \eqn{\times} chain
#' \eqn{\times} variable) of permutation-invariant scalar summaries. This is the
#' reusable primitive behind \code{\link{convergence_diagnostics}}: because it
#' returns a standard \code{draws_array}, any function in the \pkg{posterior}
#' package can be applied to it directly.
#'
#' Diagnostics are computed on invariant summaries rather than the raw
#' \code{A}/\code{D} matrices because strain labels are not identifiable across
#' iterations or chains, which makes element-wise R-hat / ESS on \code{A} or
#' \code{D} uninterpretable.
#'
#' @param results A \code{snp_slice_results} object (with \code{store_mcmc =
#'   TRUE}).
#' @param pars Character vector of parameters to include. One of \code{"logpost"},
#'   \code{"n_strains"}, \code{"kstar"}, \code{"ktrunc"}, which are all scalar per sample, or
#'   \code{"coi"} (per-host complexity of infection), which is a vector that expands to
#'   \code{coi[1]..coi[N]}).
#' @param additional_burnin Number of additional stored samples to discard from
#'   the start of each retained chain before assembling the array.
#'
#' @details
#' All chains are pooled into one array so that R-hat and ESS are computed for the
#' whole run. Without early stopping every chain retains the same number of
#' samples; when \code{gap} early-stopping fired, chains can differ in length, in
#' which case every chain is truncated to the shortest and a warning is issued.
#'
#' @return A \code{posterior::draws_array}.
#' @export
#' @examples
#' result <- load_example_results()
#' draws <- as_draws_snp_slice(result)
#' posterior::summarise_draws(draws)
as_draws_snp_slice <- function(results,
                               pars = c("logpost", "n_strains", "kstar", "ktrunc"),
                               additional_burnin = 0) {
  if (!inherits(results, "snp_slice_results")) {
    stop("Input must be an snp_slice_results object")
  }

  valid_pars <- names(invariant_extractors)
  unknown <- setdiff(pars, valid_pars)
  if (length(unknown) > 0) {
    stop("Unknown pars: ", paste(unknown, collapse = ", "),
         ". Valid values are: ", paste(valid_pars, collapse = ", "))
  }
  if (length(pars) == 0) {
    stop("pars must name at least one parameter")
  }

  # One list of retained samples per chain, reusing the existing accessors.
  n_chains <- if (is.null(results$chains)) 1L else length(results$chains)
  chain_samples <- lapply(seq_len(n_chains), function(i) {
    retained_samples(get_chain(results, i), additional_burnin)
  })

  # Chains match in length unless gap early-stopping truncated some of them.
  lengths <- vapply(chain_samples, length, integer(1))
  n_iter <- min(lengths)
  if (length(unique(lengths)) > 1L) {
    warning("Chains have unequal retained lengths (shortest = ", n_iter,
            "); this usually means gap early-stopping fired at different ",
            "iterations. Truncating every chain to the first ", n_iter,
            " samples. Rerun without gap for equal-length chains.")
    chain_samples <- lapply(chain_samples, function(s) s[seq_len(n_iter)])
  }

  # Per-chain matrix of invariant summaries: n_iter rows x n_var columns, with
  # scalar pars contributing one column and coi contributing one per host.
  chain_matrices <- lapply(chain_samples, sample_variable_matrix, pars = pars)

  var_names <- colnames(chain_matrices[[1]])
  draws <- array(
    NA_real_,
    dim = c(n_iter, n_chains, length(var_names)),
    dimnames = list(iteration = NULL, chain = NULL, variable = var_names)
  )
  for (i in seq_len(n_chains)) {
    draws[, i, ] <- chain_matrices[[i]]
  }

  posterior::as_draws_array(draws)
}

#' Build the invariant-summary matrix for one chain
#'
#' @param samples Retained samples for one chain
#' @param pars Parameters to extract
#'
#' @return Numeric matrix, one row per sample and one named column per (expanded)
#'   variable
#' @keywords internal
sample_variable_matrix <- function(samples, pars) {
  columns <- lapply(pars, function(par) {
    values <- lapply(samples, invariant_extractors[[par]])
    mat <- do.call(rbind, values)
    if (ncol(mat) == 1L) {
      colnames(mat) <- par
    } else {
      colnames(mat) <- paste0(par, "[", seq_len(ncol(mat)), "]")
    }
    mat
  })
  do.call(cbind, columns)
}

#' Convergence diagnostics (R-hat and ESS) for a SNP-Slice run
#'
#' @description
#' Reports rank-normalized split-R-hat and effective sample size (bulk and tail),
#' pooled across all chains, for permutation-invariant summaries of the run. One
#' row is returned per parameter for the whole run rather than one per chain, so
#' R-hat reflects between-chain agreement and ESS reflects the total effective
#' sample count. A single-chain run still yields a (split-)R-hat.
#'
#' @inheritParams as_draws_snp_slice
#'
#' @details
#' Constant parameters (e.g. \code{kstar} or \code{ktrunc} are frequently constant
#' across retained samples) have zero variance, for which \pkg{posterior} returns
#' \code{NaN} R-hat and ESS. That means "this parameter did not move", not that
#' the diagnostic failed. To diagnose the \code{A}/\code{D} matrices, use
#' \code{"coi"}, which is invariant to strain relabeling; element-wise R-hat / ESS
#' on \code{A}/\code{D} is intentionally not provided (see
#' \code{\link{as_draws_snp_slice}}).
#'
#' @return A data frame with one row per variable and columns \code{variable},
#'   \code{mean}, \code{median}, \code{sd}, \code{q5}, \code{q95}, \code{rhat},
#'   \code{ess_bulk}, \code{ess_tail}.
#' @export
#' @examples
#' result <- load_example_results()
#' convergence_diagnostics(result)
#' # Per-host complexity of infection
#' convergence_diagnostics(result, pars = "coi")
convergence_diagnostics <- function(results,
                                    pars = c("logpost", "n_strains", "kstar", "ktrunc"),
                                    additional_burnin = 0) {
  draws <- as_draws_snp_slice(results, pars = pars, additional_burnin = additional_burnin)
  summary <- posterior::summarise_draws(draws)
  keep <- c("variable", "mean", "median", "sd", "q5", "q95",
            "rhat", "ess_bulk", "ess_tail")
  as.data.frame(summary[, keep])
}
