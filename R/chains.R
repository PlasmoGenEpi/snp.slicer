#' Draw independent per-chain seeds
#'
#' Generate a separate random seed for each chain, based on the 
#' caller's random number generator to maintain reproducibility.
#'
#' @param n_chains Number of chains
#'
#' @details Seeds are drawn from the caller's RNG stream, so setting \code{seed} in
#' \code{\link{snp_slice}} makes the whole multi-chain run reproducible while
#' each chain still starts from a different point. Because each chain receives an
#' explicit integer seed rather than inheriting global RNG state, the same
#' seeding scheme works for a compiled kernel with its own RNG.
#'
#' @return Integer vector of length \code{n_chains}
#' @keywords internal
generate_chain_seeds <- function(n_chains) {
  sample.int(.Machine$integer.max, n_chains, replace = FALSE)
}

#' Dispatch a per-chain function, possibly in parallel
#'
#' Run the provided function on each chain, in parallel if 
#' \code{n_cores} > 1.
#'
#' @param n_chains Number of chains
#' @param n_cores Number of cores to use
#' @param fun Function of a single argument (the chain index)
#'
#' @details Uses forking (\code{parallel::mclapply}) where available and a PSOCK cluster
#' on Windows. Falls back to \code{lapply} when a single core is requested.
#'
#' @return List of chain results
#' @keywords internal
dispatch_chains <- function(n_chains, n_cores, fun) {
  n_cores <- min(n_cores, n_chains)

  if (n_cores == 1) {
    return(lapply(seq_len(n_chains), fun))
  }

  if (.Platform$OS.type == "windows") {
    cl <- parallel::makePSOCKcluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, requireNamespace("snp.slicer", quietly = TRUE))
    return(parallel::parLapply(cl, seq_len(n_chains), fun))
  }

  parallel::mclapply(seq_len(n_chains), fun, mc.cores = n_cores)
}

#' Run one or more MCMC chains
#'
#' Orchestrates \code{\link{run_chain}}: draws per-chain seeds, dispatches the
#' chains (serially or in parallel), and collects them into a single result. The
#' chain that reaches the highest MAP log posterior is reported as the best
#' chain and drives the top-level estimates.
#'
#' @param ... Parameters passed on to \code{run_chain}
#' @param n_chains Number of independent chains to run
#' @param n_cores Number of cores used to run chains simultaneously
#' @param verbose Whether to print progress. Named explicitly because this
#'   function reports on the chain set as a whole as well as passing it down.
#'
#' @return List with \code{chains} (per-chain results) and \code{best_chain}.
#'   \code{create_results_object()} keeps those two and the shared
#'   \code{parameters}; nothing else is carried through to the results object.
#' @keywords internal
run_chains <- function(..., n_chains = 1L, n_cores = 1L, verbose = TRUE) {

  seeds <- generate_chain_seeds(n_chains)

  chain_fun <- function(i) {
    run_chain(
      ...,
      verbose = verbose,
      chain_id = i,
      seed = seeds[i]
    )
  }

  if (n_chains == 1L) {
    chains <- list(chain_fun(1L))
  } else {
    if (verbose) {
      cat("Running", n_chains, "chains on", min(n_cores, n_chains), "core(s)\n")
    }
    chains <- dispatch_chains(n_chains, n_cores, chain_fun)
  }

  # mclapply reports worker errors as try-error elements rather than signalling
  failed <- vapply(chains, function(x) inherits(x, "try-error") || is.null(x$map_state),
                   logical(1))
  if (any(failed)) {
    first <- which(failed)[1]
    stop("Chain ", first, " failed: ",
         if (inherits(chains[[first]], "try-error")) {
           conditionMessage(attr(chains[[first]], "condition"))
         } else {
           "chain returned no MAP state"
         })
  }

  map_logpost <- vapply(chains, function(x) x$diagnostics$map_logpost, numeric(1))
  best_chain <- which.max(map_logpost)

  if (verbose && n_chains > 1L) {
    cat("Best chain:", best_chain, "(MAP log posterior",
        round(map_logpost[best_chain], 2), ")\n")
  }

  results <- chains[[best_chain]]
  results$chains <- chains
  results$best_chain <- best_chain
  results$parameters$n_chains <- n_chains
  results$parameters$n_cores <- n_cores

  results
}
