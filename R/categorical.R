#' Categorical Model for SNP-Slice
#'
#' @description
#' Implementation of the categorical observation model for SNP-Slice.
#' This model handles categorical observations (0, 0.5, 1) with error parameters.

#' Build categorical log-likelihood lookup table
#'
#' @param e1 False-positive error rate
#' @param e2 False-negative / mixture error rate
#' @return 3x3 matrix: rows = proportion bins, cols = observation bins
#' @keywords internal
build_categorical_llik_tab <- function(e1, e2) {
  logf0 <- log(c(1 - e1, e1, 0))
  logf1 <- log(c(0, e1, 1 - e1))
  logfmix <- log(c(e2 / 2, 1 - e2, e2 / 2))
  rbind(logf0, logfmix, logf1)
}

#' Map proportion to lookup-table row (1-based)
#' @keywords internal
categorical_prop_bin <- function(prop) {
  if (prop <= 0) {
    1L
  } else if (prop <= 0.99) {
    2L
  } else {
    3L
  }
}

#' Vectorised form of [categorical_prop_bin()]
#'
#' The exception detector in [categorical_resolve_exceptions()] must classify
#' proportions exactly as the likelihood does, so both call into the same bin
#' edges. Keeping a separate scalar and matrix form is deliberate: the scalar
#' one is called per-cell in the likelihood's inner loop.
#'
#' @param prop Numeric vector or matrix of proportions.
#' @return Integer object of the same shape, with values 1, 2 or 3.
#' @keywords internal
categorical_prop_bin_vec <- function(prop) {
  out <- ifelse(prop <= 0, 1L, ifelse(prop <= 0.99, 2L, 3L))
  if (is.matrix(prop)) matrix(out, nrow = nrow(prop), ncol = ncol(prop)) else out
}

#' Strains a host could plausibly be carrying, given its own calls
#'
#' A mixed host used to be initialized carrying every strain in the dictionary.
#' That is not a biological state -- on a 400-specimen panel it starts hosts at a
#' COI of ~289 -- and it manufactures likelihood exceptions that cannot be
#' repaired. `llik_tab` is `-Inf` at [prop bin 3, y == 0], bin 3 being
#' `prop > 0.99`; a host carrying all K strains at a locus where K-1 of them
#' carry the alternate sits at `(K-1)/K`, which is inside bin 3 once K > 100.
#' [categorical_resolve_exceptions()] repairs such a cell by adding a single
#' reference-carrying strain, which only reaches bin 2 when K <= 100, so above
#' that the repair is arithmetically incapable of clearing the exception and
#' initialization aborts at the iteration cap.
#'
#' Restricting a host to the strains that match its own homozygous calls fixes
#' this at the source: at every locus called 0 or 1 the host then carries only
#' strains with that allele, so the proportion is exactly 0 or 1 and lands in the
#' bin the observation agrees with. Heterozygous and missing loci are left
#' unconstrained, since any mixture is consistent with them.
#'
#' The result is never empty: the dictionary is built by de-duplicating the
#' hosts' own resolved genotypes, so a host's own strain always matches it at
#' every homozygous locus. `fallback` guards the degenerate case anyway.
#'
#' @param y_row Observed calls for one host (0, 1, 0.5 or NA per locus).
#' @param D Strain dictionary, strains in rows.
#' @param fallback Strain index to use if nothing matches.
#' @return Integer vector of strain indices.
#' @keywords internal
categorical_compatible_strains <- function(y_row, D, fallback) {
  observed <- which(!is.na(y_row) & (y_row == 0 | y_row == 1))
  if (length(observed) == 0L) {
    return(seq_len(nrow(D)))
  }
  keep <- which(apply(D[, observed, drop = FALSE], 1L,
                      function(d) all(d == y_row[observed])))
  if (length(keep) == 0L) fallback else keep
}

#' Map observation to lookup-table column (1-based)
#' @keywords internal
categorical_y_bin <- function(y) {
  if (y == 0) {
    1L
  } else if (y == 0.5) {
    2L
  } else {
    3L
  }
}

#' Log-likelihood for categorical model (matrix version)
#'
#' @param A Allocation matrix
#' @param D Dictionary matrix
#' @param model_obj Model object containing data
#'
#' @return Log-likelihood value
#' @keywords internal
categorical_loglikelihood_matrix <- function(A, D, model_obj) {
  start_timer("categorical_loglikelihood_matrix")

  proportions <- (A %*% D) / rowSums(A)
  llik_tab <- model_obj$llik_tab
  loglik <- 0
  for (p in seq_len(model_obj$P)) {
    loglik <- loglik + categorical_loglikelihood_vector(
      proportions[, p],
      model_obj$y[, p],
      llik_tab = llik_tab
    )
  }

  end_timer("categorical_loglikelihood_matrix")
  as.numeric(loglik)
}

#' Log-likelihood for categorical model (vector version)
#'
#' @param propvec Proportion vector
#' @param yvec Observed categorical values
#' @param rvec Not used for categorical model (kept for interface consistency)
#'
#' @return Log-likelihood value
#' @keywords internal
categorical_loglikelihood_vector <- function(propvec, yvec, rvec = NULL, llik_tab = NULL) {
  if (is.null(llik_tab)) {
    llik_tab <- build_categorical_llik_tab(0.05, 0.05)
  }

  loglik <- 0
  for (idx in seq_along(propvec)) {
    y_val <- yvec[idx]
    if (is.na(y_val)) {
      next
    }
    pb <- categorical_prop_bin(propvec[idx])
    yb <- categorical_y_bin(y_val)
    loglik <- loglik + llik_tab[pb, yb]
  }

  as.numeric(loglik)
}

#' Initialize state for categorical model
#'
#' @param model_obj Model object
#' @param threshold Threshold for identifying single infections
#'
#' @return Initialized state
#' @keywords internal
categorical_initialize_state <- function(model_obj, threshold = 0.001) {
  y <- model_obj$y
  N <- model_obj$N
  P <- model_obj$P
  
  # Replace NA or 0.5 in y with 0 or 1 (randomly)
  cate <- y
  for (i in 1:N) {
    for (j in 1:P) {
      if (is.na(y[i, j]) || y[i, j] == 0.5) {
        cate[i, j] <- (stats::runif(1) < 0.5)
      }
    }
  }
  
  # Remove duplicates to create dictionary
  ratios2dict <- remove_duplicates(cate)
  assignments <- ratios2dict$assignments
  nstrain <- length(unique(assignments))
  
  # Initialize state
  state <- list()
  state$D <- ratios2dict$D
  state$A <- matrix(0, nrow = N, ncol = nstrain)
  
  # Identify single vs mixed infections
  is_single <- apply(y, 1, function(x) all(x == 0 | x == 1))
  nsingleppl <- sum(is_single, na.rm = TRUE)
  
  if (nsingleppl >= 1) {
    # Handle single infections
    which_single <- which(is_single)
    which_mixed <- c(which(!is_single), which(is.na(is_single)))
    nsinglestrain <- length(unique(assignments[which_single]))
    
    # Assign mixed infections to the strains compatible with their own calls
    for (imixed in which_mixed) {
      state$A[imixed, categorical_compatible_strains(y[imixed, ], state$D,
                                                     assignments[imixed])] <- 1
    }

    # Assign single infections to specific strains
    for (isingle in 1:nsingleppl) {
      state$A[which_single[isingle], assignments[which_single[isingle]]] <- 1
    }
    
    # Reorder strains (single infection strains first)
    ord <- rep(NA, nstrain)
    ord[1:nsinglestrain] <- unique(assignments[which_single])
    remaining_strains <- setdiff(unique(assignments), unique(assignments[which_single]))
    if (length(remaining_strains) > 0) {
      ord[(nsinglestrain + 1):(nsinglestrain + length(remaining_strains))] <- remaining_strains
    }
    
    state$A <- state$A[, ord, drop = FALSE]
    state$D <- state$D[ord, , drop = FALSE]
    state$mixed <- which_mixed
    state$kmin <- nsinglestrain + 1
    
    # Resolve exceptions
    state$loglik <- categorical_loglikelihood_matrix(state$A, state$D, model_obj)
    iter <- 1
    while (is.infinite(state$loglik)) {
      state <- categorical_resolve_exceptions(state, model_obj)
      iter <- iter + 1
      if (iter > 10) {
        stop("Failed to initialize valid state after 10 iterations")
      }
      state$loglik <- categorical_loglikelihood_matrix(state$A, state$D, model_obj)
    }
  } else {
    # All mixed infections
    state$mixed <- 1:N
    state$kmin <- 1
    for (imixed in state$mixed) {
      state$A[imixed, categorical_compatible_strains(y[imixed, ], state$D,
                                                     assignments[imixed])] <- 1
    }
    state$loglik <- categorical_loglikelihood_matrix(state$A, state$D, model_obj)
    
    # Resolve exceptions
    iter <- 1
    while (is.infinite(state$loglik)) {
      state <- categorical_resolve_exceptions(state, model_obj)
      iter <- iter + 1
      if (iter > 10) {
        stop("Failed to initialize valid state after 10 iterations")
      }
    }
  }
  
  return(state)
}

#' Resolve exceptions for categorical model
#'
#' @param state Current state
#' @param model_obj Model object
#'
#' @return Updated state
#' @keywords internal
categorical_resolve_exceptions <- function(state, model_obj) {
  y <- model_obj$y
  
  # Calculate current proportions
  ratios <- (state$A %*% state$D) / rowSums(state$A)
  
  # When y == 0, the likelihood is -Inf for any proportion in bin 3, which is
  # `prop > 0.99` and NOT `prop == 1`. Testing equality here missed every
  # proportion in (0.99, 1) -- e.g. 124/125 strains = 0.992 -- leaving an
  # unfixable -Inf that made initialization fail after 10 no-op passes.
  exceptions2 <- which(categorical_prop_bin_vec(ratios) == 3L & y == 0, arr.ind = TRUE)
  if (nrow(exceptions2) > 0) {
    state$A[exceptions2[, 1], state$kmin] <- 1
    state$D[state$kmin, exceptions2[, 2]] <- 0
  }
  
  # Recalculate ratios
  ratios <- (state$A %*% state$D) / rowSums(state$A)
  
  # Likewise, y == 1 is -Inf for proportion bin 1, which is `prop <= 0`.
  exceptions1 <- which(y == 1 & categorical_prop_bin_vec(ratios) == 1L, arr.ind = TRUE)
  if (nrow(exceptions1) > 0) {
    state$A[exceptions1[, 1], state$kmin + 1] <- 1
    state$D[state$kmin + 1, exceptions1[, 2]] <- 1
  }
  
  # Update log-likelihood
  state$loglik <- categorical_loglikelihood_matrix(state$A, state$D, model_obj)
  
  return(state)
}
