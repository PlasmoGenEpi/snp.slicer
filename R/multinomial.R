#' Multinomial Model for SNP-Slice
#'
#' @description
#' Observation model for targets with any number of alleles. Each strain
#' carries exactly one allele per target, stored in the dictionary as an
#' integer allele code (\code{0} for slot 1, \code{1} for slot 2, and so on).
#' For a specimen carrying a set of strains, the expected allele proportions at
#' a target are the fractions of those strains carrying each allele, and the
#' observed read counts across alleles are multinomial with those proportions.
#'
#' With two alleles per target this is exactly the binomial model, so the
#' multinomial model nests it. The dictionary prior is a per-target categorical
#' distribution over alleles rather than the Bernoulli(\code{rho}) prior of the
#' biallelic models.
#'
#' The allocation update, which dominates run time, runs in compiled code via
#' \code{\link{mcmc_kernel_multinomial_cpp}} when the compiled kernel is
#' available; \code{\link{multinomial_update_a_r}} is the R reference. The
#' dictionary update is always the R Gibbs step \code{\link{multinomial_update_d_r}}.
#'
#' @section Expanded layout:
#' Internally, per-target allele counts are stored side by side in one
#' \code{N x L} matrix (\code{counts_exp}), where \code{L} is the total number
#' of allele slots across targets. \code{col_offset[p]} is the number of
#' columns before target \code{p}'s block, \code{col_locus} maps each column to
#' its target, and \code{col_allele} gives the 0-based allele code of each
#' column. [expand_dictionary()] turns an integer dictionary into the matching
#' one-hot \code{K x L} matrix, so \code{A \%*\% Dx} counts, per specimen and
#' allele slot, the strains carrying that allele.
#'
#' @name multinomial_model
#' @keywords internal
NULL

#' Expand an integer dictionary to one-hot allele-slot columns
#'
#' @param D Integer dictionary, strains in rows, one allele code per target.
#' @param model_obj Multinomial model object carrying the expanded layout.
#' @return Binary matrix with \code{nrow(D)} rows and \code{model_obj$L}
#'   columns; row \code{k} has a 1 in the slot of the allele strain \code{k}
#'   carries at each target.
#' @keywords internal
expand_dictionary <- function(D, model_obj) {
  K <- nrow(D)
  P <- ncol(D)
  Dx <- matrix(0, nrow = K, ncol = model_obj$L)
  if (K == 0L) {
    return(Dx)
  }
  idx <- cbind(
    rep(seq_len(K), times = P),
    rep(model_obj$col_offset, each = K) + as.vector(D) + 1
  )
  Dx[idx] <- 1
  Dx
}

#' Log-likelihood for multinomial model (matrix version)
#'
#' Includes the log multinomial coefficients, so values are comparable with
#' the binomial model's \code{dbinom()}-based likelihood on two-allele data.
#'
#' @param A Allocation matrix
#' @param D Integer dictionary matrix
#' @param model_obj Model object containing data
#' @return Log-likelihood value
#' @keywords internal
multinomial_loglikelihood_matrix <- function(A, D, model_obj) {
  start_timer("multinomial_loglikelihood_matrix")

  Dx <- expand_dictionary(D, model_obj)
  props <- (A %*% Dx) / rowSums(A)
  loglik <- multinomial_loglikelihood_vector(props, model_obj$counts_exp) +
    model_obj$loglik_const

  end_timer("multinomial_loglikelihood_matrix")
  as.numeric(loglik)
}

#' Log-likelihood for multinomial model (expanded-vector version)
#'
#' Kernel term only: \code{sum(y * log(prop))} over observed cells with reads.
#' The multinomial coefficient is dropped because this function is used in
#' Metropolis and Gibbs ratios where it cancels. A cell with reads for an
#' allele no assigned strain carries contributes \code{-Inf}.
#'
#' @param propvec Expected allele proportions in the expanded layout (a vector
#'   for one specimen or a matrix for several).
#' @param yvec Observed allele counts in the same layout; \code{NA} marks a
#'   missing genotype.
#' @param rvec Ignored; kept for interface consistency with the other models.
#' @return Log-likelihood value
#' @keywords internal
multinomial_loglikelihood_vector <- function(propvec, yvec, rvec = NULL) {
  obs <- !is.na(yvec) & yvec > 0
  if (!any(obs)) {
    return(0)
  }
  sum(yvec[obs] * log(propvec[obs]))
}

#' Log prior of an integer dictionary under per-target categorical priors
#'
#' @param D Integer dictionary matrix
#' @param model_obj Model object with \code{log_dict_prior_pad}, a
#'   \code{max(n_alleles) x P} matrix of log prior probabilities.
#' @return Log prior probability
#' @keywords internal
multinomial_logprior_d <- function(D, model_obj) {
  K <- nrow(D)
  if (K == 0L) {
    return(0)
  }
  idx <- cbind(as.vector(D) + 1, rep(seq_len(ncol(D)), each = K))
  sum(model_obj$log_dict_prior_pad[idx])
}

#' Draw one categorical value by inversion in index order
#'
#' One uniform, walked along the cumulative weights. The compiled kernel draws
#' the same way, so the R reference and compiled paths consume the RNG
#' identically.
#'
#' @param w Non-negative weights (need not be normalised)
#' @return 0-based index of the drawn category
#' @keywords internal
sample_categorical_inversion <- function(w) {
  u <- stats::runif(1) * sum(w)
  hit <- which(u <= cumsum(w))
  if (length(hit) == 0L) length(w) - 1L else hit[1] - 1L
}

#' Draw one dictionary row from the per-target categorical priors
#'
#' @param model_obj Model object
#' @return Integer vector of allele codes, one per target
#' @keywords internal
multinomial_sample_dict_row <- function(model_obj) {
  vapply(seq_len(model_obj$P), function(p) {
    sample_categorical_inversion(model_obj$dict_prior[[p]])
  }, integer(1))
}

#' Resolve per-target dictionary priors
#'
#' @param dict_prior \code{"empirical"} (pooled allele read fractions with a
#'   pseudocount of one per allele), \code{"uniform"}, or a list of numeric
#'   vectors, one per target, each of length \code{n_alleles[p]}.
#' @param processed_data Processed multi-allelic data.
#' @return List of probability vectors, one per target.
#' @keywords internal
resolve_dict_prior <- function(dict_prior, processed_data) {
  P <- processed_data$P
  n_alleles <- processed_data$n_alleles
  if (is.character(dict_prior) && length(dict_prior) == 1L) {
    if (dict_prior == "uniform") {
      return(lapply(n_alleles, function(m) rep(1 / m, m)))
    }
    if (dict_prior == "empirical") {
      return(lapply(seq_len(P), function(p) {
        cols <- processed_data$col_offset[p] + seq_len(n_alleles[p])
        pooled <- colSums(processed_data$counts_exp[, cols, drop = FALSE], na.rm = TRUE)
        w <- pooled + 1
        w / sum(w)
      }))
    }
    stop("dict_prior must be \"empirical\", \"uniform\", or a list of per-target probability vectors")
  }
  if (!is.list(dict_prior) || length(dict_prior) != P) {
    stop("dict_prior must be a list with one probability vector per target (", P, ")")
  }
  out <- vector("list", P)
  for (p in seq_len(P)) {
    w <- dict_prior[[p]]
    if (!is.numeric(w) || length(w) != n_alleles[p] || any(w < 0) || sum(w) <= 0) {
      stop("dict_prior[[", p, "]] must be ", n_alleles[p], " non-negative weights")
    }
    out[[p]] <- w / sum(w)
  }
  out
}

#' Build the multinomial model object
#'
#' @param processed_data Output of \code{\link{load_dataframe_multiallelic}} or
#'   the list-input path of \code{\link{preprocess_data}} with
#'   \code{model = "multinomial"}.
#' @param alpha IBP concentration parameter
#' @param dict_prior See \code{\link{resolve_dict_prior}}.
#' @return Model object list (class assigned by \code{\link{create_model}}).
#' @keywords internal
build_multinomial_model <- function(processed_data, alpha, dict_prior = "empirical") {
  dict_prior <- resolve_dict_prior(dict_prior, processed_data)
  P <- processed_data$P
  n_alleles <- processed_data$n_alleles
  max_m <- max(n_alleles)
  log_pad <- matrix(-Inf, nrow = max_m, ncol = P)
  prior_pad <- matrix(0, nrow = max_m, ncol = P)
  for (p in seq_len(P)) {
    log_pad[seq_len(n_alleles[p]), p] <- log(dict_prior[[p]])
    prior_pad[seq_len(n_alleles[p]), p] <- dict_prior[[p]]
  }

  # Log multinomial coefficients, summed over observed (specimen, target) cells.
  counts_exp <- processed_data$counts_exp
  r <- processed_data$r
  loglik_const <- sum(lgamma(r + 1), na.rm = TRUE) -
    sum(lgamma(counts_exp + 1), na.rm = TRUE)

  # Per (specimen, target): 0-based code of the most-read allele and its read
  # fraction, used to initialise the dictionary and to spot single infections.
  N <- processed_data$N
  allele_argmax <- matrix(NA_integer_, nrow = N, ncol = P)
  max_allele_frac <- matrix(NA_real_, nrow = N, ncol = P)
  n_observed_alleles <- matrix(0L, nrow = N, ncol = P)
  for (p in seq_len(P)) {
    cols <- processed_data$col_offset[p] + seq_len(n_alleles[p])
    block <- counts_exp[, cols, drop = FALSE]
    has <- !is.na(r[, p]) & r[, p] > 0
    if (any(has)) {
      allele_argmax[has, p] <- max.col(block[has, , drop = FALSE], ties.method = "first") - 1L
      max_allele_frac[has, p] <- apply(block[has, , drop = FALSE], 1L, max) / r[has, p]
      n_observed_alleles[has, p] <- rowSums(block[has, , drop = FALSE] > 0)
    }
  }

  list(
    name = "multinomial",
    y = processed_data$y,
    r = r,
    N = N,
    P = P,
    alpha = alpha,
    rho = NULL,
    L = ncol(counts_exp),
    n_alleles = n_alleles,
    col_offset = processed_data$col_offset,
    col_locus = processed_data$col_locus,
    col_allele = processed_data$col_allele,
    counts_exp = counts_exp,
    loglik_const = loglik_const,
    allele_argmax = allele_argmax,
    max_allele_frac = max_allele_frac,
    n_observed_alleles = n_observed_alleles,
    dict_prior = dict_prior,
    log_dict_prior_pad = log_pad,
    dict_prior_pad = prior_pad,
    loglikelihood_matrix = multinomial_loglikelihood_matrix,
    loglikelihood_vector = multinomial_loglikelihood_vector,
    logprior_d = multinomial_logprior_d,
    sample_dict_row = multinomial_sample_dict_row,
    initialize_state = multinomial_initialize_state,
    resolve_exceptions = multinomial_resolve_exceptions
  )
}

#' Initialize state for multinomial model
#'
#' Each specimen's most-read allele per target gives a candidate strain; the
#' de-duplicated candidates form the initial dictionary. A specimen is a single
#' infection when its top allele holds at least \code{1 - threshold} of the
#' reads at every observed target, and is assigned only its own strain. A mixed
#' specimen starts on every strain compatible with its own calls: strains
#' matching its allele at each target where it shows exactly one allele. Alleles
#' a specimen shows but no assigned strain carries are then repaired by
#' \code{\link{multinomial_resolve_exceptions}}.
#'
#' @param model_obj Model object
#' @param threshold Threshold for identifying single infections
#' @return Initialized state
#' @keywords internal
multinomial_initialize_state <- function(model_obj, threshold = 0.001) {
  N <- model_obj$N
  P <- model_obj$P

  cate <- model_obj$allele_argmax
  cate[is.na(cate)] <- 0L
  storage.mode(cate) <- "double"

  dedup <- remove_duplicate_rows(cate)
  assignments <- dedup$assignments
  nstrain <- nrow(dedup$D)

  state <- list()
  state$D <- dedup$D
  state$A <- matrix(0, nrow = N, ncol = nstrain)

  is_single <- apply(model_obj$max_allele_frac, 1L, function(x) {
    obs <- !is.na(x)
    any(obs) && all(x[obs] >= 1 - threshold)
  })
  which_single <- which(is_single)
  which_mixed <- which(!is_single)

  for (i in which_mixed) {
    state$A[i, multinomial_compatible_strains(i, state$D, model_obj, assignments[i])] <- 1
  }
  for (i in which_single) {
    state$A[i, assignments[i]] <- 1
  }

  if (length(which_single) >= 1L) {
    single_strains <- unique(assignments[which_single])
    ord <- c(single_strains, setdiff(seq_len(nstrain), single_strains))
    state$A <- state$A[, ord, drop = FALSE]
    state$D <- state$D[ord, , drop = FALSE]
    state$kmin <- length(single_strains) + 1
  } else {
    state$kmin <- 1
  }
  state$mixed <- which_mixed

  state$loglik <- multinomial_loglikelihood_matrix(state$A, state$D, model_obj)
  iter <- 1
  while (is.infinite(state$loglik)) {
    state <- multinomial_resolve_exceptions(state, model_obj)
    iter <- iter + 1
    if (iter > 10) {
      stop("Failed to initialize valid state after 10 iterations")
    }
  }

  state
}

#' Strains compatible with a specimen's single-allele calls
#'
#' @param i Specimen index
#' @param D Integer dictionary
#' @param model_obj Model object
#' @param fallback Strain index used when nothing matches
#' @return Integer vector of strain indices
#' @keywords internal
multinomial_compatible_strains <- function(i, D, model_obj, fallback) {
  fixed <- which(model_obj$n_observed_alleles[i, ] == 1L)
  if (length(fixed) == 0L) {
    return(seq_len(nrow(D)))
  }
  target <- model_obj$allele_argmax[i, fixed]
  keep <- which(apply(D[, fixed, drop = FALSE], 1L, function(d) all(d == target)))
  if (length(keep) == 0L) fallback else keep
}

#' Resolve exceptions for multinomial model
#'
#' An exception is a specimen with reads for an allele that none of its
#' assigned strains carries, which gives the likelihood \code{-Inf}. Each is
#' repaired by assigning the specimen an existing strain carrying that allele
#' at that target (preferring updatable strains, \code{k >= kmin}) or, if
#' there is none, appending a new strain equal to the specimen's most-read
#' genotype with that one allele substituted.
#'
#' @param state Current state
#' @param model_obj Model object
#' @return Updated state
#' @keywords internal
multinomial_resolve_exceptions <- function(state, model_obj) {
  Dx <- expand_dictionary(state$D, model_obj)
  props <- (state$A %*% Dx) / rowSums(state$A)
  y <- model_obj$counts_exp
  exceptions <- which(!is.na(y) & y > 0 & props == 0, arr.ind = TRUE)

  for (e in seq_len(nrow(exceptions))) {
    i <- exceptions[e, 1]
    col <- exceptions[e, 2]
    p <- model_obj$col_locus[col]
    m <- model_obj$col_allele[col]

    carriers <- which(state$D[, p] == m)
    if (length(carriers) > 0L) {
      updatable <- carriers[carriers >= state$kmin]
      k <- if (length(updatable) > 0L) updatable[1] else carriers[1]
      state$A[i, k] <- 1
    } else {
      new_row <- model_obj$allele_argmax[i, ]
      new_row[is.na(new_row)] <- 0L
      new_row[p] <- m
      state$D <- rbind(state$D, as.numeric(new_row))
      state$A <- cbind(state$A, 0)
      state$A[i, ncol(state$A)] <- 1
    }
  }

  state$loglik <- multinomial_loglikelihood_matrix(state$A, state$D, model_obj)
  state
}

#' Allocation update for the multinomial model (R kernel)
#'
#' Reuses \code{\link{slice_update_a_r}} on the expanded layout: the integer
#' dictionary is swapped for its one-hot expansion and the observation matrix
#' for the expanded allele counts, so the generic per-cell Metropolis step
#' evaluates the multinomial likelihood unchanged.
#'
#' @param state Current state
#' @param model_obj Model object
#' @return Updated state
#' @keywords internal
multinomial_update_a_r <- function(state, model_obj) {
  D_int <- state$D
  state$D <- expand_dictionary(D_int, model_obj)
  model_obj$y <- model_obj$counts_exp
  model_obj$r <- model_obj$counts_exp
  state <- slice_update_a_r(state, model_obj)
  state$D <- D_int
  state
}

#' Log weights for the allele options of one dictionary cell
#'
#' Pure function behind the Gibbs step in
#' \code{\link{multinomial_update_d_r}}. Setting the cell to allele \code{m}
#' adds the strain's carriers to slot \code{m} only, so relative to a constant
#' shared by every option the log conditional of \code{m} is
#' \code{log_prior[m] + sum_h y[h, m] * (log(base[h, m] + 1) - log(base[h, m]))}.
#' The normalisation by each specimen's strain count is the same for every
#' option and drops out.
#'
#' A slot with reads in some carrying specimen and no other carrier there
#' (\code{base == 0} with \code{y > 0}) makes every other option impossible:
#' that allele gets weight 0 on the log scale and the rest \code{-Inf}. In a
#' valid state only the current allele can be such a slot.
#'
#' @param base Leave-strain-out carrier counts for the carrying specimens
#'   (specimens x alleles at this target).
#' @param yb Observed allele counts for those specimens at this target.
#' @param a_h The strain's allocation value for each of those specimens.
#' @param log_prior Log prior over the target's alleles.
#' @return Numeric vector of log weights, one per allele.
#' @keywords internal
multinomial_cell_logweights <- function(base, yb, a_h, log_prior) {
  has_reads <- !is.na(yb) & yb > 0
  required <- which(colSums(has_reads & base <= 0) > 0)
  if (length(required) > 0L) {
    lw <- rep(-Inf, length(log_prior))
    if (length(required) == 1L) {
      lw[required] <- 0
    }
    return(lw)
  }
  vapply(seq_along(log_prior), function(m) {
    obs <- has_reads[, m]
    if (!any(obs)) {
      return(log_prior[m])
    }
    c0 <- base[obs, m]
    log_prior[m] + sum(yb[obs, m] * (log(c0 + a_h[obs]) - log(c0)))
  }, numeric(1))
}

#' Dictionary update for the multinomial model (R kernel)
#'
#' Exact Gibbs step over each updatable dictionary cell: for strain \code{k}
#' and target \code{p}, the allele is redrawn from its full conditional given
#' every other cell, the current allocations, and the per-target prior. Only
#' specimens carrying strain \code{k} enter the likelihood ratio, and only
#' the slot each option adds carriers to (see
#' \code{\link{multinomial_cell_logweights}}); the specimen-by-slot carrier
#' count matrix is kept incrementally. A cell with a single feasible allele is
#' set to it without a draw, and a strain no specimen carries is redrawn from
#' the prior.
#'
#' @param state Current state
#' @param model_obj Model object
#' @return Updated state
#' @keywords internal
multinomial_update_d_r <- function(state, model_obj) {
  A <- state$A
  D <- state$D
  M <- A %*% expand_dictionary(D, model_obj)
  y <- model_obj$counts_exp

  # seq() rather than kmin:kstar: when every strain belongs to a single
  # infection kmin exceeds kstar and the colon would count downward.
  if (state$kstar < state$kmin) {
    return(state)
  }
  for (k in seq(state$kmin, state$kstar)) {
    hosts <- which(A[, k] > 0)
    if (length(hosts) == 0L) {
      D[k, ] <- multinomial_sample_dict_row(model_obj)
      next
    }
    a_h <- A[hosts, k]
    for (p in seq_len(model_obj$P)) {
      cols <- model_obj$col_offset[p] + seq_len(model_obj$n_alleles[p])
      cur <- D[k, p]
      base <- M[hosts, cols, drop = FALSE]
      base[, cur + 1] <- base[, cur + 1] - a_h
      lw <- multinomial_cell_logweights(
        base, y[hosts, cols, drop = FALSE], a_h,
        model_obj$log_dict_prior_pad[seq_along(cols), p]
      )
      finite <- which(is.finite(lw))
      if (length(finite) == 0L) {
        next
      }
      new <- if (length(finite) == 1L) {
        finite - 1L
      } else {
        sample_categorical_inversion(exp(lw - max(lw)))
      }
      if (new != cur) {
        D[k, p] <- new
        M[hosts, cols[cur + 1]] <- M[hosts, cols[cur + 1]] - a_h
        M[hosts, cols[new + 1]] <- M[hosts, cols[new + 1]] + a_h
      }
    }
  }

  state$D <- D
  state
}

#' De-duplicate rows of a matrix of allele codes
#'
#' Like \code{\link{remove_duplicates}} but keyed on exact row equality, so it
#' works for integer allele codes as well as binary patterns.
#'
#' @param dd Numeric matrix, one candidate strain per row
#' @return List with \code{assignments} (row -> unique-row index) and \code{D}
#'   (the unique rows, in order of first appearance)
#' @keywords internal
remove_duplicate_rows <- function(dd) {
  rownames(dd) <- NULL
  keys <- apply(dd, 1L, paste, collapse = ",")
  unique_keys <- unique(keys)
  list(
    assignments = match(keys, unique_keys),
    D = dd[match(unique_keys, keys), , drop = FALSE]
  )
}
