# Global variables for data loading and dplyr/tidyr NSE
utils::globalVariables(c(
  "example_snp_data",
  "target_id", "target_value", "specimen_id", "target_count",
  "total_count", "target_idx"
))

#' Validate input data for SNP-Slice
#'
#' @param data Input data to validate
#' @param model Model type to validate against
#' @param ... Passed through; used for column-name args when validating categorical data.frames
#'
#' @return TRUE if valid, otherwise throws error
#' @keywords internal
validate_input_data <- function(data, model, ...) {
  
  # Check if data is NULL or empty
  if (is.null(data) || (is.list(data) && length(data) == 0)) {
    stop("Data cannot be NULL or empty")
  }
  
  # Handle different data types (check data.frame before list, since data.frames are lists in R)
  if (is.character(data)) {
    # File path - validate file exists
    if (!file.exists(data)) {
      stop("File not found: ", data)
    }
    return(TRUE)
  }

  if (is.data.frame(data)) {
    params <- list(
      target_id_col = "target_id",
      target_value_col = "target_value",
      specimen_id_col = "specimen_id",
      target_count_col = "target_count"
    )
    overrides <- list(...)
    for (nm in names(overrides)) if (nm %in% names(params)) params[[nm]] <- overrides[[nm]]

    # All long-format data.frames (categorical and count) need the same four columns.
    # Categorical uses target_count too: counts_to_categorical derives 0/0.5/1/NA from
    # whether each allele's summed count is > 0.
    required <- c(params$target_id_col, params$target_value_col, params$specimen_id_col, params$target_count_col)
    missing <- setdiff(required, names(data))
    if (length(missing) > 0) {
      stop("Long-format data.frame must have columns: ", paste(required, collapse = ", "), "; missing: ", paste(missing, collapse = ", "))
    }
    if (nrow(data) == 0) {
      stop("Long-format data.frame cannot be empty")
    }
    count_col <- params$target_count_col
    if (is.numeric(data[[count_col]]) && any(data[[count_col]] < 0, na.rm = TRUE)) {
      stop("target_count cannot be negative")
    }

    # Structural checks shared by all long-format data.frames (categorical or count)
    validate_long_dataframe(data, params)
    return(TRUE)
  }

  if (is.matrix(data)) {
    # Direct matrix/data.frame input
    if (nrow(data) == 0 || ncol(data) == 0) {
      stop("Data matrix cannot be empty")
    }
    
    if (model == "categorical") {
      # Check for valid categorical values
      valid_values <- c(0, 0.5, 1, NA)
      if (!all(as.matrix(data) %in% valid_values, na.rm = TRUE)) {
        stop("Categorical data must contain only values 0, 0.5, 1, or NA")
      }
    } else {
      # Check for numeric data
      if (!is.numeric(as.matrix(data))) {
        stop("Data must be numeric for non-categorical models")
      }
      if (any(as.matrix(data) < 0, na.rm = TRUE)) {
        stop("Data cannot contain negative values")
      }
    }
    return(TRUE)
  }
  
  if (is.list(data)) {
    # List format for read count data
    if (model == "categorical") {
      stop("Categorical model expects matrix/data.frame, not list")
    }

    # Check for required elements
    if (!all(c("read1", "read0") %in% names(data))) {
      stop("List data must contain 'read1' and 'read0' elements")
    }

    # Validate read1 and read0
    read1 <- data$read1
    read0 <- data$read0

    if (!is.matrix(read1) || !is.matrix(read0)) {
      stop("read1 and read0 must be matrices")
    }

    if (!identical(dim(read1), dim(read0))) {
      stop("read1 and read0 must have identical dimensions")
    }

    if (any(read1 < 0, na.rm = TRUE) || any(read0 < 0, na.rm = TRUE)) {
      stop("Read counts cannot be negative")
    }

    return(TRUE)
  }

  stop("Unsupported data type. Expected matrix, data.frame, list, or file path")
}

#' Structural validation for long-format data.frames
#'
#' Checks shared by every long-format input (categorical and count models): each
#' allele must appear at most once per specimen and target, and at least one
#' target must be biallelic (exactly two observed alleles) so the model has
#' variation to resolve strains. Monomorphic targets are allowed alongside
#' biallelic ones; targets with more than two alleles are dropped downstream by
#' \code{\link{load_dataframe}}. The caller (\code{\link{validate_input_data}})
#' guarantees the key columns are present before this is called.
#'
#' @param data Long-format data.frame.
#' @param params List of resolved column names with elements
#'   \code{specimen_id_col}, \code{target_id_col}, \code{target_value_col}.
#' @return Invisibly \code{TRUE}; otherwise throws an error.
#' @keywords internal
validate_long_dataframe <- function(data, params) {
  sid <- params$specimen_id_col
  tid <- params$target_id_col
  tval <- params$target_value_col
  tcol <- params$target_count_col

  # No duplicate allele rows: each (specimen, target, value) at most once.
  if (anyDuplicated(data[, c(sid, tid, tval), drop = FALSE]) > 0) {
    stop("Input data contains duplicate (", sid, ", ", tid, ", ", tval,
         ") rows; each allele must appear at most once per specimen and target")
  }

  # At least one biallelic locus (exactly two observed alleles) must survive the
  # downstream <=2-allele filter, or there is no variation to model.
  alleles_per_target <- tapply(data[[tval]], data[[tid]], function(v) length(unique(v)))
  if (length(alleles_per_target) == 0 || !any(alleles_per_target == 2)) {
    stop("Input data has no biallelic loci (targets with exactly two observed alleles)")
  }

  # At least one biallelic locus must have positive total reads
  biallelic_targets <- names(alleles_per_target)[alleles_per_target == 2]
  biallelic_rows <- data[[tid]] %in% biallelic_targets
  total_by_biallelic <- tapply(data[[tcol]][biallelic_rows], data[[tid]][biallelic_rows], sum, na.rm = TRUE)
  if (!any(total_by_biallelic > 0)) {
    stop("All biallelic loci have zero or missing total reads")
  }

  invisible(TRUE)
}

#' Validate algorithm parameters
#'
#' @param alpha IBP concentration parameter
#' @param rho Dictionary sparsity parameter
#' @param threshold Single infection threshold
#'
#' @return TRUE if valid, otherwise throws error
#' @keywords internal
validate_parameters <- function(alpha, rho, threshold) {
  
  if (!is.numeric(alpha) || alpha <= 0) {
    stop("alpha must be a positive number")
  }
  
  if (!is.null(rho) && (!is.numeric(rho) || rho < 0 || rho > 1)) {
    stop("rho must be a number between 0 and 1")
  }
  
  if (!is.numeric(threshold) || threshold < 0 || threshold > 1) {
    stop("threshold must be a number between 0 and 1")
  }
  
  return(TRUE)
}

#' Validate MCMC settings
#'
#' @param n_sample Number of post-burn-in iterations to retain
#' @param n_burnin Number of iterations discarded before sampling begins
#' @param gap Early stopping threshold
#'
#' @return TRUE if valid, otherwise throws error
#' @keywords internal
validate_mcmc_settings <- function(n_sample, n_burnin, gap, n_chains = 1, n_cores = 1) {

  if (!is.numeric(n_sample) || n_sample <= 0 || n_sample != as.integer(n_sample)) {
    stop("n_sample must be a positive integer")
  }

  if (!is.null(n_burnin)) {
    if (!is.numeric(n_burnin) || n_burnin < 0 || n_burnin != as.integer(n_burnin)) {
      stop("n_burnin must be a non-negative integer")
    }
  }

  if (!is.null(gap)) {
    if (!is.numeric(gap) || gap <= 0 || gap != as.integer(gap)) {
      stop("gap must be a positive integer")
    }
  }

  if (!is.numeric(n_chains) || length(n_chains) != 1 || n_chains <= 0 ||
        n_chains != as.integer(n_chains)) {
    stop("n_chains must be a positive integer")
  }

  if (!is.numeric(n_cores) || length(n_cores) != 1 || n_cores <= 0 ||
        n_cores != as.integer(n_cores)) {
    stop("n_cores must be a positive integer")
  }

  return(TRUE)
}

#' Convert read count matrices to categorical observation matrix
#'
#' Maps reference (read0) and alternate (read1) counts to categorical states:
#' 0 (ref only), 1 (alt only), 0.5 (both), or NA (zero total).
#'
#' @param read0 Matrix of reference allele counts (specimens x targets).
#' @param read1 Matrix of alternate allele counts (same dimensions as read0).
#' @return Numeric matrix of same shape with values in \{0, 0.5, 1, NA\}, dimnames preserved.
#' @keywords internal
counts_to_categorical <- function(read0, read1) {
  if (!identical(dim(read0), dim(read1))) {
    stop("read0 and read1 must have identical dimensions")
  }
  if (any(read0 < 0, na.rm = TRUE) || any(read1 < 0, na.rm = TRUE)) {
    stop("Read counts cannot be negative")
  }
  total <- read0 + read1
  y <- matrix(NA_real_, nrow = nrow(read0), ncol = ncol(read0), dimnames = dimnames(read0))
  y[total > 0 & read0 > 0 & read1 == 0] <- 0
  y[total > 0 & read0 == 0 & read1 > 0] <- 1
  y[total > 0 & read0 > 0 & read1 > 0] <- 0.5
  y
}

#' Load data from a dataframe
#'
#' @param data Input dataframe with columns:
#' \describe{
#'   \item{specimen_id}{Specimen ID}
#'   \item{target_id}{Target ID}
#'   \item{target_value}{Target value (allele)}
#'   \item{target_count}{Target count}
#' }
#'   Input is assumed to have passed \code{\link{validate_input_data}}: the four
#'   columns are present, no \code{(specimen, target, value)} row is duplicated, and
#'   at least one target is biallelic.
#' @param model Model type
#' @param target_id_col Name of the target ID column
#' @param target_value_col Name of the target value column
#' @param specimen_id_col Name of the specimen ID column
#' @param target_count_col Name of the target count column
#'
#' @return Processed data list with y, r, and metadata
#'
#' @details The alleles for each target are sorted and given indices such that 
#'   target_idx 1 (read0) is the minor allele and target_idx 2 (read1) 
#'   the major. For monomorphic targets, with only one observed allele, both 
#'   slots carry that same label and read1 is zero (because no second allele 
#'   was observed). A specimen with no reads at a target (i.e., total_count 0) 
#'   is encoded as \code{NA} for both slots, thus distinguishing a missing 
#'   genotype from a homozygous call.
#'
#'   When \code{model == "categorical"}, long-format data.frames are handled by
#'   \code{load_dataframe_categorical}, which builds read0/read1 matrices from the same
#'   layout and then converts counts to categorical observations.  For categorical
#'   data the allele order is determined by descending \code{target_value} (not by
#'   count), so the lexicographically larger allele label occupies target_idx 1
#'   (read0).  A specimen observed with only that allele yields \code{y = 0}, with
#'   only the other allele \code{y = 1}, with both \code{y = 0.5}, and with no reads
#'   \code{y = NA}.
#'
#' @keywords internal
load_dataframe <- function(
  data, model, target_id_col = "target_id", target_value_col = "target_value", specimen_id_col = "specimen_id", target_count_col = "target_count") {
  data_renamed <- data |>
    dplyr::rename(
      target_id = !!target_id_col,
      target_value = !!target_value_col,
      specimen_id = !!specimen_id_col,
      target_count = !!target_count_col
    ) |>
    dplyr::select(target_id, target_value, specimen_id, target_count)

  # Order alleles per target, remove loci with more than two alleles, 
  # and assign target indices
  target_alleles <- data_renamed |>
    dplyr::group_by(target_id, target_value) |>
    dplyr::summarize(total_count = sum(target_count), .groups = "drop") |>
    dplyr::group_by(target_id) |>
    dplyr::filter(dplyr::n() <= 2) |>
    dplyr::arrange(
      # SNP-Slice assumes the minor allele is in the first slot
      if (model != "categorical") total_count else 0, 
      dplyr::desc(target_value), 
      .by_group = TRUE
    ) |>
    dplyr::mutate(target_idx = dplyr::row_number()) |>
    dplyr::select(-total_count) |>
    dplyr::ungroup()

  # Pivot
  target_alleles_wide <- target_alleles |>
    tidyr::pivot_wider(names_from = target_idx, values_from = target_value, names_prefix = "a")
  # Monomorphic targets have no target_idx 2, so duplicate the target_idx 1 label into a2
  target_alleles_wide$a2 <- dplyr::coalesce(target_alleles_wide$a2, target_alleles_wide$a1)

  # Observed counts per (specimen, target, target_idx)
  counts_by_idx <- data_renamed |>
    # inner_join maps each observation to exactly one target_idx and drops >2-allele loci.
    dplyr::inner_join(target_alleles, by = c("target_id", "target_value")) |>
    dplyr::select(specimen_id, target_id, target_idx, target_count)

  # Complete the specimen x target x target_idx grid.
  completed_data <- counts_by_idx |>
    # Fill in second target_idx of monomorphic loci with count of 0
    tidyr::complete(specimen_id, target_id, target_idx = 1:2, fill = list(target_count = 0)) |>
    dplyr::group_by(specimen_id, target_id) |>
    # Set count to NA for missing loci
    dplyr::mutate(total_count = sum(target_count)) |>
    dplyr::ungroup() |>
    dplyr::mutate(target_count = ifelse(total_count == 0, NA, target_count)) |>
    dplyr::select(-total_count) |>
    dplyr::arrange(specimen_id, target_id)

  read0_df <- completed_data |> dplyr::filter(target_idx == 1)
  read1_df <- completed_data |> dplyr::filter(target_idx == 2)

  specimen_ids <- unique(read0_df$specimen_id)
  target_ids <- unique(read0_df$target_id)
  nspecs <- length(specimen_ids)
  ntars <- length(target_ids)

  read0_mat <- matrix(read0_df$target_count, nrow = nspecs, ncol = ntars, dimnames = list(specimen_ids, target_ids), byrow = TRUE)
  read1_mat <- matrix(read1_df$target_count, nrow = nspecs, ncol = ntars, dimnames = list(specimen_ids, target_ids), byrow = TRUE)

  # Per-column allele labels (target_idx 2 == target_idx 1 for monomorphic loci).
  r0_values_ordered <- target_alleles_wide$a1[match(target_ids, target_alleles_wide$target_id)]
  r1_values_ordered <- target_alleles_wide$a2[match(target_ids, target_alleles_wide$target_id)]

  return(list(
    y = read0_mat,
    r = read0_mat + read1_mat,
    N = nspecs,
    P = ntars,
    model = model,
    data_type = "read_counts",
    target_ids = target_ids,
    specimen_ids = specimen_ids,
    r0_values = r0_values_ordered,
    r1_values = r1_values_ordered
  ))

}

#' Load long-format data.frame for categorical model
#'
#' Uses the same join/completion logic as \code{\link{load_dataframe}} to build
#' read0/read1 matrices, then converts counts to categorical observations (0, 0.5, 1, NA)
#' via \code{\link{counts_to_categorical}}.
#'
#' @param data Input dataframe (same column semantics as \code{load_dataframe}).
#' @param model Model type (must be \code{"categorical"}).
#' @param target_id_col,target_value_col,specimen_id_col,target_count_col
#'   Same as \code{\link{load_dataframe}}.
#' @return Processed data list with \code{data_type = "categorical"}, \code{y} the categorical
#'   matrix, and \code{r = NULL}.
#' @keywords internal
load_dataframe_categorical <- function(
  data, model, target_id_col = "target_id",
  target_value_col = "target_value", specimen_id_col = "specimen_id",
  target_count_col = "target_count", ...) {
  raw <- load_dataframe(
    data, model,
    target_id_col = target_id_col,
    target_value_col = target_value_col,
    specimen_id_col = specimen_id_col,
    target_count_col = target_count_col
  )
  read0_mat <- raw$y
  read1_mat <- raw$r - raw$y
  y_cat <- counts_to_categorical(read0_mat, read1_mat)
  list(
    y = y_cat,
    r = NULL,
    N = raw$N,
    P = raw$P,
    model = model,
    data_type = "categorical",
    target_ids = raw$target_ids,
    specimen_ids = raw$specimen_ids,
    r0_values = raw$r0_values,
    r1_values = raw$r1_values
  )
}

#' Preprocess data for SNP-Slice
#'
#' @param data Input data
#' @param model Model type
#'
#' @return Processed data list with y, r, and metadata
#' @keywords internal
preprocess_data <- function(data, model, ...) {

  if (is.data.frame(data)) {
    if (model == "categorical") {
      return(load_dataframe_categorical(data, model, ...))
    }
    return(load_dataframe(data, model, ...))
  }
 
  if (is.matrix(data)) {
    # Direct matrix input
    if (is.null(rownames(data))) {
      rownames(data) <- paste0("specimen_", seq_len(nrow(data)))
    }
    if (is.null(colnames(data))) {
      colnames(data) <- paste0("target_", seq_len(ncol(data)))
    }
    
    if (model == "categorical") {
      # Categorical data - no r matrix needed
      return(list(
        y = data,
        r = NULL,
        N = nrow(data),
        P = ncol(data),
        data_type = "categorical",
        model = model,
        target_ids = colnames(data),
        specimen_ids = rownames(data),
        r0_values = rep("ref", ncol(data)),
        r1_values = rep("alt", ncol(data))
      ))
    } else {
      # For other models, assume y contains proportions or counts
      # This is a simplified approach - in practice, you'd want more sophisticated handling
      stop("Matrix input for non-categorical models requires read count data in list format")
    }
  }

  if (is.list(data)) {
    # Read count data
    read1 <- as.matrix(data$read1)
    read0 <- as.matrix(data$read0)
    
    y <- read1
    r <- read1 + read0

    if (is.null(rownames(y))) {
      rownames(y) <- paste0("specimen_", seq_len(nrow(y)))
    }
    if (is.null(colnames(y))) {
      colnames(y) <- paste0("target_", seq_len(ncol(y)))
    }
    
    return(list(
      y = y,
      r = r,
      N = nrow(y),
      P = ncol(y),
      data_type = "read_counts",
      model = model,
      target_ids = colnames(y),
      specimen_ids = rownames(y),
      r0_values = rep("ref", ncol(y)),
      r1_values = rep("alt", ncol(y))
    ))
  }

  if (is.character(data)) {
    # Categorical files (e.g. example_cat.txt) are wide: host_id, strain_id, site_1..site_N
    if (model == "categorical" || grepl("_cat\\.txt$", data)) {
      mat <- read_snp_data(data, format = "categorical")
      return(preprocess_data(mat, model, ...))
    }
    df <- readr::read_tsv(data)
    return(load_dataframe(df, model, ...))
  }
  
  stop("Unable to preprocess data")
}

#' Read SNP data from file
#'
#' @param file_path Path to data file
#' @param format Data format ("auto", "read_counts", "categorical")
#'
#' @return Processed data list
#' @keywords internal
read_snp_data <- function(file_path, format = "auto") {
  
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  # Try to determine format from file extension or content
  if (format == "auto") {
    if (grepl("_cat\\.txt$", file_path)) {
      format <- "categorical"
    } else if (grepl("_read[01]\\.txt$", file_path)) {
      format <- "read_counts"
    } else {
      # Default to read counts
      format <- "read_counts"
    }
  }
  
  if (format == "categorical") {
    # Read categorical data
    data <- utils::read.delim(file_path, stringsAsFactors = FALSE)
    
    # Remove ID columns if present
    id_cols <- c("host_id", "strain_id")
    data <- data[, !names(data) %in% id_cols, drop = FALSE]
    
    return(as.matrix(data))
  }
  
  if (format == "read_counts") {
    # For read count data, we need both read1 and read0 files
    # This is a simplified implementation
    stop("Reading read count data from files requires both read1 and read0 files")
  }
  
  stop("Unsupported format: ", format)
}

#' Create results object
#'
#' @param mcmc_result MCMC results
#' @param model_obj Model object
#' @param processed_data Processed data
#'
#' @return snp_slice_results object
#' @keywords internal
create_results_object <- function(mcmc_result, model_obj, processed_data) {

  chains <- lapply(mcmc_result$chains, function(chain) {
    chain_maps <- map_estimates(chain$map_state)
    list(
      chain_id = chain$chain_id,
      seed = chain$seed,
      allocation_matrix = chain_maps$allocation_matrix,
      dictionary_matrix = chain_maps$dictionary_matrix,
      mcmc_samples = chain$samples,
      diagnostics = chain$diagnostics,
      convergence = chain$convergence,
      performance = chain$performance
    )
  })

  # Create results object
  results <- list(
    chains = chains,
    best_chain = mcmc_result$best_chain,
    parameters = mcmc_result$parameters,
    model_info = list(
      model = model_obj$name,
      data_type = processed_data$data_type,
      processed_data = processed_data
    )
  )

  class(results) <- "snp_slice_results"

  return(results)
}

#' Extract MAP allocation and dictionary matrices from a chain state
#'
#' @param map_state MAP state of a chain
#'
#' @return List with \code{allocation_matrix} and \code{dictionary_matrix},
#'   restricted to strains with at least one host
#' @keywords internal
map_estimates <- function(map_state) {
  active_strains <- colSums(map_state$A) > 0
  list(
    allocation_matrix = map_state$A[, active_strains, drop = FALSE],
    dictionary_matrix = map_state$D[active_strains, , drop = FALSE]
  )
}

#' Load Example Analysis Results
#'
#' Loads pre-computed SNP-Slice analysis results from the example data.
#' These results were generated using the negative binomial model with 1000
#' burn-in iterations followed by 1000 retained samples.
#'
#' @return A \code{snp_slice_results} object containing the analysis results
#' @export
#' @examples
#' # Load the pre-computed results
#' result <- load_example_results()
#' print(result)
load_example_results <- function() {
  results_file <- system.file("extdata", "example_analysis_results.rds", package = "snp.slicer")
  
  if (!file.exists(results_file)) {
    warning("Example analysis results not found. Running analysis now (this may take a few minutes)...")
    
    # Load example data
    utils::data(example_snp_data, package = "snp.slicer", envir = environment())
    example_data <- example_snp_data
    data <- list(read0 = example_data$read0, read1 = example_data$read1)
    
    # Run analysis
    set.seed(123)
    result <- snp_slice(data,
                       model = "negative_binomial",
                       n_sample = 1000,
                       n_burnin = 1000,
                       store_mcmc = TRUE,
                       verbose = FALSE)
    
    return(result)
  }
  
  readRDS(results_file)
}

#' Example SNP Data
#'
#' Example SNP data for testing and demonstration purposes.
#'
#' @format A list containing:
#' \describe{
#'   \item{read0}{Matrix of reference allele read counts}
#'   \item{read1}{Matrix of alternate allele read counts}
#' }
#' @source Simulated data for package testing
#' @keywords datasets
"example_snp_data"

#' Example Read Count Data - Reference Alleles
#'
#' Example reference allele read count data for testing and demonstration purposes.
#'
#' @format A matrix of reference allele read counts with 200 hosts and 96 SNPs
#' @source Simulated data for package testing
#' @keywords datasets
"example_read0"

#' Example Read Count Data - Alternate Alleles
#'
#' Example alternate allele read count data for testing and demonstration purposes.
#'
#' @format A matrix of alternate allele read counts with 200 hosts and 96 SNPs
#' @source Simulated data for package testing
#' @keywords datasets
"example_read1"

#' Calculate Allele Frequencies from MCMC Results
#'
#' Calculates allele frequencies for a collection of SNPs treated as a single allele.
#' The function takes SNP indices and calculates the frequency of each possible
#' allele combination across all individuals based on their strain allocations.
#'
#' @param results A \code{snp_slice_results} object containing MCMC results
#' @param snp_indices A vector of SNP indices to treat as a single allele
#' @param use_map Logical, whether to use MAP estimates (TRUE) or sample from MCMC (FALSE)
#' @param n_samples Number of MCMC samples to use if use_map = FALSE (default: 100)
#' @param interval Numeric in (0, 1). Credible interval width when using MCMC (e.g. 0.95).
#' @param allele_sep Separator for allele strings (default: "|")
#' @param additional_burnin Number of additional stored samples to discard from
#'   the start of the retained chain. The retained chain is already
#'   post-burn-in, so this defaults to 0.
#' @param chain Which chain of a multi-chain run to use. \code{NULL} (default)
#'   uses the top-level estimates, i.e. the chain with the highest MAP log
#'   posterior (\code{results$best_chain}). Give an index to analyse a specific
#'   chain instead; see \code{\link{compare_chains}}.
#'
#' @return The structure depends on \code{use_map}. \describe{
#'   \item{MAP (\code{use_map = TRUE})}{Data frame with columns: \code{allele} (string representation of the allele, e.g. \code{"ref|alt|ref"} for 3 SNPs), \code{frequency} (proportion of total parasites with this allele; sums to 1), \code{count} (number of parasites with this allele in the MAP allocation), \code{total_parasites} (total parasites in the MAP allocation; same for every row).}
#'   \item{MCMC (\code{use_map = FALSE})}{Data frame with columns: \code{allele}, \code{frequency} (posterior mean proportion), \code{frequency_sd} (posterior SD of proportion across samples), \code{frequency_lower} and \code{frequency_upper} (credible interval, e.g. 2.5\% and 97.5\%), \code{mean_count} (mean parasites with this allele per MCMC sample; does not scale with \code{n_samples}), \code{n_samples}. Attribute \code{mean_total_parasites}: mean total parasites per MCMC sample (same for all alleles).}
#' }
#'
#' @details With \code{use_map = FALSE}, counts are summarized as per-sample means rather than sums,
#'   so \code{mean_count} and \code{mean_total_parasites} are interpretable regardless of \code{n_samples}.
#'   Frequency uncertainty (\code{frequency_sd}, \code{frequency_lower}, \code{frequency_upper}) is
#'   computed from the distribution of allele frequencies across MCMC samples.
#'
#' @export
#' @examples
#' # Load example results
#' result <- load_example_results()
#' 
#' # Calculate allele frequencies for SNPs 1, 5, and 10 (MAP)
#' allele_freqs <- calculate_allele_frequencies(result, c(1, 5, 10))
#' print(allele_freqs)
#'
#' # With MCMC: posterior mean, SD, credible interval, and per-sample mean count
#' if (!is.null(get_chain(result)$mcmc_samples)) {
#'   allele_freqs_mcmc <- calculate_allele_frequencies(result, c(1, 5, 10), use_map = FALSE, n_samples = 50)
#'   print(allele_freqs_mcmc)
#' }
calculate_allele_frequencies <- function(results, snp_indices, use_map = TRUE, n_samples = 100, interval = 0.95, allele_sep = "|", additional_burnin = 0, chain = NULL) {
  if (length(unique(snp_indices)) != length(snp_indices)) {
    stop("snp_indices must be unique")
  }

  if (inherits(results, "snp_slice_results")) {
    results <- resolve_chain(results, chain)
  }

  if (is.character(snp_indices)) {
    found_snps <- match(snp_indices, results$model_info$processed_data$target_ids)
    if (any(is.na(found_snps))) {
      stop("Some SNPs not found in the data: ", paste(snp_indices[is.na(found_snps)], collapse = ", "))
    }
    snp_indices <- found_snps
  }
  
  # Validate inputs
  if (!inherits(results, "snp_slice_results")) {
    stop("results must be a snp_slice_results object")
  }
  
  if (length(snp_indices) == 0) {
    stop("snp_indices must be a non-empty vector")
  }
  
  if (any(snp_indices < 1) || any(snp_indices > ncol(results$dictionary_matrix))) {
    stop("snp_indices must be valid SNP positions (1 to ", ncol(results$dictionary_matrix), ")")
  }
  
  if (!is.logical(use_map) || length(use_map) != 1) {
    stop("use_map must be a single logical value")
  }
  
  if (!is.numeric(n_samples) || n_samples < 1) {
    stop("n_samples must be a positive integer")
  }

  if (!is.numeric(interval) || interval <= 0 || interval >= 1) {
    stop("interval must be a number between 0 and 1 (exclusive)")
  }
  
  r0_values <- results$model_info$processed_data$r0_values[snp_indices]
  r1_values <- results$model_info$processed_data$r1_values[snp_indices] 

  if (use_map) {
    # Use MAP estimates: single allocation, return count and total_parasites
    A <- results$allocation_matrix
    D <- results$dictionary_matrix
    allele_counts <- calculate_allele_counts_single(A, D, snp_indices, r0_values, r1_values, sep = allele_sep)
    total_parasites <- sum(allele_counts$count)
    result_df <- data.frame(
      allele = allele_counts$allele,
      frequency = allele_counts$count / total_parasites,
      count = allele_counts$count,
      total_parasites = total_parasites,
      stringsAsFactors = FALSE
    )
  } else {
    # MCMC: per-sample frequencies and counts, then summarize (sample-size invariant)
    if (is.null(results$mcmc_samples)) {
      stop("MCMC samples not available. Set use_map = TRUE or run snp_slice with store_mcmc = TRUE")
    }
    samples <- retained_samples(results, additional_burnin)
    n_total_samples <- length(samples)
    if (n_samples > n_total_samples) {
      warning("Requested ", n_samples, " samples but only ", n_total_samples, " retained. Using all retained samples.")
      n_samples <- n_total_samples
    }
    sample_indices <- sample(seq_len(n_total_samples), n_samples, replace = FALSE)
    allele_counts_list <- lapply(sample_indices, function(i) {
      sample_data <- samples[[i]]
      calculate_allele_counts_single(sample_data$A, sample_data$D, snp_indices, r0_values, r1_values, sep = allele_sep)
    })
    result_df <- summarize_allele_frequencies_mcmc(allele_counts_list, interval = interval, n_samples = n_samples)
  }
  
  result_df <- result_df[order(result_df$frequency, decreasing = TRUE), ]
  return(result_df)
}

#' Calculate allele frequencies for multiple target sets
#'
#' For each target set, computes allele frequencies from MCMC/MAP results and
#' returns a list of frequency tables (one per set). Each set is a vector of
#' target indices or target names.
#'
#' @param results A \code{snp_slice_results} object containing MCMC results.
#' @param target_sets List of vectors; each element is target indices (integer)
#'   or target names (character) defining one set. If the list is named, those
#'   names are used for the returned list.
#' @param ... Arguments passed on to \code{\link{calculate_allele_frequencies}}.
#'
#' @return A named list of data frames, one per target set. List names come from
#'   \code{names(target_sets)} or \code{"set_1"}, \code{"set_2"}, etc. Each data frame
#'   has the same structure as the return value of \code{\link{calculate_allele_frequencies}}:
#'   with MAP, columns \code{allele}, \code{frequency}, \code{count}, \code{total_parasites};
#'   with MCMC, columns \code{allele}, \code{frequency}, \code{frequency_sd},
#'   \code{frequency_lower}, \code{frequency_upper}, \code{mean_count}, \code{n_samples},
#'   and attribute \code{mean_total_parasites}. See that function's help for the meaning
#'   of each column.
#'
#' @export
#' @examples
#' result <- load_example_results()
#' target_sets <- list(locus_a = c(1, 5), locus_b = c(10))
#' freqs <- calculate_allele_frequencies_by_sets(result, target_sets)
#' print(freqs$locus_a)
calculate_allele_frequencies_by_sets <- function(results, target_sets, ...) {
  if (!is.list(target_sets)) {
    stop("target_sets must be a list")
  }
  if (length(target_sets) == 0) {
    stop("target_sets must contain at least one set")
  }
  for (i in seq_along(target_sets)) {
    set <- target_sets[[i]]
    if (!is.vector(set) || length(set) == 0) {
      stop("Each element of target_sets must be a non-empty vector of indices or target names")
    }
  }
  out <- lapply(target_sets, function(set) {
    calculate_allele_frequencies(results, snp_indices = set, ...)
  })
  if (!is.null(names(target_sets))) {
    names(out) <- names(target_sets)
  } else {
    names(out) <- paste0("set_", seq_along(target_sets))
  }
  out
}

#' Summarize allele frequencies from per-sample counts (MCMC path).
#' Builds per-sample frequencies, then returns posterior mean, SD, credible interval,
#' and per-sample mean count. Sample-size invariant.
#' @keywords internal
summarize_allele_frequencies_mcmc <- function(allele_counts_list, interval = 0.95, n_samples = length(allele_counts_list)) {
  all_alleles <- unique(unlist(lapply(allele_counts_list, function(x) x$allele)))
  n_samp <- length(allele_counts_list)
  total_per_sample <- vapply(allele_counts_list, function(x) sum(x$count), numeric(1L))

  # Count matrix: rows = alleles, columns = samples
  count_matrix <- matrix(0, nrow = length(all_alleles), ncol = n_samp, dimnames = list(all_alleles, NULL))
  for (j in seq_len(n_samp)) {
    sample_df <- allele_counts_list[[j]]
    for (i in seq_len(nrow(sample_df))) {
      a <- sample_df$allele[i]
      count_matrix[a, j] <- sample_df$count[i]
    }
  }

  # Per-sample frequency (0 where total is 0)
  freq_matrix <- count_matrix
  for (j in seq_len(n_samp)) {
    tot <- total_per_sample[j]
    if (tot > 0) {
      freq_matrix[, j] <- count_matrix[, j] / tot
    }
  }

  probs <- c((1 - interval) / 2, 1 - (1 - interval) / 2)
  frequency <- rowMeans(freq_matrix)
  frequency_sd <- apply(freq_matrix, 1L, stats::sd)
  frequency_lower <- apply(freq_matrix, 1L, stats::quantile, probs = probs[1L], names = FALSE, na.rm = TRUE)
  frequency_upper <- apply(freq_matrix, 1L, stats::quantile, probs = probs[2L], names = FALSE, na.rm = TRUE)
  mean_count <- rowMeans(count_matrix)
  mean_total_parasites <- mean(total_per_sample)

  result_df <- data.frame(
    allele = all_alleles,
    frequency = frequency,
    frequency_sd = frequency_sd,
    frequency_lower = frequency_lower,
    frequency_upper = frequency_upper,
    mean_count = mean_count,
    n_samples = n_samp,
    stringsAsFactors = FALSE
  )
  attr(result_df, "mean_total_parasites") <- mean_total_parasites
  result_df
}

#' Calculate allele counts for a single sample
#' @keywords internal
calculate_allele_counts_single <- function(A, D, snp_indices, r0_values, r1_values, sep = "|") {
  
  n_snps <- length(snp_indices)
  n_strains <- nrow(D)
  n_individuals <- nrow(A)
  rs <- Map(c, r1_values, r0_values)
  
  # Get the SNP values for each strain at the specified indices
  strain_snps <- D[, snp_indices, drop = FALSE]
  
  # Generate all possible allele combinations. Monomorphic loci have r0 == r1, so the
  # cartesian product yields duplicate strings; collapse them so each allele appears once.
  all_combinations <- expand.grid(rs)
  allele_strings <- unique(apply(all_combinations, 1, paste, collapse = sep))

  # Initialize count vector
  allele_counts <- rep(0, length(allele_strings))
  names(allele_counts) <- allele_strings
  
  # For each individual, calculate their allele composition
  for (i in seq_len(n_individuals)) {
    # Get strains present in this individual
    individual_strains <- which(A[i, ] > 0)
    
    if (length(individual_strains) > 0) {
      # For each strain in this individual
      for (strain_idx in individual_strains) {
        # Get the allele for this strain
        strain_allele <- strain_snps[strain_idx, ]
        strain_allele_string <- paste(
          Map(\(snp, variants) variants[snp + 1], strain_allele, rs), 
          collapse = sep
        )
        strain_allele_idx <- match(strain_allele_string, allele_strings)

        # Add the count for this strain to the allele count
        allele_counts[strain_allele_idx] <- allele_counts[strain_allele_idx] + A[i, strain_idx]
      }
    }
  }
  
  # Return as data frame
  return(data.frame(
    allele = names(allele_counts),
    count = allele_counts
  ))
}

#' Aggregate allele counts across multiple MCMC samples
#' @keywords internal
aggregate_allele_counts <- function(allele_counts_list) {
  
  # Get all unique alleles
  all_alleles <- unique(unlist(lapply(allele_counts_list, function(x) x$allele)))
  
  # Initialize aggregated counts
  aggregated_counts <- rep(0, length(all_alleles))
  names(aggregated_counts) <- all_alleles
  
  # Sum counts across all samples
  for (sample_counts in allele_counts_list) {
    for (i in seq_len(nrow(sample_counts))) {
      allele <- sample_counts$allele[i]
      count <- sample_counts$count[i]
      aggregated_counts[allele] <- aggregated_counts[allele] + count
    }
  }
  
  # Return as data frame
  return(data.frame(
    allele = names(aggregated_counts),
    count = aggregated_counts,
    stringsAsFactors = FALSE
  ))
}
