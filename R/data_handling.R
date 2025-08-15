# Global variables for data loading
utils::globalVariables("example_snp_data")

#' Validate input data for SNP-Slice
#'
#' @param data Input data to validate
#' @param model Model type to validate against
#'
#' @return TRUE if valid, otherwise throws error
#' @keywords internal
validate_input_data <- function(data, model) {
  
  # Check if data is NULL or empty
  if (is.null(data) || (is.list(data) && length(data) == 0)) {
    stop("Data cannot be NULL or empty")
  }
  
  # Handle different data types
  if (is.character(data)) {
    # File path - validate file exists
    if (!file.exists(data)) {
      stop("File not found: ", data)
    }
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

  if(is.data.frame(data)) {
    return(TRUE)
    # TODO: validate data frame
  }
  
  stop("Unsupported data type. Expected matrix, data.frame, list, or file path")
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
  
  if (!is.numeric(rho) || rho < 0 || rho > 1) {
    stop("rho must be a number between 0 and 1")
  }
  
  if (!is.numeric(threshold) || threshold < 0 || threshold > 1) {
    stop("threshold must be a number between 0 and 1")
  }
  
  return(TRUE)
}

#' Validate MCMC settings
#'
#' @param n_mcmc Number of MCMC iterations
#' @param burnin Burn-in period
#' @param gap Early stopping threshold
#'
#' @return TRUE if valid, otherwise throws error
#' @keywords internal
validate_mcmc_settings <- function(n_mcmc, burnin, gap) {
  
  if (!is.numeric(n_mcmc) || n_mcmc <= 0 || n_mcmc != as.integer(n_mcmc)) {
    stop("n_mcmc must be a positive integer")
  }
  
  if (!is.null(burnin)) {
    if (!is.numeric(burnin) || burnin < 0 || burnin != as.integer(burnin)) {
      stop("burnin must be a non-negative integer")
    }
    if (burnin >= n_mcmc) {
      stop("burnin must be less than n_mcmc")
    }
  }
  
  if (!is.null(gap)) {
    if (!is.numeric(gap) || gap <= 0 || gap != as.integer(gap)) {
      stop("gap must be a positive integer")
    }
  }
  
  return(TRUE)
}

#' Load data from a dataframe
#'
#' @param data Input dataframe with columns:
#' \describe{
#'   \item{specimen_id}{Specimen ID}
#'   \item{target_id}{Target ID}
#'   \item{target_value}{Target value}
#'   \item{target_count}{Target count (optional)}
#' }
#' For each target, there are at exactly 2 taget values observed. If there is only one, 
#' the second value is set to unknown_target_value. Target count is required if model is
#' not "categorical".
#' @param model Model type
#' @param unknown_target_value Value to use for unknown targets
#' @param target_id_col Name of the target ID column
#' @param target_value_col Name of the target value column
#' @param specimen_id_col Name of the specimen ID column
#' @param target_count_col Name of the target count column
#'
#' @return Processed data list with y, r, and metadata
#' @keywords internal
load_dataframe <- function(
  data, model, unknown_target_value = "?", target_id_col = "target_id", target_value_col = "target_value", specimen_id_col = "specimen_id", target_count_col = "target_count") {
  # Data needs to be completed so that every target has 2 values observed and every specimen has 
  # a value for every target. If only one target value is observed for a specimen, the second value
  # is set to 0. If neither target value is observed for a specimen, both values are set to NA.
 data_renamed <- data |>
    dplyr::rename(
      target_id = !!target_id_col, 
      target_value = !!target_value_col, 
      specimen_id = !!specimen_id_col, 
      target_count = !!target_count_col
    ) |>
    dplyr::select(target_id, target_value, specimen_id, target_count)
  
  complete_targets <- data_renamed |>
    dplyr::select(target_id, target_value) |>
    dplyr::distinct() |>
    dplyr::group_by(target_id) |>
    dplyr::summarize(target_value = list(target_value)) |>
    # remove targets with more than 2 variants
    dplyr::filter(sapply(target_value, length) <= 2) |>
    dplyr::mutate(
      target_value = ifelse(
        sapply(target_value, length) == 1, 
        sapply(target_value, function(x) c(x, unknown_target_value)),
        target_value
      ),
    ) |>
    tidyr::unnest(target_value)
  
  
  completed_data <- data_renamed |>
    dplyr::right_join(complete_targets) |>
    tidyr::complete(
      specimen_id, tidyr::nesting(target_id, target_value), 
      fill = list(target_count = 0)
    ) |>
    dplyr::filter(!is.na(specimen_id)) |>
    dplyr::group_by(specimen_id, target_id) |>
    dplyr::mutate(total_count = sum(target_count)) |>
    dplyr::ungroup() |>
    dplyr::mutate(target_count = ifelse(total_count == 0, NA, target_count)) |>
    dplyr::select(-total_count) |>
    dplyr::group_by(specimen_id, target_id) |>
    dplyr::arrange(dplyr::desc(target_value)) |>
    dplyr::mutate(target_idx = seq(2)) |>
    dplyr::ungroup() |>
    dplyr::arrange(specimen_id, target_id)

  read0_df <- completed_data |>
    dplyr::filter(target_idx == 1) |>
    dplyr::arrange(specimen_id, target_id)

  read1_df <- completed_data |>
    dplyr::filter(target_idx == 2) |>
    dplyr::arrange(specimen_id, target_id)
  
  specimen_ids <- read0_df$specimen_id |> unique()
  target_ids <- read0_df$target_id |> unique()

  nspecs <- length(specimen_ids)
  ntars <- length(target_ids)

  read0_mat <- matrix(read0_df$target_count, nrow = nspecs, ncol = ntars, dimnames = list(specimen_ids, target_ids), byrow = TRUE)
  read1_mat <- matrix(read1_df$target_count, nrow = nspecs, ncol = ntars, dimnames = list(specimen_ids, target_ids), byrow = TRUE)

  # get the target_values in order of the columns of read0_mat
  r0_values_ordered <- read0_df$target_value[match(colnames(read0_mat), read0_df$target_id)]

  # get the target_values in order of the columns of read1_mat
  r1_values_ordered <- read1_df$target_value[match(colnames(read1_mat), read1_df$target_id)]


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


#' Preprocess data for SNP-Slice
#'
#' @param data Input data
#' @param model Model type
#'
#' @return Processed data list with y, r, and metadata
#' @keywords internal
preprocess_data <- function(data, model, ...) {

  if (is.data.frame(data)) {
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
        N = nrow(y),
        P = ncol(y),
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
  
  # Extract MAP estimates
  map_state <- mcmc_result$map_state
  
  # Get active strains only
  active_strains <- colSums(map_state$A) > 0
  allocation_matrix <- map_state$A[, active_strains, drop = FALSE]
  dictionary_matrix <- map_state$D[active_strains, , drop = FALSE]
  
  # Create results object
  results <- list(
    allocation_matrix = allocation_matrix,
    dictionary_matrix = dictionary_matrix,
    mcmc_samples = mcmc_result$samples,
    diagnostics = mcmc_result$diagnostics,
    parameters = mcmc_result$parameters,
    model_info = list(
      model = model_obj$name,
      processed_data = processed_data
    ),
    convergence = mcmc_result$convergence
  )
  
  class(results) <- "snp_slice_results"
  
  return(results)
}

#' Load Example Analysis Results
#'
#' Loads pre-computed SNP-Slice analysis results from the example data.
#' These results were generated using the negative binomial model with 200 MCMC iterations.
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
                       n_mcmc = 200,
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
#' @param allele_sep Separator for allele strings (default: "|")
#'
#' @return A data frame with columns:
#'   \item{allele}{String representation of the allele (e.g., "A|T|T" for 3 SNPs)}
#'   \item{frequency}{Proportion of total parasites with this allele}
#'   \item{count}{Absolute count of parasites with this allele}
#'   \item{total_parasites}{Total number of parasites across all individuals}
#'
#' @export
#' @examples
#' # Load example results
#' result <- load_example_results()
#' 
#' # Calculate allele frequencies for SNPs 1, 5, and 10
#' allele_freqs <- calculate_allele_frequencies(result, c(1, 5, 10))
#' print(allele_freqs)
calculate_allele_frequencies <- function(results, snp_indices, use_map = TRUE, n_samples = 100, allele_sep = "|") {
  if (length(unique(snp_indices)) != length(snp_indices)) {
    stop("snp_indices must be unique")
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
  
  # Get the number of SNPs in the allele
  n_snps <- length(snp_indices)
  r0_values <- results$model_info$processed_data$r0_values[snp_indices]
  r1_values <- results$model_info$processed_data$r1_values[snp_indices] 

  if (use_map) {
    # Use MAP estimates
    A <- results$allocation_matrix
    D <- results$dictionary_matrix
       
    # Calculate allele frequencies for MAP estimate
    allele_counts <- calculate_allele_counts_single(A, D, snp_indices, r0_values, r1_values, sep = allele_sep)
    
  } else {
    # Sample from MCMC results
    if (is.null(results$mcmc_samples)) {
      stop("MCMC samples not available. Set use_map = TRUE or run snp_slice with store_mcmc = TRUE")
    }
    
    # Sample from MCMC results
    n_total_samples <- length(results$mcmc_samples)
    if (n_samples > n_total_samples) {
      warning("Requested ", n_samples, " samples but only ", n_total_samples, " available. Using all available samples.")
      n_samples <- n_total_samples
    }
    
    # Randomly sample from MCMC results
    sample_indices <- sample(seq_len(n_total_samples), n_samples, replace = FALSE)
    
    # Calculate allele frequencies across samples
    allele_counts_list <- lapply(sample_indices, function(i) {
      sample_data <- results$mcmc_samples[[i]]
      calculate_allele_counts_single(sample_data$A, sample_data$D, snp_indices, r0_values, r1_values, sep = allele_sep)
    })
    
    # Aggregate counts across samples
    allele_counts <- aggregate_allele_counts(allele_counts_list)
  }
  
  # Convert to frequency data frame
  total_parasites <- sum(allele_counts$count)
  
  result_df <- data.frame(
    allele = allele_counts$allele,
    frequency = allele_counts$count / total_parasites,
    count = allele_counts$count,
    total_parasites = total_parasites,
    stringsAsFactors = FALSE
  )
  
  # Sort by frequency (descending)
  result_df <- result_df[order(result_df$frequency, decreasing = TRUE), ]
  
  return(result_df)
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
  
  # Generate all possible allele combinations
  all_combinations <- expand.grid(rs)
  allele_strings <- apply(all_combinations, 1, paste, collapse = sep)
  
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
        strain_allele_idx <- which(allele_strings == strain_allele_string)
        
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
