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
  
  if (is.matrix(data) || is.data.frame(data)) {
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

#' Preprocess data for SNP-Slice
#'
#' @param data Input data
#' @param model Model type
#'
#' @return Processed data list with y, r, and metadata
#' @keywords internal
preprocess_data <- function(data, model) {
  
  if (is.character(data)) {
    # Read from file
    data <- read_snp_data(data)
  }
  
  if (is.matrix(data) || is.data.frame(data)) {
    # Direct matrix input
    y <- as.matrix(data)
    
    if (model == "categorical") {
      # Categorical data - no r matrix needed
      return(list(
        y = y,
        r = NULL,
        N = nrow(y),
        P = ncol(y),
        data_type = "categorical",
        model = model
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
    
    return(list(
      y = y,
      r = r,
      N = nrow(y),
      P = ncol(y),
      data_type = "read_counts",
      model = model
    ))
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
      N = processed_data$N,
      P = processed_data$P,
      data_type = processed_data$data_type
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
  results_file <- system.file("extdata", "example_analysis_results.rds", package = "snp.slice")
  
  if (!file.exists(results_file)) {
    warning("Example analysis results not found. Running analysis now (this may take a few minutes)...")
    
    # Load example data
    data(example_snp_data, package = "snp.slice")
    data <- list(read0 = example_snp_data$read0, read1 = example_snp_data$read1)
    
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
