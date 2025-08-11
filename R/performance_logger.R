#' Performance Logger for SNP-Slice
#'
#' @description
#' Utility functions for tracking performance and identifying bottlenecks

# Performance tracking environment
performance_env <- new.env()

#' Initialize performance tracking
#'
#' @keywords internal
init_performance_tracking <- function() {
  performance_env$performance_log <- list()
  performance_env$performance_timers <- list()
}

#' Start timing a function
#'
#' @param function_name Name of the function being timed
#' @keywords internal
start_timer <- function(function_name) {
  if (is.null(performance_env$performance_timers)) {
    init_performance_tracking()
  }
  performance_env$performance_timers[[function_name]] <- Sys.time()
}

#' End timing a function
#'
#' @param function_name Name of the function being timed
#' @keywords internal
end_timer <- function(function_name) {
  if (is.null(performance_env$performance_timers) || is.null(performance_env$performance_timers[[function_name]])) {
    return()
  }
  
  elapsed <- as.numeric(difftime(Sys.time(), performance_env$performance_timers[[function_name]], units = "secs"))
  
  if (is.null(performance_env$performance_log[[function_name]])) {
    performance_env$performance_log[[function_name]] <- list(
      total_time = 0,
      call_count = 0,
      avg_time = 0
    )
  }
  
  performance_env$performance_log[[function_name]]$total_time <- performance_env$performance_log[[function_name]]$total_time + elapsed
  performance_env$performance_log[[function_name]]$call_count <- performance_env$performance_log[[function_name]]$call_count + 1
  performance_env$performance_log[[function_name]]$avg_time <- performance_env$performance_log[[function_name]]$total_time / performance_env$performance_log[[function_name]]$call_count
  
  # Clean up timer
  performance_env$performance_timers[[function_name]] <- NULL
}

#' Get performance summary
#'
#' @return Performance summary data frame
#' @export
get_performance_summary <- function() {
  if (is.null(performance_env$performance_log) || length(performance_env$performance_log) == 0) {
    return(data.frame())
  }
  
  summary_data <- data.frame(
    function_name = names(performance_env$performance_log),
    total_time = sapply(performance_env$performance_log, function(x) x$total_time),
    call_count = sapply(performance_env$performance_log, function(x) x$call_count),
    avg_time = sapply(performance_env$performance_log, function(x) x$avg_time),
    stringsAsFactors = FALSE
  )
  
  # Sort by total time
  summary_data <- summary_data[order(summary_data$total_time, decreasing = TRUE), ]
  
  return(summary_data)
}

#' Clear performance log
#'
#' @export
clear_performance_log <- function() {
  performance_env$performance_log <- NULL
  performance_env$performance_timers <- NULL
}

#' Print performance summary
#'
#' @export
print_performance_summary <- function() {
  summary <- get_performance_summary()
  if (nrow(summary) == 0) {
    cat("No performance data available\n")
    return()
  }
  
  cat("Performance Summary:\n")
  cat("===================\n")
  for (i in 1:nrow(summary)) {
    row <- summary[i, ]
    cat(sprintf("%-25s: %8.3f s total, %6d calls, %8.4f s avg\n", 
                row$function_name, row$total_time, row$call_count, row$avg_time))
  }
  cat("===================\n")
}
