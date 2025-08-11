#' Log prior for allocation matrix A
#'
#' @param A Allocation matrix
#' @param mu Stick-breaking weights
#' @param alpha IBP concentration parameter
#'
#' @return Log prior probability
#' @keywords internal
logprior_a <- function(A, mu, alpha) {
  start_timer("logprior_a")
  
  # IBP prior for allocation matrix
  lpr <- sum(t(A) * log(mu))
  lpr <- lpr + sum(t(1 - A) * log(1 - mu))
  
  end_timer("logprior_a")
  return(lpr)
}

#' Log prior for stick-breaking weights mu
#'
#' @param mu Stick-breaking weights
#' @param alpha IBP concentration parameter
#' @param N Number of hosts
#'
#' @return Log prior probability
#' @keywords internal
logprior_mu <- function(mu, alpha, N) {
  start_timer("logprior_mu")
  
  # Stick-breaking prior for mu
  lpr <- length(mu) * log(alpha) + alpha * log(tail(mu, 1)) - sum(log(mu))
  
  end_timer("logprior_mu")
  return(lpr)
}

#' Log prior for dictionary matrix D
#'
#' @param D Dictionary matrix
#' @param rho Dictionary sparsity parameter
#'
#' @return Log prior probability
#' @keywords internal
logprior_d <- function(D, rho) {
  start_timer("logprior_d")
  
  # Independent Bernoulli prior for dictionary
  sumD <- sum(D)
  lpr <- sumD * log(rho) + (prod(dim(D)) - sumD) * log(1 - rho)
  
  end_timer("logprior_d")
  return(lpr)
}
