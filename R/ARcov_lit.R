#' Function to generate an ARIMA(1,1,0) variance-covariance matrix for the Litterman method with parameter rho s.t. |rho| < 1. 
#'  
#' @param rho
#' @param n
#' @keywords internal 
#' @importFrom Rdpack reprompt
#' @importFrom stats lm rbinom rnorm


ARcov_lit <- function(rho, n) {
  
  diags <- list(rep(1, times = n), rep(-1, times = n-1))
  Delta_t <- bandSparse(n, k = 0:1, diagonals = diags, symmetric = FALSE)
  Delta <- t(Delta_t) 
  
  diags2 <- list(rep(1, times = n), rep(-rho, times = n-1))
  H_r_t <-  bandSparse(n, k = 0:1, diagonals = diags2, symmetric = FALSE)
  H_r <- t(H_r_t) 
  
  Sigma <- solve(Delta_t %*% H_r_t %*% H_r %*% Delta)
  
  Sigma <- forceSymmetric(Sigma)
  
  return(Sigma)
  
}


