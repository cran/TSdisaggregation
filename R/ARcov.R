#' Function to generate an AR(1) variance-covariance matrix with parameter rho s.t. |rho| < 1. 
#'  
#' @param rho
#' @param n
#' @keywords internal 
#' @importFrom Rdpack reprompt
#' @importFrom stats lm rbinom rnorm


ARcov <- function(rho,n) {
  
  sqnc <- rho^seq(0, n, by = 1)
  Omega <- toeplitz(sqnc[1: n])
  sig <- (1/(1-rho^2)) * Omega
  return(sig)
  
}
