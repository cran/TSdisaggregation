#' BIC score 
#' 
#' This function calculates the BIC score that has been shown to work better than ordinary BIC in 
#' high-dimensional scenarios. It uses the variance estimator given in \insertCite{yu2019estimating;textual}{TSdisaggregation}.
#' 
#' @param X           Aggregated indicator series matrix that has been GLS rotated. 
#' @param Y           Low-frequency response vector that has been GLS rotated. 
#' @param covariance  Aggregated AR covariance matrix. 
#' @param beta        Estimate of beta from LARS algorithm for a certain lambda.
#' @keywords internal 
#' @references 
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @importFrom stats lm rbinom rnorm

hdBIC <- function(X, Y, covariance, beta) {
  
  n_l <- length(Y)
  support <- sum(beta != 0)
  u <- Y-X%*%beta
  log.lik <- -(n_l-support)/2 - n_l/2*log(2*pi/(n_l-support)*crossprod(u)) - log(det(covariance))/2
  BIC <- -2*(log.lik) + log(n_l)*support
  
  return(BIC)
  
}