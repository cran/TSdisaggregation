#' Refit LASSO estimate into GLS 
#' 
#' This function reduces the bias in LASSO estimates by re-fitting the
#' support back into GLS. 
#' 
#' @param X     Aggregated indicator series matrix that has been GLS rotated. 
#' @param Y     Low-frequency response vector that has been GLS rotated. 
#' @param beta  Estimate of beta from LARS algorithm for a certain lambda.
#' @keywords internal 
#' @importFrom Rdpack reprompt
#' @importFrom stats lm rbinom rnorm

refit <- function(X, Y, beta) {
  
  p = ncol(X)
  active <- which(beta != 0 )
  X_active <- X[,active]
  lm_fit <- lm(Y ~ 0 + X_active)
  betahat_ols <- lm_fit$coefficients
  betahat_debias <- rep(0,p)
  betahat_debias[active] <- betahat_ols
  
  return(betahat_debias)
  
}