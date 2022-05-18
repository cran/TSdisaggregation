#' Function to calculate the BIC score from sparse temporal disaggregation.  
#' 
#' Used in disaggregation.R to find estimates of the optimal rho parameter. 
#' 
#' @param Y  		The low-frequency response series (n_l x 1 matrix).
#' @param X  		The aggregated high-frequency indicator series (n_l x p matrix).
#' @param vcov Aggregated variance-covariance matrix of AR(1) residuals. 
#' @keywords chow lin litterman temporal disaggregation
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt	
#' @importFrom stats lm rbinom rnorm


sptd_BIC <- function(Y,X,vcov) {
  
  n_l = dim(Y)[1]
  
  # Simplification and Cholesky factorization of the Sigma 
  
  Uchol <- chol(vcov)
  Lchol <- t(Uchol)
  
  # Preconditioning the variables
  
  X_F <- solve(Lchol) %*% X
  Y_F <- solve(Lchol) %*% Y
  
  
  # Fit LARS algorithm to the data 
  lars.fit <- lars(X_F, Y_F, intercept = F, normalize = F)
  betamat <- lars.fit$beta 
  
  # Don't allow support to be bigger than n_l/2
  npath <- k.index(betamat, n_l)
  
  # Find BIC for each re-fitted betahat 
  beta_refit <- list()
  BIC <- c()
  BIC[1] <- hdBIC(X_F, Y_F, vcov, betamat[1,])
  beta_refit[[1]] <- betamat[1,]
  
  for(lam in 2:npath) {
    
    beta_refit[[lam]] <- refit(X_F, Y_F, betamat[lam, ])
    BIC[lam] <- hdBIC(X_F, Y_F, vcov, beta_refit[[lam]])
  }
  
  return(min(BIC))
  
}