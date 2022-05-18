#' Function to do Chow-Lin temporal disaggregation from \insertCite{chow1971best;textual}{TSdisaggregation} and Litterman. 
#' 
#' Used in disaggregation.R to find estimates given the optimal rho parameter. 
#' 
#' @param Y  		The low-frequency response series (n_l x 1 matrix).
#' @param X  		The high-frequency indicator series (n x p matrix).
#' @param rho   The AR(1) residual parameter (strictly between -1 and 1).
#' @param aggMat 	Aggregation matrix according to 'first', 'sum', 'average', 'last' (default is 'sum').
#' @param aggRatio Aggregation ratio e.g. 4 for annual-to-quarterly, 3 for quarterly-to-monthly (default is 4). 
#' @param litterman TRUE to use litterman vcov. FALSE for Chow-Lin vcov. Default is FALSE. 
#' @return y	Estimated high-frequency response series (n x 1 matrix).
#' @return betaHat	Estimated coefficient vector (p x 1 matrix).
#' @return u_l	Estimated aggregate residual series (n_l x 1 matrix). 
#' @keywords chow lin litterman temporal disaggregation
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt	
#' @importFrom stats lm rbinom rnorm


chowlin <- function(Y, X, rho, aggMat, aggRatio, litterman = FALSE) {
  
  n_l = dim(Y)[1]
  n = dim(X)[1]
  p = dim(X)[2]
  nfull = aggRatio*n_l
  extr = n - nfull # number of extrapolations
  
  
  # Generate the aggregation matrix C
  
  if(aggMat == 'sum'){
    
    C <- kronecker(diag(n_l), matrix(data = 1, nrow = 1, ncol = aggRatio))
    C <- cbind(C, matrix(0L, n_l, extr))
    
  }else if(aggMat == 'avg'){
    
    C <- kronecker(diag(n_l), matrix(data = 1/aggRatio, nrow = 1, ncol = aggRatio))
    C <- cbind(C, matrix(0L, n_l, extr))
    
  }else if(aggMat == 'first'){
    
    C <- kronecker(diag(n_l), matrix(data = c(1, rep(0, times = aggRatio-1)), nrow = 1, ncol = aggRatio))
    C <- cbind(C, matrix(0L, n_l, extr))
    
  }else if(aggMat == 'last'){
    
    C <- kronecker(diag(n_l), matrix(data = c(rep(0, times = aggRatio-1), 1), nrow = 1, ncol = aggRatio))
    C <- cbind(C, matrix(0L, n_l, extr))
    
  }
  
  X_l = C %*% X
  
  if(litterman) {
    vcov = ARcov_lit(rho, n)
  }else {
    vcov = ARcov(rho, n)
  }
  
  # Simplification and Cholesky factorization of the Sigma 
  
  vcov_agg = forceSymmetric(C %*% vcov %*% t(C))
  Uchol <- chol(vcov_agg)
  Lchol <- t(Uchol)
  
  # Preconditioning the variables
  
  X_F <- solve(Lchol) %*% X_l
  Y_F <- solve(Lchol) %*% Y
  
  # Estimate betaHat_0 using GLS assuming Sigma with rho
  
  betaHat <- solve(t(X_F) %*% X_F) %*% t(X_F) %*% Y_F 

  
  # The distribution matrix
  
  D <- vcov %*% t(C) %*% solve(vcov_agg)
  
  # Obtain the residuals using betaHat_1
  
  u_l <- Y - X_l %*% betaHat  
  
  # Generate the high-frequency series
  
  y <- X %*% betaHat + (D %*% u_l)
  
  output = list('y' = y, 'betaHat' = betaHat, 'u_l' = u_l)
  

  return(output)
  
}