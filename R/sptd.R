#' Function to do sparse temporal disaggregation from \insertCite{mosley2021sparse;textual}{TSdisaggregation}. 
#' 
#' Used in disaggregation.R to find estimates given the optimal rho parameter. 
#' 
#' @param Y  		The low-frequency response series (n_l x 1 matrix).
#' @param X  		The high-frequency indicator series (n x p matrix).
#' @param rho   The AR(1) residual parameter (strictly between -1 and 1).
#' @param aggMat 	Aggregation matrix according to 'first', 'sum', 'average', 'last' (default is 'sum').
#' @param aggRatio Aggregation ratio e.g. 4 for annual-to-quarterly, 3 for quarterly-to-monthly (default is 4). 
#' @param adaptive TRUE to use adaptive lasso penalty. FALSE for lasso penalty. Default is FALSE. 
#' @return y	Estimated high-frequency response series (n x 1 matrix).
#' @return betaHat	Estimated coefficient vector (p x 1 matrix).
#' @return u_l	Estimated aggregate residual series (n_l x 1 matrix). 
#' @keywords sparse lasso temporal disaggregation
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt	
#' @importFrom stats lm rbinom rnorm



sptd <- function(Y, X, rho, aggMat, aggRatio, adaptive = FALSE) {
  
  
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
  
  # Simplification and Cholesky factorization of the Sigma 
  
  vcov = ARcov(rho, n)
  vcov_agg = C %*% vcov %*% t(C)
  
  Uchol <- chol(vcov_agg)
  Lchol <- t(Uchol)
  
  # Preconditioning the variables
  
  X_F <- solve(Lchol) %*% X_l
  Y_F <- solve(Lchol) %*% Y
  
  
  # Fit LARS algorithm to the data 
  lars.fit <- lars(X_F, Y_F, intercept = F, normalize = F)
  betamat <- lars.fit$beta 
  
  # Don't allow support to be bigger than n_l/2
  npath <- k.index(betamat, n_l)
  
  # Find BIC for each re-fitted betahat 
  beta_refit <- list()
  BIC <- c()
  BIC[1] <- hdBIC(X_F, Y_F, vcov_agg, betamat[1,])
  beta_refit[[1]] <- betamat[1,]
  
  for(lam in 2:npath) {
    
    beta_refit[[lam]] <- refit(X_F, Y_F, betamat[lam, ])
    BIC[lam] <- hdBIC(X_F, Y_F, vcov_agg, beta_refit[[lam]])
  }
  
  min_bic_idx <- which.min(BIC)
  betaHat <- beta_refit[[min_bic_idx]]
  
  if(adaptive) {
    
    # get adaptive lasso weight
    ada_weight <- abs(betaHat)
    
    # Scale the data A*|betahat_lasso|
    X_F_scaled <- scale(X_F, center = F, scale = 1/ada_weight)
    
    # Apply lars to scaled data 
    lars.fit <- lars(X_F_scaled, Y_F, intercept = F, normalize = F)
    betamat <- lars.fit$beta 
    
    # Scale beta estimates back to original 
    betamat_scaled <- scale(betamat, center = F, scale = 1/ada_weight)
    
    # Tune with BIC 
    npath <- k.index(betamat_scaled, n_l)
    beta_refit <- list()
    BIC <- c()
    BIC[1] <- hdBIC(X_F, Y_F, vcov_agg, betamat_scaled[1,])
    beta_refit[[1]] <-  betamat_scaled[1,]
    
    for(lam in 2:npath) {
      
      beta_refit[[lam]] <- refit(X_F, Y_F, betamat_scaled[lam, ])
      BIC[lam] <- hdBIC(X_F, Y_F, vcov_agg, beta_refit[[lam]])
      
    }
    
    # Store best BIC, betahat and lambdahat 
    min_bic_idx <- which.min(BIC)
    betaHat <- beta_refit[[min_bic_idx]]
  }
  
  # The distribution matrix
  
  D <- vcov %*% t(C) %*% solve(vcov_agg)
  
  # Obtain the residuals using betaHat_1
  
  u_l <- Y - X_l %*% betaHat  
  
  # Generate the high-frequency series
  
  y <- X %*% betaHat + (D %*% u_l)
  
  output = list('y' = y, 'betaHat' = betaHat, 'u_l' = u_l)
  
  return(output)
  
}