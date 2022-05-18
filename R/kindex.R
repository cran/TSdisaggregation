#' Index of support for LARS algorithm when in high-dimensions 
#' 
#' This function prevents the support of beta becoming greater than n_l/2.
#' This heuristic approach prevents erratic values of BIC when in high-dimensions. 
#' 
#' @param matrix
#' @param n_l
#' @keywords internal 
#' @importFrom Rdpack reprompt
#' @importFrom stats lm rbinom rnorm

k.index <- function(matrix, n_l) {
  
  count <- apply(matrix, 1, function(x) {sum(x != 0)})
  if(max(count) > n_l/2) {
    kindex <- min(which(count > n_l/2))
  }
  else {
    kindex <- nrow(matrix)
  }
  
  return(kindex)
}