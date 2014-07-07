#' @title Mixture of Beta and Uniform Distribution
#' @description Generate random observations from a mixture of 
#' beta and uniform distribution. It's primarily used for demostrating 
#' use of other functions in this package. 
#' @param n number of observations to generate.
#' @param alpha the first parameter of beta distribution.
#' @param beta the second parameter of beta distribution.
#' @param gamma probability of a observation coming from uniform distribution.
#' @export
#' @details The mixture distribution is gamma*U(0,1) + (1-gamma)*Beta(alpha,beta).
#' @examples
#' \dontrun{
#'    rbunif(100,alpha=1,beta=29,gamma=0.7)
#' }
rbunif = function(n,alpha,beta,gamma){
  sam = rep(0,n)
  for(i in 1:n){
    if(rbinom(1,1,gamma)){
      sam[i] = runif(1,0,1) 
    }else{
      sam[i] = rbeta(1,alpha,beta)
    }
  }
  sam
}
