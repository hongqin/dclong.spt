#' @title False Discover Rate (FDR)
#'
#' @description Function \code{fdr} calculate FDR for specified cut-offs. 
#' If the pvalues is specified as the cutoffs, then you get the qvalues.
#'
#' @param cutoff a vector of cut-offs. 
#' For a given cut-off, the null 
#' hypotheses with pvalues less or equal to the cut-off are rejected.
#'
#' @param p the vector of pvalues. 
#' 
#' @param m0.hat an estimate of the number of true null hypotheses. 
#' 
#' @param delta the error tolerance in comparing pvalues with cut-offs. 
#' To specify an appropriate error tolerance is important if there are lots of 
#' pvalues that are very close to a cut-off. 
#' An extreme but happen-often (e.g. when calculating qvalues) case
#' is that a cut-off is inside the support of sequential permutation pvalues. 
#' A value no greater than \code{h/(100*n^2)} is recommended in this case. 
#' Please refer to \code{\link{spt}} for documentation of \code{h} and \code{n}.
#' @export
fdr = function(cutoff, p, m0.hat, delta){
  p.sorted = sort(p)
  p.sorted*m0.hat/(1:length(p.sorted)) -> p.adjusted
  cutoff.length = length(cutoff)
  result = rep(0,cutoff.length)
  for(i in 1:cutoff.length){
    index = which(p.sorted>cutoff[i]-delta)
    if(length(index)==0){
      result[i] = 1
    }else{
      result[i] = min(p.adjusted[index])      
    }
  }
  result
}


