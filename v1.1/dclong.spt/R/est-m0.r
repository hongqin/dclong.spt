#' @title Estimate the Number of True Null Hypotheses
#'
#' @description Estimating the number of true null hypotheses (m0) 
#' is critical for estiamting the false discover rate (FDR). 
#' The functions listed here offers different methods for estimating m0 
#' for both regular pvalues and sequential permutation test
#' pvalues. 
#' 
#' @return the estimate of the number of true null hypotheses.  
#' @seealso \code{\link{fdr}} for estimating false discover rate and calculating qvalues.
#' @references 
#' Storey JD and Tibshirani R. (2003) Statistical significance for genome-wide studies.
#' 
#'
#' Dan Nettleton, J. T. Gene Hwang, Rico A. Caldo and Roger P. Wise (2006) 
#' Estimating the Number of True Null Hypotheses from a Histogram of p Values.
#'
#'
#' Tim Bancroft, Chuanlong Du and Dan Nettleton (2012) Estimation of False Discovery 
#' Rate Using Sequential Permutation P -Values.

#' @examples
#' \dontrun{
#'   #simulate p vlaues 
#'   p = rbunif(10000,beta=29)
#'   #estimate m0 using Storey's method
#'   m0.storey(p)
#'   #estimate m0 using Nettleton's method for regular pvalues
#'   m0.nettleton(p)
#'   #load data 
#'   data(leukemia)
#'   #sequential permutation pvalues
#'   spt.mean(leukemia,5,10,1000)[,1] -> p
#'   #estimate m0 using Nettleton's method for sequential 
#'   #permutation test pvalues
#'   m0.nettleton(p,control=list(h=10,n=1000))
#' }
#' @rdname est-m0
#' @export
#' @param p a vector of pvalues.
#' @param lambda points chosen for fitting spline which can be considered 
#' as tuning parameters in estimating m0. These points must be in [0,1] and 
#' a defualt value seq(0,0.95,0.05) is used. For more information, please 
#' see Storey and Tibshirani's (PNAS, 2003).

m0.storey = function(p,lambda=seq(0,0.95,0.05)){
  m = length(p) 
  lambda.len = length(lambda)
  pi0 = rep(0,lambda.len)
  for(i in 1:lambda.len) {
    pi0[i] <- mean(p>=lambda[i]) / (1-lambda[i])
  }
  spi0 <- smooth.spline(lambda,pi0,df=3)
  pi0 <- max(predict(spi0,x=1)$y,0)
  pi0 <- min(pi0,1)
  return(pi0*m)
} 

#' @rdname est-m0
#' @export
#' @param bins the number of bins in Nettleton's histogram based method for 
#' estimating m0 (see Nettleton et al. (2006) JABES 11, 337-356 
#' and Bancroft et al. (2012)).
#' @param control either \code{NULL} or a list. If \code{control} is NULL, 
#' then Nettleton's method for uniformlly distribution 
#' pvalues under null hypothesis is used. If \code{control} is a list  
#' (must contain at least two variables \code{h} which is the number of significant
#' test statistics to hit before early termination in sequential permutation
#' test) and \code{n} which is the number of permutations sampled),
#' Nettleton's method for non-uniformly pvalues is used. 
#' You can include an extra variable \code{delta} in list \code{control} which 
#' is the error tolerance in counting pvalues. A default value \code{h/(10*n^2)} 
#' is used. 
m0.nettleton = function(p,bins=20,control){
	if(is.null(control)){
    return(m0.nettleton.regular(p=p,B=bins))
	}
  h = control$h
  n = control$n
  if(is.null(h)||is.null(n)){
    stop('argument "control" has an illegal format.')
  }
  delta = control$delta
  if(is.null(delta)){
    delta = h/(10*n^2)
  }
  return(m0.nettleton.sequential(p,bcutoff(p,h,n,delta,bins)))
}

#' @noRd
#' Nettleton's method for unifromly distributed pvalues under null 
#' hypothesis, which is equivalent to mosig et al. Genetics 157: 1683-1698.
#' For mor information, see Nettleton et al. (2006) JABES 11: 337-356.
m0.nettleton.regular = function(p,B)
{
  m <- length(p)
  index = 1:B
  bin <- c(-0.1, index/B)
  bin[B+1] = 2#make the last cut-off one bigger to avoid numerical problems
  bin.counts=rep(0,B)
  for(i in index){
    bin.counts[i]=sum((p>bin[i])&(p<=bin[i+1]))
  }
  tail.means <- rev(cumsum(rev(bin.counts))/index)
  tail.means[min(which(bin.counts - tail.means <= 0))] * B
}

#' @noRd 
#' Nettleton's method for pvalues obtained using sequential permutation
#' test. 
m0.nettleton.sequential = function(p,cutoff){
  probs = cutoff[,3]
  counts = cutoff[,4]
  revCumSum.prob = rev(cumsum(rev(probs)))
  revCumSum.count = rev(cumsum(rev(counts)))
  fraction.prob = probs/revCumSum.prob
  fraction.count = counts/revCumSum.count
  flag = min(which(fraction.count<=fraction.prob))
  revCumSum.count[flag]/revCumSum.prob[flag]
}

#' @noRd
#' Calculate cut-offs for estimating m0 based on non-uniformly distributed
#' pvalues under null hypotheses. The cut-offs are decided based on the 
#' rule of making approximately equl probability for each bin. 
#' For more detains, see Bancroft et al. (2012). 
bcutoff = function(p,h,n,delta,bins){
  prob = 1/bins
  pvals.support(h,n) -> pvals.theory
  pvals = pvals.theory[,1]
  cutoff = NULL
  currentCut = 0
  while(currentCut<1-delta){
    index = which(pvals>currentCut+prob-delta)
    if(length(index)==0){
      currentCut = 1
    }else{
      currentCut = min(pvals[index])
    }
    if(currentCut<1-delta){
      cutoff = c(cutoff,currentCut)
    }else{
      cutoff = c(cutoff,1)
    }
  }
  cutoff.length = length(cutoff)
  cutoff = cbind(rep(0,cutoff.length),cutoff,rep(0,cutoff.length),rep(0,cutoff.length))
  for(i in 2:cutoff.length){
    cutoff[i,1] = cutoff[i-1,2]
  }
  for(i in 1:cutoff.length){
    cutoff[i,3] = cutoff[i,2] - cutoff[i,1]
    cutoff[i,4] = sum(p>cutoff[i,1]+delta&p<cutoff[i,2]+delta)
  }
  colnames(cutoff) = c("lower","upper","prob","count")
  cutoff
}

#' @noRd
#' Calculate theorectical sequential pvalues and corresponding probabilities.
pvals.support= function(h,n){
  pvals = c((1:(h-1))/n,h/(n:h))
  prob = diff(c(0,pvals))
  cbind(pvals,prob)
}