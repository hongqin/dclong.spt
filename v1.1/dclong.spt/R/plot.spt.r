#' @title Plot Sequential Permutation Test Pvalues
#'
#' @description The function \code{plot.spt} plot the observed sequential permutation pvalues
#' together with the theoretical ones when all null hypotheses are true, 
#' so that you can have a general idea about whether there
#' are lots of significant tests. It's analogous to the usual histogram of regular 
#' pvalues.
#'
#' @details If the plot shows a lack of big pvalues but overabundance of small pvalues,
#' one expect there to be some siginificant tests. 
#' 
#' @usage \method{plot}{spt}(x,plim=1,...)
#'
#' @param x on object of class ``spt''.
#' @param plim the upper limit of pvalues to display. 
#' @param ... some extra parameters than can be passed to function \code{plot}.
#'
#' @author Dan Nettleton and Chuanlong Du.
#' @export
#' @examples
#' \dontrun{
#'    #load data
#'    data(leukemia)
#'    spt.mean(leukemia,5,10,1000) -> spt.mean.out
#'    plot.spt(spt.mean.out,col="red")
#' }
#' @S3method plot spt
plot.spt=function(x,plim=1,...)
{
  p = x$p
  h = x$h
  n = x$n
  #-----------------------------------
  s = c((1:h)/n,h/((n-1):h))
  tp = table(p)
  counts = rep(0,n)
  ps = as.numeric(names(tp))
  for(i in 1:length(ps)){
    x = which.min(abs(ps[i]-s))
    counts[x] = tp[i] 
  }
  probs = counts/length(p)
  y = diff(c(0,s))
  plot(s[s<=plim],probs[s<=plim],pch=1,type="l",
       ylim=c(0,max(c(probs[s<=plim]),y[s<=plim])),xlab="p-value",ylab="Probability",...)  
  for(i in 1:length(s)){
    lines(rep((s[s<=plim])[i],2),c(0,(y[s<=plim])[i]),lwd=1,col=1,lty=1)
  }
  out=cbind(s,probs)
  invisible(out)
}
