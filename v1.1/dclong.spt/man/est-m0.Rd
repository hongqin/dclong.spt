\name{m0.storey}
\alias{m0.nettleton}
\alias{m0.storey}
\title{Estimate the Number of True Null Hypotheses}
\usage{
  m0.storey(p, lambda = seq(0, 0.95, 0.05))

  m0.nettleton(p, bins = 20, control)
}
\arguments{
  \item{p}{a vector of pvalues.}

  \item{lambda}{points chosen for fitting spline which can
  be considered as tuning parameters in estimating m0.
  These points must be in [0,1] and a defualt value
  seq(0,0.95,0.05) is used. For more information, please
  see Storey and Tibshirani's (PNAS, 2003).}

  \item{bins}{the number of bins in Nettleton's histogram
  based method for estimating m0 (see Nettleton et al.
  (2006) JABES 11, 337-356 and Bancroft et al. (2012)).}

  \item{control}{either \code{NULL} or a list. If
  \code{control} is NULL, then Nettleton's method for
  uniformlly distribution pvalues under null hypothesis is
  used. If \code{control} is a list (must contain at least
  two variables \code{h} which is the number of significant
  test statistics to hit before early termination in
  sequential permutation test) and \code{n} which is the
  number of permutations sampled), Nettleton's method for
  non-uniformly pvalues is used. You can include an extra
  variable \code{delta} in list \code{control} which is the
  error tolerance in counting pvalues. A default value
  \code{h/(10*n^2)} is used.}
}
\value{
  the estimate of the number of true null hypotheses.
}
\description{
  Estimating the number of true null hypotheses (m0) is
  critical for estiamting the false discover rate (FDR).
  The functions listed here offers different methods for
  estimating m0 for both regular pvalues and sequential
  permutation test pvalues.
}
\examples{
\dontrun{
  #simulate p vlaues
  p = rbunif(10000,beta=29)
  #estimate m0 using Storey's method
  m0.storey(p)
  #estimate m0 using Nettleton's method for regular pvalues
  m0.nettleton(p)
  #load data
  data(leukemia)
  #sequential permutation pvalues
  spt.mean(leukemia,5,10,1000)[,1] -> p
  #estimate m0 using Nettleton's method for sequential
  #permutation test pvalues
  m0.nettleton(p,control=list(h=10,n=1000))
}
}
\references{
  Storey JD and Tibshirani R. (2003) Statistical
  significance for genome-wide studies.

  Dan Nettleton, J. T. Gene Hwang, Rico A. Caldo and Roger
  P. Wise (2006) Estimating the Number of True Null
  Hypotheses from a Histogram of p Values.

  Tim Bancroft, Chuanlong Du and Dan Nettleton (2012)
  Estimation of False Discovery Rate Using Sequential
  Permutation P -Values.
}
\seealso{
  \code{\link{fdr}} for estimating false discover rate and
  calculating qvalues.
}

