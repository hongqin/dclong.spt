\name{fdr}
\alias{fdr}
\title{False Discover Rate (FDR)}
\usage{
  fdr(cutoff, p, m0.hat, delta)
}
\arguments{
  \item{cutoff}{a vector of cut-offs. For a given cut-off,
  the null hypotheses with pvalues less or equal to the
  cut-off are rejected.}

  \item{p}{the vector of pvalues.}

  \item{m0.hat}{an estimate of the number of true null
  hypotheses.}

  \item{delta}{the error tolerance in comparing pvalues
  with cut-offs. To specify an appropriate error tolerance
  is important if there are lots of pvalues that are very
  close to a cut-off. An extreme but happen-often (e.g.
  when calculating qvalues) case is that a cut-off is
  inside the support of sequential permutation pvalues. A
  value no greater than \code{h/(100*n^2)} is recommended
  in this case. Please refer to \code{\link{spt}} for
  documentation of \code{h} and \code{n}.}
}
\description{
  Function \code{fdr} calculate FDR for specified cut-offs.
  If the pvalues is specified as the cutoffs, then you get
  the qvalues.
}

