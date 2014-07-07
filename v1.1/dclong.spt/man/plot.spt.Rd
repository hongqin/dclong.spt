\name{plot.spt}
\alias{plot.spt}
\title{Plot Sequential Permutation Test Pvalues}
\usage{
  \method{plot}{spt}(x,plim=1,...)
}
\arguments{
  \item{x}{on object of class ``spt''.}

  \item{plim}{the upper limit of pvalues to display.}

  \item{...}{some extra parameters than can be passed to
  function \code{plot}.}
}
\description{
  The function \code{plot.spt} plot the observed sequential
  permutation pvalues together with the theoretical ones
  when all null hypotheses are true, so that you can have a
  general idea about whether there are lots of significant
  tests. It's analogous to the usual histogram of regular
  pvalues.
}
\details{
  If the plot shows a lack of big pvalues but overabundance
  of small pvalues, one expect there to be some
  siginificant tests.
}
\examples{
\dontrun{
   #load data
   data(leukemia)
   spt.mean(leukemia,5,10,1000) -> spt.mean.out
   plot.spt(spt.mean.out,col="red")
}
}
\author{
  Dan Nettleton and Chuanlong Du.
}

