\name{spt.mean}
\alias{spt.mean}
\title{Sequential Permutation Test}
\usage{
  spt.mean(data, size1, h, n)
}
\arguments{
  \item{data}{a matrix with each row being a data set for a
  test (gene).}

  \item{size1}{the size of the first group.}

  \item{h}{the number of significant test statistics to hit
  before early termination.}

  \item{n}{the maximum number of permutations to use,
  including the observed one, i.e., at most n - 1
  permutations will be sampled.}
}
\value{
  an object of \code{\link[base]{class}} ``spt''. An object
  of class ``spt'' is a list containing at least the
  following components: \item{p}{pvalues of sequential
  permutation tests.} \item{h}{see description in the
  Arguments section.} \item{n}{see description in the
  Arguments section.} The object return by function
  \code{spt.corr} also contains the following component(s):
  \item{n.perms}{number of sequential permutations sampled
  for each test.} Note that though \code{spt.mean} performs
  sequential permutation tests, you can get results for
  regular permutation tests by setting \code{h >= n}.
}
\description{
  Performs (multiple) sequential permutation tests for no
  difference between two treatment groups.
}
\details{
  For each data set (a row in \code{data}) containing two
  groups, a sequential permutation test is done for H0:
  there's difference between the two groups VS no
  difference. The absolute difference between the means of
  the two groups is the test statistic used.
}
\examples{
\dontrun{
   #load data
   data(leukemia)
   spt.mean(leukemia,5,10,1000) -> spt.mean.out
   head(spt.mean.out)
}
}
\seealso{
  \code{\link{spt.corr}} for sequential permutations test
  for no correlation between a response variable and a
  bunch of covariates.
}

