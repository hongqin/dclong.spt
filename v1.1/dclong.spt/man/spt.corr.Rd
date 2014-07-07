\name{spt.corr}
\alias{spt.corr}
\title{Sequential Permutation Test}
\usage{
  spt.corr(x, y, h, n)
}
\arguments{
  \item{x}{a matrix/vector with columns containing
  covariates.}

  \item{y}{a matrix/vector with columns containing response
  variables.}

  \item{h}{the number of significant test statistics to hit
  before early termination.}

  \item{n}{the maximum number of permutations to use,
  including the observed one, i.e., at most n - 1
  permutations will be sampled.}
}
\value{
  an object of \code{\link[base]{class}} ``spt''. An object
  of class ``spt'' is a list containing at least the
  following components: \item{p}{pvalue of sequential
  permutation test for each gene/test.} \item{h}{see
  description in the Arguments section.} \item{n}{see
  description in the Arguments section.} The object return
  by function \code{spt.corr} also contains the following
  component(s): \item{max.ac}{the maximum absolute
  correlation (which is the test statistic) for each gene.}
  \item{max.index}{the index of the covariate where the
  maximum absolute correlations is achieved for each
  gene/test.} Note that though \code{spt.corr} performs
  sequential permutation test, you can get results for
  regular permutation test by setting \code{h >= n}.
}
\description{
  Performs (multiple) sequential permutation tests for no
  correlations between response variables and covariates.
}
\details{
  For each response varaible \code{y0} (a column in
  \code{y}), a sequential permutation test is done for H0:
  there's correlation between \code{y0} and covariates
  (columns in \code{x}) VS Ha: no correlation. The maximum
  absolute correlation is the test statistic used.
}
\examples{
\dontrun{
#load data
data(marker)
data(barley)
#sequntial permutation test for no correlation between gene expression
#and the markers (it might take a while)
spt.corr(t(marker),t(barley),10,1000)-> spt.corr.out
head(spt.corr.out)
}
}
\seealso{
  \code{\link{spt.mean}} sequential permutation tests for
  no difference between two treatment groups.
}

