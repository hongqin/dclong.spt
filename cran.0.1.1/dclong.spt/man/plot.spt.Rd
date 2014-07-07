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
# ---------- a real example in Bancroft, Du and Nettleton (2013) --------------
# read in C++ source code
tsf_code = read_code("test_statistics.cpp")
# genrate function for sequential permutation test
# for no location difference between two samples
if(!exists("cxxspt.mean")){
    cxxwrapper(tsf_code, "abs_mean_diff0", "abs_mean_diff", c("n1", "", "int")) -> cxxspt.mean
}
# download data
if(!file.exists('spt_data.rda')){
    download.file('http://dclong.github.io/media/spt/spt_data.rda', 'spt_data.rda')
}
load('spt_data.rda')
cxxspt.mean(leukemia, 10, 1000, n1=5) -> cxxspt.mean.out
plot(cxxspt.mean.out, col="red")
}
\author{
  Dan Nettleton and Chuanlong Du.
}

