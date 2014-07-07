\name{cxxwrapper}
\alias{cxxwrapper}
\title{C++ Code for Sequential Permutation Test}
\usage{
  cxxwrapper(tsf_code, tsf0, tsf, epList)
}
\arguments{
  \item{tsf_code}{the C++ code for user-defined test
  statistic functions. Similar to the R code, you have to
  write 2 C++ functions. One of them calculate the observed
  test statistic (based on original data), the other
  permutes data and calculates the test statistic again.}

  \item{tsf0}{the name of the C++ function which calculates
  the observed test statistic (based on original data).
  This function must take a (first) argument of type
  \code{arma::rowvec}. Extra parameters passing to this
  function are allowed. See \code{epList} for more on
  passing extra parameters to \code{tsf0}. For convenience
  of user, "abs_mean_diff0" and "max_abs_corr0" are defined
  in "rcpp/test_statistics.cpp" under the home directory of
  this package.}

  \item{tsf}{the name of the C++ function which permutes
  data and then calculate the test statistic again.
  \code{tsf0} and \code{tsf} must be defined in the way
  such that bigger test statistics are more extreme. This
  function must take a (first) argument of type
  \code{arma::rowvec}. Extra parameters passing to this
  function are allowed. See \code{epList} for more on
  passing extra parameters to \code{tsf}. For convenience
  of user, "abs_mean_diff" and "max_abs_corr" are defined
  in "rcpp/test_statistics.cpp" under the home directory of
  this package. Using \code{abs_mean_diff0} and
  \code{abs_mean_diff}, one can do sequential permutation
  test for no loation difference between two samples. Using
  \code{max_abs_corr0} and \code{max_abs_corr}, one can do
  sequential permutation test for no correlation between a
  reponse variable and a bunch of covariates.}

  \item{epList}{a string matrix (or string vector if only 1
  extra parameter) containing information of extra
  parameters to be passed to the user-defined functions
  \code{tsf0} and \code{tsf}. Each row of \code{epList}
  containing information of an extra parameter. The first
  column of \code{epList} contains names of extra
  parameters. The names of extra parameters must be valid
  C++ and R variable names. The names cannot be
  \code{data}, \code{h} or \code{n} which has already been
  used in the underlying C++ code. These names cannot
  conflict with the names of the user-defined test
  statistic functions \code{tsf0} and \code{tsf}. The
  second column of \code{epList} contains the data
  structure of the extra parameters. Currently you can only
  use matirx, vector and scalor as extra parameters in the
  user-defined test statistic functions. A specification
  "matrix" is transformed to \code{arma::Mat}, and "vector"
  is transformed to \code{arma::Col}. All other
  specifications are transformed to a scalor (i.e., no data
  structure is used). The third column of \code{epList}
  contains the type of data in the data structure. It must
  be valid C++ primitive types (e.g., double, float, int
  and bool). The following are some examples of
  specification of extra parameters. \code{c("x", "matrix",
  "double")}: 1 extra parameter \code{x} with type
  \code{arma::Mat<double>}. \code{c("x", "vector",
  "double")}: 1 extra parameter \code{x} with type
  \code{arma::Col<double>}. \code{c("x", "", "int")}: 1
  extra parameter \code{x} with type \code{int}.}
}
\value{
  an R function which calls C++ code for doing sequential
  permutation test. The parameters of the returned R
  function has the same meaning with parameters in the R
  function \code{\link{spt}}.
}
\description{
  The function \code{cxxwrapper} allows user to run C++
  code for sequential permutation test easily. Similar to
  what you can do with the R code in this package, you can
  write C++ functions to permute data and calculate test
  statistics, and then use the function \code{cxxwrapper}
  to compile your C++ code and generate an R function
  calling the compiled code. To solve the problem of
  passing data between R and C++, the underlying C++ code
  limits the data structure that you can use.
  \code{arma::mat} (a matrix with double elements) is
  currently the only supported data structure for storing
  data corresponding to all genes. Each row of the matrix
  contains data corresponding to one gene.
}
\examples{
# read source code
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
# generate function for sequential permutation test
# for no correlation between a response variable and a bunch of covariates
if(!exists("cxxspt.corr")){
    cxxwrapper(tsf_code, "max_abs_corr0", "max_abs_corr", c("X", "matrix", "double")) -> cxxspt.corr
}
# center and standardize columns of marker
X = scale(t(marker))
# center and standardized rows of barley
Y = t(scale(t(barley)))
# for the purpose of illustration, only test the first 100 genes since it's time consuming
cxxspt.corr(Y[1:100,], 10, 1000, X=X) -> cxxspt.corr.out
}
\seealso{
  \code{\link{spt}} which is a pure R version function for
  doing sequential permutation test.
}

