\name{read_code}
\alias{read_code}
\title{Read in C++ code}
\usage{
  read_code(title)
}
\arguments{
  \item{title}{title of the file containing (C++) code.}
}
\description{
  The function \code{read_code} read in C++ code in the
  "rcpp" folder under the home directory of this package.
}
\examples{
read_code('spt.cpp')
read_code('test_statistics.cpp')
}

