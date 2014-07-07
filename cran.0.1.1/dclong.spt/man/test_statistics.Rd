\name{abs_mean_diff0}
\alias{abs_mean_diff}
\alias{abs_mean_diff0}
\alias{max_abs_corr}
\alias{max_abs_corr0}
\title{User-defined Test Statistic Functions}
\usage{
  abs_mean_diff0(data, n1)

  abs_mean_diff(data, n1)

  max_abs_corr0(y, X)

  max_abs_corr(y, X)
}
\arguments{
  \item{data}{the data corresponding to a gene/test. In
  this case, it's a vector containing observations from two
  samples. The first \code{n1} observations belong to the
  first sample, and other observations correspond to the
  second sample.}

  \item{n1}{size of the first sample.}

  \item{y}{the data corresponding to a gene/test, which is
  a vector in this case. It must be centered and
  standardized.}

  \item{X}{a matrix with columns being covariates. Each of
  the columns/covariates must be centered and
  standardized.}
}
\description{
  The functions end with "0" calculate the observed test
  statistics. The corresponding versions that doesn't end
  with "0" permutes data and then recalculate the test
  statistics.
}

