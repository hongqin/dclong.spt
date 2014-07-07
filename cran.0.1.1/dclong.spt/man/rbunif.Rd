\name{rbunif}
\alias{rbunif}
\title{Mixture of Beta and Uniform Distribution}
\usage{
  rbunif(n, alpha, beta, gamma)
}
\arguments{
  \item{n}{number of observations to generate.}

  \item{alpha}{the first parameter of beta distribution.}

  \item{beta}{the second parameter of beta distribution.}

  \item{gamma}{probability of a observation coming from
  uniform distribution.}
}
\description{
  Generate random observations from a mixture of beta and
  uniform distribution. It's primarily used for
  demostrating use of other functions in this package.
}
\details{
  The mixture distribution is gamma*U(0,1) +
  (1-gamma)*Beta(alpha,beta).
}
\examples{
rbunif(100, alpha=1, beta=29, gamma=0.7)
}

