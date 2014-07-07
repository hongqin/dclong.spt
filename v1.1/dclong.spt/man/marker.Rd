\docType{data}
\name{marker}
\alias{marker}
\title{Barley Marker}
\format{The dataset "barley" is a matrix with 22841 rows and 75 columns.}
\source{
  \url{http://www.plexdb.org/modules/PD_browse/experiment_browser.php}.
}
\description{
  Biologists genetically mutated/changed the genotypes of
  barley. They could not change everywhere, so they changed
  378 positions on the chromosome of barley. In the map,
  "A" and "B" are two types (sort of open and close).
  Because they know where the mutations are, they called
  them "markers" (so that if a barley with a certain
  genotype has a higher expression level, then you may
  infer and say, oh that may be caused by the 145th marker,
  etc.). The map has 7 chromosomes of barley, 1H, 2H, ...,
  7H. These numbers are locations of markers on the
  chromosomes, like coordinates. There are some missing
  values in the original map, a naive method was used to
  interpolate the missing values and produced this dataset
  "barley". Note that this data set is used for
  illustration in Bancoft et al. (2012).
}
\references{
  Tim Bancroft, Chuanlong Du and Dan Nettleton (2012)
  Estimation of False Discovery Rate Using Sequential
  Permutation P -Values.
}
\keyword{data}
\keyword{datasets}
\keyword{set}

